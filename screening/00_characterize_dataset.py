#!/usr/bin/env python3
"""Annotate a compound table with the molecular properties the pipeline filters on.

Step 1 of the screening workflow. Reads a vendor export (CSV or XLSX) holding
SMILES and compound IDs, validates every structure with RDKit, and writes a CSV
carrying molecular weight, heavy-atom count, logP, hydrogen-bond donors and
acceptors, rotatable bonds and a PAINS flag.

04_analyze_ligand_predictions.py filters on these columns; without them the
size filters that keep ligands inside Boltz-2's supported range silently do
nothing.
"""
import argparse
import os
import sys

import pandas as pd

SMILES_COLUMN_CANDIDATES = ("smiles", "SMILES", "Smiles", "canonical_smiles")
ID_COLUMN_CANDIDATES = (
    "molport_number", "Molport ID", "MolPort ID", "Molport Id",
    "compound_id", "id", "ID", "chembl_id",
)


def calc_properties(smiles, catalog):
    """Molecular properties for one SMILES, or None if RDKit rejects it."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, Crippen

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return {
        "mw": Descriptors.MolWt(mol),
        "heavy_atoms": Lipinski.HeavyAtomCount(mol),
        "num_atoms": mol.GetNumAtoms(),
        "logp": Crippen.MolLogP(mol),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "pains": catalog.HasMatch(mol),
    }


def pick_column(header, requested, candidates, what):
    if requested:
        if requested not in header:
            raise SystemExit(
                f"ERROR: {what} column '{requested}' not in the input.\n"
                f"       Available columns: {', '.join(header)}"
            )
        return requested
    for cand in candidates:
        if cand in header:
            return cand
    raise SystemExit(
        f"ERROR: no {what} column found. Looked for: {', '.join(candidates)}.\n"
        f"       Available columns: {', '.join(header)}\n"
        f"       Name it explicitly with --{what}-column."
    )


def read_table(path):
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xls", ".xlsm"):
        return pd.read_excel(path)
    if ext in (".tsv", ".tab"):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def plot_distributions(df, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle("Property distributions", fontsize=16)

    panels = [
        (axes[0, 0], "mw", "Molecular weight (Da)", "steelblue", True),
        (axes[0, 1], "logp", "logP", "green", True),
        (axes[0, 2], "hbd", "Hydrogen-bond donors", "orange", False),
        (axes[1, 0], "hba", "Hydrogen-bond acceptors", "purple", False),
        (axes[1, 1], "rotatable_bonds", "Rotatable bonds", "brown", False),
    ]
    for ax, col, label, colour, continuous in panels:
        if continuous:
            ax.hist(df[col], bins=50, color=colour, alpha=0.7, edgecolor="black")
            ax.axvline(df[col].median(), color="red", linestyle="--",
                       label=f"Median: {df[col].median():.2f}")
            ax.legend()
        else:
            ax.hist(df[col], bins=range(0, int(df[col].max()) + 2),
                    color=colour, alpha=0.7, edgecolor="black")
        ax.set_xlabel(label)
        ax.set_ylabel("Count")

    axes[1, 2].scatter(df["mw"], df["logp"], alpha=0.3, s=10)
    axes[1, 2].set_xlabel("Molecular weight (Da)")
    axes[1, 2].set_ylabel("logP")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(
        description="Annotate a compound table with molecular properties and a PAINS flag",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  python 00_characterize_dataset.py --input library.xlsx --out-dir data/
  python 00_characterize_dataset.py --input library.csv --out-dir data/ --plot

Accepts .csv, .tsv and .xlsx. SMILES and ID columns are auto-detected; vendor
exports name them differently ('SMILES' vs 'smiles', 'Molport ID' vs
'molport_number'), so override with --smiles-column / --id-column if needed.
""",
    )
    ap.add_argument("--input", required=True, help="Compound table (CSV, TSV or XLSX)")
    ap.add_argument("--out-dir", required=True, help="Output directory")
    ap.add_argument("--out-name", default=None,
                    help="Output filename (default: <input stem>_properties.csv)")
    ap.add_argument("--smiles-column", help="SMILES column (default: auto-detect)")
    ap.add_argument("--id-column", help="Compound ID column (default: auto-detect)")
    ap.add_argument("--id-out-column", default="molport_number",
                    help="Name of the ID column in the output (default: molport_number)")
    ap.add_argument("--plot", action="store_true",
                    help="Also write a property-distribution figure")
    args = ap.parse_args()

    try:
        from rdkit.Chem import FilterCatalog
        from rdkit import RDLogger
    except ImportError:
        raise SystemExit("ERROR: RDKit is not installed. pip install rdkit")
    RDLogger.DisableLog("rdApp.*")

    if not os.path.isfile(args.input):
        raise SystemExit(f"ERROR: input not found: {args.input}")

    df = read_table(args.input)
    header = list(df.columns)
    smi_col = pick_column(header, args.smiles_column, SMILES_COLUMN_CANDIDATES, "smiles")
    id_col = pick_column(header, args.id_column, ID_COLUMN_CANDIDATES, "id")
    print(f"Read {len(df)} row(s) from {args.input} (smiles: '{smi_col}', id: '{id_col}')")

    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)

    rows = []
    n_failed = 0
    for i, row in enumerate(df.itertuples(index=False), start=1):
        if i % 1000 == 0:
            print(f"  {i}/{len(df)}")
        smiles = getattr(row, smi_col.replace(" ", "_"), None) if hasattr(row, smi_col.replace(" ", "_")) else None
        if smiles is None:
            smiles = df.iloc[i - 1][smi_col]
        if not isinstance(smiles, str) or not smiles.strip():
            n_failed += 1
            continue
        props = calc_properties(smiles.strip(), catalog)
        if props is None:
            n_failed += 1
            continue
        props[args.id_out_column] = df.iloc[i - 1][id_col]
        props["smiles"] = smiles.strip()
        rows.append(props)

    if not rows:
        raise SystemExit(
            "ERROR: RDKit accepted none of the structures. "
            f"Check that '{smi_col}' really holds SMILES."
        )

    out = pd.DataFrame(rows)
    cols = [args.id_out_column, "smiles"] + [c for c in out.columns
                                             if c not in (args.id_out_column, "smiles")]
    out = out[cols]

    os.makedirs(args.out_dir, exist_ok=True)
    name = args.out_name or f"{os.path.splitext(os.path.basename(args.input))[0]}_properties.csv"
    out_csv = os.path.join(args.out_dir, name)
    out.to_csv(out_csv, index=False)

    print(f"\nValid: {len(out)}/{len(df)}" + (f"  ({n_failed} rejected)" if n_failed else ""))
    print(f"PAINS flagged: {int(out['pains'].sum())} "
          f"({out['pains'].sum() / len(out) * 100:.1f}%)")
    for col, unit in [("mw", " Da"), ("logp", ""), ("heavy_atoms", ""),
                      ("hbd", ""), ("hba", ""), ("rotatable_bonds", "")]:
        print(f"  {col:<16} median {out[col].median():>7.1f}{unit}"
              f"   range {out[col].min():.1f} - {out[col].max():.1f}")
    print(f"\nWrote {out_csv}")

    if args.plot:
        fig_path = os.path.splitext(out_csv)[0] + "_distributions.png"
        plot_distributions(out, fig_path)
        print(f"Wrote {fig_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
