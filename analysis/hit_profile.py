#!/usr/bin/env python3
"""Profile a hit list: chemical space, internal diversity, shared binding mode.

Two questions you want answered about a ranking before ordering anything:

  Are these one scaffold or many?      -> Morgan fingerprints, PCA, and the
                                          nearest-neighbour Tanimoto of each
                                          hit within the list
  Do they bind the same way?           -> which pocket residues each hit
                                          contacts, from PLIP reports

A hit list that collapses to one cluster with nearest-neighbour similarities
near 1.0 is usually the scoring function rewarding a single motif, not a
discovery. A hit list where every compound touches a different set of residues
usually means the poses are unreliable.

The interaction half is optional and needs PLIP reports, e.g. the ones
constraints/07_generate_constraint_yamls.py leaves in its --plip-output folder.
"""
import argparse
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

SMILES_COLUMN_CANDIDATES = ("smiles", "SMILES", "canonical_smiles", "Smiles")
ID_COLUMN_CANDIDATES = ("molport_number", "compound_id", "id", "ID", "chembl_id")
SCORE_COLUMN_CANDIDATES = ("combined_score", "boltz_pred_value_mean", "pred_score",
                           "prediction_score", "score")

INTERACTION_GROUPS = [
    ("hydrogen_bonds", "hydrogen_bond"),
    ("hydrophobic_interactions", "hydrophobic_interaction"),
    ("pi_stacks", "pi_stack"),
    ("pi_cation_interactions", "pi_cation_interaction"),
    ("salt_bridges", "salt_bridge"),
    ("water_bridges", "water_bridge"),
    ("halogen_bonds", "halogen_bond"),
]

RESIDUE_RE = re.compile(r"^\s*(\d+)\s*([A-Za-z]?)\s*$")


def pick_column(header, requested, candidates, what, required=True):
    if requested:
        if requested not in header:
            raise SystemExit(
                f"ERROR: {what} column '{requested}' not in the hit list.\n"
                f"       Available columns: {', '.join(header)}"
            )
        return requested
    for cand in candidates:
        if cand in header:
            return cand
    if required:
        raise SystemExit(
            f"ERROR: no {what} column found. Looked for: {', '.join(candidates)}.\n"
            f"       Available columns: {', '.join(header)}\n"
            f"       Name it explicitly with --{what}-col."
        )
    return None


def fingerprint_matrix(smiles_list):
    """Morgan fingerprints, radius 2, 2048 bits, as a dense uint8 matrix."""
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator

    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    rows, keep, failed = [], [], []
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        if mol is None:
            failed.append(i)
            continue
        arr = np.zeros((2048,), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(gen.GetFingerprint(mol), arr)
        rows.append(arr)
        keep.append(i)
    if not rows:
        raise SystemExit("ERROR: RDKit could not parse a single SMILES in the hit list.")
    return np.vstack(rows), keep, failed


def nearest_neighbour_tanimoto(matrix):
    """For each compound, the highest Tanimoto to any *other* compound.

    Computed on the dense bit matrix: |A & B| / (|A| + |B| - |A & B|).
    """
    m = matrix.astype(np.float32)
    inter = m @ m.T
    counts = m.sum(axis=1)
    union = counts[:, None] + counts[None, :] - inter
    with np.errstate(divide="ignore", invalid="ignore"):
        tan = np.where(union > 0, inter / union, 0.0)
    np.fill_diagonal(tan, -1.0)          # exclude self
    return tan.max(axis=1)


def parse_plip_report(xml_path):
    """Interaction-type counts and contact residues from one PLIP report."""
    root = ET.parse(xml_path).getroot()

    counts = {}
    for group, item in INTERACTION_GROUPS:
        node = root.find(f".//interactions/{group}")
        counts[group] = len(node.findall(item)) if node is not None else 0

    residues = []
    for node in root.findall(".//bs_residue"):
        if (node.attrib.get("contact") or "").lower() != "true":
            continue
        match = RESIDUE_RE.match(node.text or "")
        if not match:
            continue
        resnum, chain = match.group(1), match.group(2) or ""
        aa = node.attrib.get("aa") or ""
        residues.append(f"{aa}{resnum}{chain}")
    counts["contact_residues"] = sorted(set(residues))
    return counts


def find_report(plip_dir, compound_id):
    """Locate the PLIP report belonging to a compound.

    Folder layouts differ between runs, so both '<dir>/<id>/report.xml' and
    any '*report*.xml' under a folder whose name contains the ID are accepted.
    """
    direct = os.path.join(plip_dir, str(compound_id), "report.xml")
    if os.path.isfile(direct):
        return direct
    for entry in os.scandir(plip_dir):
        if entry.is_dir() and str(compound_id) in entry.name:
            for name in sorted(os.listdir(entry.path)):
                if name.endswith(".xml") and "report" in name:
                    return os.path.join(entry.path, name)
    return None


def plot_chemical_space(df, variance, out_path, score_label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    sc = ax1.scatter(df["pc1"], df["pc2"], c=df["_score"], cmap="viridis_r",
                     s=45, alpha=0.85, edgecolors="black", linewidths=0.3)
    ax1.set_xlabel(f"PC1 ({variance[0] * 100:.1f} % variance)")
    ax1.set_ylabel(f"PC2 ({variance[1] * 100:.1f} % variance)")
    ax1.set_title("Chemical space of the hit list\n(Morgan fingerprints, radius 2)")
    fig.colorbar(sc, ax=ax1, label=score_label)

    ax2.hist(df["nn_tanimoto"], bins=25, range=(0, 1), color="steelblue",
             alpha=0.8, edgecolor="black")
    median = df["nn_tanimoto"].median()
    ax2.axvline(median, color="red", linestyle="--",
                label=f"median {median:.2f}")
    ax2.set_xlabel("Tanimoto to the most similar other hit")
    ax2.set_ylabel("Count")
    ax2.set_title("Internal diversity\n(higher = the list repeats itself)")
    ax2.legend()

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_interaction_heatmap(matrix_df, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    height = max(4.0, 0.28 * len(matrix_df))
    width = max(6.0, 0.42 * len(matrix_df.columns))
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(matrix_df, cmap="Blues", cbar=False, linewidths=0.4,
                linecolor="white", ax=ax, vmin=0, vmax=1)
    ax.set_xlabel("Pocket residue")
    ax.set_ylabel("")
    ax.set_title("Contacts per compound (PLIP)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(
        description="Chemical space, internal diversity and shared contacts of a hit list",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  python analysis/hit_profile.py --hits results/ca2/ranking/ligands_ranked_top100.csv \\
      --out-dir results/profile

  # with the PLIP reports from the constraint pipeline
  python analysis/hit_profile.py --hits my_hits.csv --out-dir results/profile \\
      --plip-dir results/constraints/ca2/plip

Reads a hit list, not a screening library: a few dozen to a few hundred rows.
The PCA is fitted on the hit list itself, so the axes mean nothing across runs.
""",
    )
    ap.add_argument("--hits", required=True, help="Ranked hit list (CSV)")
    ap.add_argument("--out-dir", required=True, help="Output directory")
    ap.add_argument("--plip-dir", help="Folder of PLIP reports, one per compound (optional)")
    ap.add_argument("--smiles-col", help="SMILES column (default: auto-detect)")
    ap.add_argument("--id-col", help="Compound ID column (default: auto-detect)")
    ap.add_argument("--score-col", help="Score column for colouring (default: auto-detect)")
    ap.add_argument("--top-n", type=int, default=None,
                    help="Use only the first N rows of the hit list")
    ap.add_argument("--min-residue-frequency", type=int, default=1,
                    help="Drop residues contacted by fewer than N compounds (default: 1)")
    args = ap.parse_args()

    if not os.path.isfile(args.hits):
        raise SystemExit(f"ERROR: hit list not found: {args.hits}")

    df = pd.read_csv(args.hits)
    if args.top_n:
        df = df.head(args.top_n)
    header = list(df.columns)
    smi_col = pick_column(header, args.smiles_col, SMILES_COLUMN_CANDIDATES, "smiles")
    id_col = pick_column(header, args.id_col, ID_COLUMN_CANDIDATES, "id")
    score_col = pick_column(header, args.score_col, SCORE_COLUMN_CANDIDATES,
                            "score", required=False)

    print(f"Loaded {len(df)} hit(s) from {args.hits}")
    print(f"  smiles: '{smi_col}'  id: '{id_col}'  score: '{score_col or 'none, using rank'}'")

    os.makedirs(args.out_dir, exist_ok=True)

    # --- Chemical space and diversity ---
    matrix, keep, failed = fingerprint_matrix(df[smi_col].tolist())
    if failed:
        print(f"  {len(failed)} SMILES could not be parsed and were dropped")
    kept = df.iloc[keep].reset_index(drop=True)

    from sklearn.decomposition import PCA
    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(matrix.astype(np.float32))

    out = pd.DataFrame({
        id_col: kept[id_col].values,
        "pc1": coords[:, 0],
        "pc2": coords[:, 1],
        "nn_tanimoto": nearest_neighbour_tanimoto(matrix),
    })
    if score_col:
        out["_score"] = kept[score_col].values
        score_label = score_col
    else:
        out["_score"] = np.arange(1, len(out) + 1)
        score_label = "rank in the input"

    space_png = os.path.join(args.out_dir, "chemical_space.png")
    plot_chemical_space(out, pca.explained_variance_ratio_, space_png, score_label)

    csv_out = out.rename(columns={"_score": score_col or "rank"})
    csv_out.to_csv(os.path.join(args.out_dir, "chemical_space.csv"), index=False)

    nn = out["nn_tanimoto"]
    print(f"\nInternal diversity over {len(out)} hits")
    print(f"  nearest-neighbour Tanimoto: median {nn.median():.3f}, "
          f"mean {nn.mean():.3f}, max {nn.max():.3f}")
    print(f"  {int((nn >= 0.8).sum())} hit(s) have a near-duplicate in the list (Tc >= 0.8)")
    print(f"  PC1 + PC2 explain {pca.explained_variance_ratio_[:2].sum() * 100:.1f} % of the variance")
    print(f"Wrote {space_png}")

    # --- Interactions ---
    if not args.plip_dir:
        print("\nNo --plip-dir given, skipping the interaction profile.")
        return 0

    if not os.path.isdir(args.plip_dir):
        raise SystemExit(f"ERROR: PLIP folder not found: {args.plip_dir}")

    records, residue_sets, missing = [], {}, []
    for compound_id in kept[id_col]:
        report = find_report(args.plip_dir, compound_id)
        if report is None:
            missing.append(compound_id)
            continue
        try:
            parsed = parse_plip_report(report)
        except ET.ParseError as exc:
            print(f"  WARNING: unreadable report for {compound_id}: {exc}", file=sys.stderr)
            continue
        residue_sets[compound_id] = parsed.pop("contact_residues")
        parsed[id_col] = compound_id
        records.append(parsed)

    if not records:
        raise SystemExit(
            f"ERROR: no PLIP report matched a compound ID under {args.plip_dir}.\n"
            f"       Expected '<id>/report.xml' or a folder whose name contains the ID."
        )
    if missing:
        shown = ", ".join(str(m) for m in missing[:3])
        print(f"\n{len(missing)} compound(s) had no PLIP report, e.g. {shown}")

    summary = pd.DataFrame(records)
    summary["contact_residues"] = summary[id_col].map(
        lambda c: ";".join(residue_sets.get(c, []))
    )
    cols = [id_col] + [c for c in summary.columns if c != id_col]
    summary = summary[cols]
    summary.to_csv(os.path.join(args.out_dir, "interaction_summary.csv"), index=False)

    all_residues = sorted(
        {r for rs in residue_sets.values() for r in rs},
        key=lambda r: int(re.sub(r"\D", "", r) or 0),
    )
    matrix_df = pd.DataFrame(
        [[1 if r in residue_sets[c] else 0 for r in all_residues] for c in residue_sets],
        index=list(residue_sets), columns=all_residues,
    )
    if args.min_residue_frequency > 1:
        keep_cols = matrix_df.columns[matrix_df.sum() >= args.min_residue_frequency]
        dropped = len(matrix_df.columns) - len(keep_cols)
        matrix_df = matrix_df[keep_cols]
        print(f"  dropped {dropped} residue(s) contacted by fewer than "
              f"{args.min_residue_frequency} compounds")

    heat_png = os.path.join(args.out_dir, "interaction_heatmap.png")
    plot_interaction_heatmap(matrix_df, heat_png)

    freq = matrix_df.sum().sort_values(ascending=False)
    print(f"\nContacts over {len(matrix_df)} compound(s) with a report")
    print(f"  {len(all_residues)} distinct residues contacted")
    print("  most frequent:")
    for residue, n in freq.head(5).items():
        print(f"    {residue:<10} {n:>4} / {len(matrix_df)}  ({n / len(matrix_df) * 100:.0f} %)")
    print(f"Wrote {heat_png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
