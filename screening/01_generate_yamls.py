#!/usr/bin/env python3
import argparse
import csv
import os
import sys


def read_fasta_sequence(fasta_path: str) -> str:
    if not os.path.isfile(fasta_path):
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")
    seq = []
    with open(fasta_path, "r", encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith(">"):
                continue
            seq.append(s)
    if not seq:
        raise ValueError(f"No sequence lines in FASTA: {fasta_path}")
    return "".join(seq).upper()


def yaml_quote(s: str) -> str:
    # Minimal YAML-safe double-quote escaping
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def ensure_outdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


MAX_CHAIN_ID_LEN = 5


def shorten_chain_id(raw_id: str) -> str:
    """Boltz chain names are limited (Lightning stores them in fixed-length keys).
    To avoid KeyError (e.g. Z17024 -> Z1702), we truncate after MAX_CHAIN_ID_LEN characters.
    """
    if len(raw_id) <= MAX_CHAIN_ID_LEN:
        return raw_id
    return raw_id[:MAX_CHAIN_ID_LEN]


def make_yaml_text(sequence: str, ligand_id: str, smiles: str, chain_id: str, msa_path: str | None) -> str:
    # Build YAML string manually to avoid external deps
    lines = []
    lines.append("version: 1")
    lines.append("")
    lines.append("sequences:")
    lines.append("  - protein:")
    lines.append(f"      id: [{chain_id}]")
    lines.append(f"      sequence: {sequence}")
    if msa_path:
        lines.append(f"      msa: {msa_path}")
    lines.append("")
    lines.append("  - ligand:")
    lines.append(f"      id: [{yaml_quote(ligand_id)}]")
    lines.append(f"      smiles: {yaml_quote(smiles)}")
    lines.append("")
    lines.append("properties:")
    lines.append("  - affinity:")
    lines.append(f"      binder: {yaml_quote(ligand_id)}")
    lines.append("")
    return "\n".join(lines)


def read_ism_file(ism_path: str) -> list:
    """Read .ism file and return list of (smiles, id, chembl_id) tuples"""
    if not os.path.isfile(ism_path):
        raise FileNotFoundError(f"ISM file not found: {ism_path}")

    entries = []
    with open(ism_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            # Format: SMILES ID CHEMBL_ID (space-separated)
            parts = line.split()
            if len(parts) >= 2:
                smiles = parts[0]
                mol_id = parts[1]
                chembl_id = parts[2] if len(parts) >= 3 else ""
                entries.append((smiles, mol_id, chembl_id))

    return entries


SMILES_COLUMN_CANDIDATES = ("smiles", "SMILES", "canonical_smiles", "Smiles")
ID_COLUMN_CANDIDATES = ("molport_number", "compound_id", "id", "ID", "chembl_id", "name")


def pick_column(header: list, requested: str | None, candidates: tuple, what: str) -> str:
    """Return the column to use, either the requested one or the first known name."""
    if requested:
        if requested not in header:
            raise SystemExit(
                f"ERROR: {what} column '{requested}' not in CSV.\n"
                f"       Available columns: {', '.join(header)}"
            )
        return requested
    for cand in candidates:
        if cand in header:
            return cand
    raise SystemExit(
        f"ERROR: no {what} column found in CSV. Looked for: {', '.join(candidates)}.\n"
        f"       Available columns: {', '.join(header)}\n"
        f"       Name it explicitly with --{what}-column."
    )


def read_smiles_csv(csv_path: str, smiles_column: str | None, id_column: str | None) -> list:
    """Read a compound table and return (smiles, id, "") tuples.

    The screening input is a CSV, not the .ism format inherited from the
    benchmark pipeline. Rows without a SMILES string are skipped and counted.
    """
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    entries = []
    skipped = 0
    with open(csv_path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise SystemExit(f"ERROR: CSV has no header row: {csv_path}")
        header = list(reader.fieldnames)
        smi_col = pick_column(header, smiles_column, SMILES_COLUMN_CANDIDATES, "smiles")
        id_col = pick_column(header, id_column, ID_COLUMN_CANDIDATES, "id")
        print(f"Reading {csv_path} (smiles: '{smi_col}', id: '{id_col}')")

        for row in reader:
            smiles = (row.get(smi_col) or "").strip()
            mol_id = (row.get(id_col) or "").strip()
            if not smiles or not mol_id:
                skipped += 1
                continue
            entries.append((smiles, mol_id, ""))

    if skipped:
        print(f"  skipped {skipped} row(s) without a SMILES string or an ID")
    if not entries:
        raise SystemExit(f"ERROR: no usable rows in {csv_path}")
    return entries


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Boltz-2 YAML generator (one YAML per ligand, affinity property enabled)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  # screening: one YAML per compound from a CSV
  python 01_generate_yamls.py --smiles-csv data/demo/demo_100.csv \\
      --fasta data/fasta/P00918.fasta --msa-path data/msa/P00918.a3m \\
      --outdir output/ca2/yamls

  # benchmark style: actives and decoys from .ism files
  python 01_generate_yamls.py --actives actives.ism --decoys decoys.ism \\
      --fasta data/fasta/P00918.fasta --outdir output/yamls

Ligand IDs inside the YAML are sequential (L1, L2, ... for actives and CSV
input, Z1, Z2, ... for decoys) and truncated to 5 characters, a Boltz chain-name
limit. The compound ID goes into the filename, which is what
03_make_boltz_predictions_summary.py matches on.
""",
    )
    ap.add_argument("--smiles-csv", help="Path to a compound table (CSV) with SMILES and ID columns")
    ap.add_argument("--smiles-column", help="Name of the SMILES column (default: auto-detect)")
    ap.add_argument("--id-column", help="Name of the compound ID column (default: auto-detect)")
    ap.add_argument("--actives", help="Path to actives .ism file")
    ap.add_argument("--decoys", help="Path to decoys .ism file")
    ap.add_argument("--fasta", required=True, help="Path to FASTA with protein sequence")
    ap.add_argument("--outdir", required=True, help="Output directory for YAML files")
    ap.add_argument("--limit", type=int, default=None, help="Max number of YAMLs to create")
    ap.add_argument("--chain-id", default="A", help="Protein chain id (default: A)")
    ap.add_argument("--prefix", default="boltz", help="Filename prefix (default: boltz)")
    ap.add_argument("--msa-path", help="Optional path to MSA file (e.g. data/msa/P00918.a3m)")
    args = ap.parse_args()

    if not args.actives and not args.decoys and not args.smiles_csv:
        ap.error("One of --smiles-csv, --actives or --decoys must be provided")

    ensure_outdir(args.outdir)
    sequence = read_fasta_sequence(args.fasta)

    n = 0
    active_counter = 0
    decoy_counter = 0

    # Process CSV input and actives through the same branch: both are "L" ligands
    if args.smiles_csv:
        entries = read_smiles_csv(args.smiles_csv, args.smiles_column, args.id_column)
    elif args.actives:
        entries = read_ism_file(args.actives)
    else:
        entries = []

    if entries:
        for smiles, mol_id, chembl_id in entries:
            if not smiles:
                continue
            
            active_counter += 1
            raw_ligand_id = f"L{active_counter}"
            ligand_id = shorten_chain_id(raw_ligand_id)
            
            yaml_text = make_yaml_text(
                sequence=sequence,
                ligand_id=ligand_id,
                smiles=smiles,
                chain_id=args.chain_id,
                msa_path=args.msa_path,
            )
            
            # Use CHEMBL ID for actives, fallback to mol_id or ligand_id
            out_stub = chembl_id if chembl_id else (mol_id if mol_id else raw_ligand_id)
            out_name = f"{args.prefix}_{out_stub}.yaml"
            out_path = os.path.join(args.outdir, out_name)
            with open(out_path, "w", encoding="utf-8") as out:
                out.write(yaml_text)
            
            n += 1
            if args.limit and n >= args.limit:
                break

    # Process decoys
    if args.decoys and (not args.limit or n < args.limit):
        entries = read_ism_file(args.decoys)
        for smiles, mol_id, chembl_id in entries:
            if not smiles:
                continue
            
            decoy_counter += 1
            raw_ligand_id = f"Z{decoy_counter}"
            ligand_id = shorten_chain_id(raw_ligand_id)
            
            yaml_text = make_yaml_text(
                sequence=sequence,
                ligand_id=ligand_id,
                smiles=smiles,
                chain_id=args.chain_id,
                msa_path=args.msa_path,
            )
            
            # Use original ID in filename for traceability
            out_stub = mol_id if mol_id else raw_ligand_id
            out_name = f"{args.prefix}_{out_stub}.yaml"
            out_path = os.path.join(args.outdir, out_name)
            with open(out_path, "w", encoding="utf-8") as out:
                out.write(yaml_text)
            
            n += 1
            if args.limit and n >= args.limit:
                break

    if n == 0:
        print(
            "ERROR: no YAML files written — every input row was empty.",
            file=sys.stderr,
        )
        return 1

    if args.smiles_csv:
        print(f"Wrote {n} YAML file(s) to {args.outdir}")
    else:
        print(
            f"Wrote {n} YAML file(s) to {args.outdir} "
            f"({active_counter} actives, {decoy_counter} decoys)"
        )

    if not args.msa_path:
        print(
            "NOTE: no --msa-path given. The YAMLs carry no MSA and Boltz-2 will "
            "need --use_msa_server, or run 02_update_msa.py before inference.",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())

