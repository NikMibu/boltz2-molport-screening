#!/usr/bin/env python3
"""Collect Boltz-2 affinity predictions and join them onto a compound table.

Reads a Boltz-2 output tree, extracts the affinity metrics from every
``affinity_*.json``, and writes the input CSV back out with two extra columns:
``boltz_prob_mean`` and ``boltz_pred_value_mean``.

Prediction folders are matched to table rows by compound ID. Boltz names each
folder after the input YAML, so ``boltz_<ID>.yaml`` becomes a folder
``boltz_<ID>``; the ID is recovered from the folder name and looked up in the
column given by ``--id-column``.
"""
import argparse
import os
import json
import sys
import numpy as np
import pandas as pd


def safe_mean_from_prefix(d: dict, prefix: str) -> float:
    vals = [float(v) for k, v in d.items() if k.startswith(prefix) and isinstance(v, (int, float))]
    return float(np.mean(vals)) if vals else float("nan")


def normalize_id(compound_id: str) -> str:
    """Normalize a compound ID for matching (uppercase, strip whitespace)."""
    return str(compound_id).strip().upper()


def find_predictions_dir(root: str) -> str:
    """Find the predictions directory (may be nested)."""
    # Direct predictions folder
    direct = os.path.join(root, "predictions")
    if os.path.isdir(direct):
        return direct

    # Nested: boltz_output/boltz_results_*/predictions
    for entry in os.scandir(root):
        if entry.is_dir():
            nested = os.path.join(entry.path, "predictions")
            if os.path.isdir(nested):
                return nested

    return root


def id_from_folder(name: str, known_ids: set) -> str | None:
    """Recover the compound ID from a Boltz prediction folder name.

    Boltz derives the folder name from the input YAML stem, so a YAML written as
    ``<prefix>_<ID>.yaml`` yields ``<prefix>_<ID>``. Prefixes vary between runs
    (``boltz_``, ``exp_1_1_``, none at all) and IDs may themselves contain
    underscores, so candidates are tested against the IDs actually present in
    the table rather than guessed from the shape of the name.
    """
    candidates = [name]
    parts = name.split("_")
    # Every suffix of the split, longest first: recovers IDs containing "_"
    for i in range(1, len(parts)):
        candidates.append("_".join(parts[i:]))
    # Individual segments, for IDs that carry no underscore
    candidates.extend(parts)

    for cand in candidates:
        key = normalize_id(cand)
        if key in known_ids:
            return key
    return None


def main():
    ap = argparse.ArgumentParser(
        description="Join Boltz-2 affinity predictions onto a compound table",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
example:
  python 03_make_boltz_predictions_summary.py \\
      --results-dir output/ca2/boltz_output \\
      --input data/my_compounds.csv \\
      --id-column molport_number \\
      --out-dir results/data
""",
    )
    ap.add_argument("--results-dir", required=True, help="Path to Boltz output folder")
    ap.add_argument(
        "--input",
        "--molport-properties",
        dest="input",
        required=True,
        help="Compound table (CSV) to annotate; must contain the ID column",
    )
    ap.add_argument(
        "--id-column",
        default="molport_number",
        help="Column holding the compound IDs used in the YAML filenames (default: molport_number)",
    )
    ap.add_argument("--out-dir", required=True, help="Output directory for the annotated CSV")
    ap.add_argument(
        "--out-name",
        default=None,
        help="Output filename (default: <input stem>_with_predictions.csv)",
    )
    args = ap.parse_args()

    if not os.path.isfile(args.input):
        print(f"ERROR: compound table not found: {args.input}", file=sys.stderr)
        return 2

    print(f"Loading compound table from: {args.input}")
    df = pd.read_csv(args.input)

    if args.id_column not in df.columns:
        print(
            f"ERROR: column '{args.id_column}' not found in {args.input}.\n"
            f"       Available columns: {', '.join(df.columns)}\n"
            f"       Set the right one with --id-column.",
            file=sys.stderr,
        )
        return 2

    # Create normalized lookup key
    df["_lookup_key"] = df[args.id_column].apply(normalize_id)
    known_ids = set(df["_lookup_key"])

    # Scan Boltz predictions
    if not os.path.isdir(args.results_dir):
        print(f"ERROR: results directory not found: {args.results_dir}", file=sys.stderr)
        return 2

    pred_root = find_predictions_dir(args.results_dir)
    print(f"Scanning predictions in: {pred_root}")

    predictions = {}
    n_folders = 0
    n_unmatched = 0
    n_no_affinity = 0
    unmatched_examples = []

    for entry in os.scandir(pred_root):
        if not entry.is_dir():
            continue
        n_folders += 1

        lookup_key = id_from_folder(entry.name, known_ids)
        if lookup_key is None:
            n_unmatched += 1
            if len(unmatched_examples) < 3:
                unmatched_examples.append(entry.name)
            continue

        # Find affinity JSON
        aff_candidates = [p for p in os.listdir(entry.path) if p.startswith("affinity_") and p.endswith(".json")]
        if not aff_candidates:
            n_no_affinity += 1
            continue

        aff_path = os.path.join(entry.path, aff_candidates[0])
        try:
            with open(aff_path) as f:
                d = json.load(f)
        except Exception as exc:
            print(f"  WARNING: could not read {aff_path}: {exc}", file=sys.stderr)
            continue

        predictions[lookup_key] = {
            "boltz_prob_mean": safe_mean_from_prefix(d, "affinity_probability_binary"),
            "boltz_pred_value_mean": safe_mean_from_prefix(d, "affinity_pred_value"),
        }

    print(f"Scanned {n_folders} prediction folder(s), matched {len(predictions)}")
    if n_no_affinity:
        print(f"  {n_no_affinity} folder(s) carried no affinity_*.json — incomplete predictions")
    if n_unmatched:
        print(
            f"  {n_unmatched} folder(s) could not be matched to a row in "
            f"'{args.id_column}', e.g. {', '.join(unmatched_examples)}"
        )

    if not predictions:
        print(
            "ERROR: no prediction folder matched a compound ID.\n"
            f"       Folder names are compared against the '{args.id_column}' column; "
            "check --id-column and --results-dir.\n"
            "       Writing an all-empty prediction table would hide this, so nothing was written.",
            file=sys.stderr,
        )
        return 1

    # Add prediction columns
    df["boltz_prob_mean"] = df["_lookup_key"].map(lambda k: predictions.get(k, {}).get("boltz_prob_mean", np.nan))
    df["boltz_pred_value_mean"] = df["_lookup_key"].map(lambda k: predictions.get(k, {}).get("boltz_pred_value_mean", np.nan))

    # Drop temporary lookup key
    df = df.drop(columns=["_lookup_key"])

    # Save to new file
    os.makedirs(args.out_dir, exist_ok=True)
    out_name = args.out_name or f"{os.path.splitext(os.path.basename(args.input))[0]}_with_predictions.csv"
    out_csv = os.path.join(args.out_dir, out_name)
    df.to_csv(out_csv, index=False)
    print(f"Created new file: {out_csv}")
    print(f"  - Total rows: {len(df)}")
    print(f"  - Rows with predictions: {df['boltz_prob_mean'].notna().sum()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
