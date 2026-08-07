#!/usr/bin/env python3
"""Exact SMILES overlap between a hit list and a reference set.

Complements novelty_check.py: this counts identical strings, novelty_check.py
measures structural similarity. The hit list is not shipped with this
repository; pass your own via --hits.
"""
import argparse
import pandas as pd
import os

DEFAULT_KNOWN = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "reference", "known_ca2_inhibitors.csv",
)


def load_smiles(path):
    df = pd.read_csv(path)
    if "smiles" in df.columns:
        col = "smiles"
    elif "SMILES" in df.columns:
        col = "SMILES"
    else:
        raise KeyError(f"{path} needs a 'smiles' or 'SMILES' column")
    return df, df[col].astype(str)


def main():
    p = argparse.ArgumentParser(
        description="Count exact SMILES matches between a hit list and a reference set",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--hits", required=True,
                   help="CSV with the hit list; needs a 'smiles' or 'SMILES' column")
    p.add_argument("--known", default=DEFAULT_KNOWN,
                   help="CSV with the reference set of known inhibitors")
    p.add_argument("--out-overlap", default="smiles_overlap.csv",
                   help="Where to write the overlapping SMILES, if any")
    args = p.parse_args()

    known_df, known_smiles = load_smiles(args.known)
    hits_df, hits_smiles = load_smiles(args.hits)

    known_set = set(known_smiles)
    hits_set = set(hits_smiles)
    overlap = known_set & hits_set

    print(f"Known set: {len(known_set)} unique SMILES")
    print(f"Hits set:  {len(hits_set)} unique SMILES")
    print(f"Overlap:   {len(overlap)} SMILES")

    if overlap:
        overlap_df = pd.DataFrame({"smiles": sorted(overlap)})
        overlap_df.to_csv(args.out_overlap, index=False)
        print(f"Saved overlap SMILES to {os.path.abspath(args.out_overlap)}")


if __name__ == "__main__":
    main()

