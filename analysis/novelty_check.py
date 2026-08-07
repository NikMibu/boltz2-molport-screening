#!/usr/bin/env python3
"""Structural novelty check for a hit list.

For every hit, computes the maximum Tanimoto similarity (Morgan, radius 2,
2048 bits) against a reference set of known inhibitors. A low maximum means
the compound resembles nothing in the reference set, i.e. it is structurally
novel.

The hit list is not shipped with this repository; pass your own via --hits.
The reference set under data/reference/ is derived from ChEMBL and public.
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs

RDLogger.DisableLog("rdApp.*")

# Novelty bins. Kept fixed so results stay comparable across runs.
BINS = [0, 0.3, 0.5, 0.7, 0.85, 1.0]
LABELS = [
    "Novel (<0.3)",
    "Low similarity (0.3-0.5)",
    "Moderate (0.5-0.7)",
    "Similar (0.7-0.85)",
    "Very similar (>=0.85)",
]
NOVELTY_THRESHOLD = 0.3

DEFAULT_KNOWN = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "reference", "known_ca2_inhibitors.csv",
)


def parse_args():
    p = argparse.ArgumentParser(
        description="Tanimoto novelty check of a hit list against known inhibitors",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--hits", required=True,
                   help="CSV with the hit list; needs a 'smiles' or 'SMILES' column")
    p.add_argument("--known", default=DEFAULT_KNOWN,
                   help="CSV with the reference set of known inhibitors")
    p.add_argument("--out-dir", default=".",
                   help="Directory for the CSV and the plots")
    p.add_argument("--score-col", default=None,
                   help="Score column to plot against novelty; auto-detected if omitted")
    p.add_argument("--radius", type=int, default=2, help="Morgan fingerprint radius")
    p.add_argument("--n-bits", type=int, default=2048, help="Morgan fingerprint length")
    return p.parse_args()


def smiles_column(df, path):
    for col in ("smiles", "SMILES", "canonical_smiles"):
        if col in df.columns:
            return col
    raise KeyError(f"{path} needs a 'smiles' or 'SMILES' column, found: {list(df.columns)}")


def fingerprint(smiles, radius, n_bits):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    hits = pd.read_csv(args.hits)
    known = pd.read_csv(args.known)
    hits_col = smiles_column(hits, args.hits)
    known_col = smiles_column(known, args.known)
    print(f"Hits:      {len(hits)} rows from {args.hits}")
    print(f"Reference: {len(known)} rows from {args.known}")

    score_col = args.score_col
    if score_col is None:
        score_col = next(
            (c for c in ("boltz_pred_value_mean", "pred_score", "prediction_score")
             if c in hits.columns),
            None,
        )
    if args.score_col and args.score_col not in hits.columns:
        raise KeyError(f"--score-col '{args.score_col}' not in {args.hits}")

    hit_fps, hit_smiles, hit_scores = [], [], []
    for _, row in hits.iterrows():
        smi = str(row[hits_col])
        fp = fingerprint(smi, args.radius, args.n_bits)
        if fp is not None:
            hit_fps.append(fp)
            hit_smiles.append(smi)
            hit_scores.append(row[score_col] if score_col else np.nan)

    known_fps = []
    known_smiles = set()
    for smi in known[known_col].astype(str):
        fp = fingerprint(smi, args.radius, args.n_bits)
        if fp is not None:
            known_fps.append(fp)
            known_smiles.add(smi)

    print(f"Valid fingerprints: {len(hit_fps)} hits, {len(known_fps)} reference")
    if not hit_fps or not known_fps:
        raise SystemExit("ERROR: no valid fingerprints on one side, nothing to compare")

    max_sim = np.array([
        max(DataStructs.BulkTanimotoSimilarity(fp, known_fps)) for fp in hit_fps
    ])

    print("\nMax Tanimoto to the reference set")
    print(f"  mean {max_sim.mean():.3f} | median {np.median(max_sim):.3f} "
          f"| min {max_sim.min():.3f} | max {max_sim.max():.3f}")

    hist, _ = np.histogram(max_sim, bins=BINS)
    print("\nNovelty classification")
    for label, count in zip(LABELS, hist):
        print(f"  {label:26s} {count:4d} / {len(max_sim)} ({100 * count / len(max_sim):5.1f} %)")

    print(f"\nMost novel:  Tanimoto {max_sim.min():.3f}  {hit_smiles[int(np.argmin(max_sim))]}")
    print(f"Least novel: Tanimoto {max_sim.max():.3f}  {hit_smiles[int(np.argmax(max_sim))]}")

    results = pd.DataFrame({
        "SMILES": hit_smiles,
        "Max_Tanimoto_to_Known": max_sim,
    })
    if score_col:
        results[score_col] = hit_scores
    results = results.sort_values("Max_Tanimoto_to_Known")
    csv_path = os.path.join(args.out_dir, "novelty_check_results.csv")
    results.to_csv(csv_path, index=False)
    print(f"\nWrote {csv_path}")

    # --- Plots ---
    sns.set_style("whitegrid")
    plot_df = results.copy()
    plot_df["Novelty_Class"] = pd.cut(
        plot_df["Max_Tanimoto_to_Known"], bins=BINS, labels=LABELS, include_lowest=True
    )
    class_counts = plot_df["Novelty_Class"].value_counts().reindex(LABELS, fill_value=0)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].hist(plot_df["Max_Tanimoto_to_Known"], bins=20,
                 color="steelblue", edgecolor="black", alpha=0.8)
    axes[0].axvline(NOVELTY_THRESHOLD, color="darkred", linestyle="--", linewidth=2,
                    label=f"Novelty threshold ({NOVELTY_THRESHOLD})")
    axes[0].set_xlabel("Max Tanimoto to reference set")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Structural novelty distribution")
    axes[0].legend()

    ax = sns.barplot(x=class_counts.index, y=class_counts.values,
                     hue=class_counts.index, palette="viridis", legend=False, ax=axes[1])
    for i, count in enumerate(class_counts.values):
        ax.text(i, count + max(class_counts.max() * 0.02, 0.2), f"{count}",
                ha="center", va="bottom", fontsize=10)
    axes[1].set_xlabel("Novelty class")
    axes[1].set_ylabel("Count")
    axes[1].set_title("Novelty classification of the hit list")
    axes[1].tick_params(axis="x", rotation=20)
    for tick in axes[1].get_xticklabels():
        tick.set_horizontalalignment("right")

    plt.tight_layout()
    overview = os.path.join(args.out_dir, "novelty_overview.png")
    plt.savefig(overview, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Wrote {overview}")

    if score_col and plot_df[score_col].notna().any():
        plt.figure(figsize=(8, 5))
        scatter = plt.scatter(
            plot_df["Max_Tanimoto_to_Known"], plot_df[score_col],
            c=plot_df["Max_Tanimoto_to_Known"], cmap="coolwarm",
            edgecolors="black", linewidths=0.4, alpha=0.8,
        )
        plt.axvline(NOVELTY_THRESHOLD, color="darkred", linestyle="--",
                    linewidth=2, alpha=0.8)
        plt.xlabel("Max Tanimoto to reference set")
        plt.ylabel(score_col)
        plt.title(f"{score_col} vs. structural novelty")
        plt.colorbar(scatter, label="Max Tanimoto")
        plt.tight_layout()
        scatter_path = os.path.join(args.out_dir, "novelty_vs_pred_score.png")
        plt.savefig(scatter_path, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Wrote {scatter_path}")

    # --- Summary ---
    # Both figures are computed here. The exact-match count used to be a
    # hard-coded literal, which silently went stale whenever the input changed.
    novel_count = int(np.sum(max_sim < NOVELTY_THRESHOLD))
    exact_matches = sum(1 for smi in hit_smiles if smi in known_smiles)
    n = len(max_sim)

    print("\nSummary")
    print(f"  {exact_matches}/{n} hits are an exact SMILES match to the reference set")
    print(f"  {novel_count}/{n} hits ({100 * novel_count / n:.0f} %) are structurally novel "
          f"(Tanimoto < {NOVELTY_THRESHOLD} to every reference compound)")


if __name__ == "__main__":
    main()
