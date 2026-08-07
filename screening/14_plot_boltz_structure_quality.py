#!/usr/bin/env python3
"""Plot Boltz structure-quality metrics per isoform."""

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


METRICS = [
    "confidence_score",
    "ptm",
    "iptm",
    "ligand_iptm",
    "complex_plddt",
    "complex_iplddt",
    "complex_pde",
    "complex_ipde",
]

KEY_METRICS = ["confidence_score", "ptm", "iptm", "ligand_iptm", "complex_plddt"]


def extract_isoform(path: Path) -> str:
    """Extract isoform from any '*_predictions' path segment."""
    for part in path.parts:
        match = re.match(r"(ca\d+)_predictions$", part, flags=re.IGNORECASE)
        if match:
            return match.group(1).lower()
    return "unknown"


def extract_molport_and_model(conf_path: Path) -> tuple[str, int | None]:
    """Extract MolPort ID and model index from confidence file name."""
    name = conf_path.name
    match = re.match(r"confidence_.*?_(MolPort-[\d-]+)_model_(\d+)\.json$", name, flags=re.IGNORECASE)
    if match:
        return match.group(1), int(match.group(2))
    fallback = re.search(r"(MolPort-[\d-]+)", name, flags=re.IGNORECASE)
    molport = fallback.group(1) if fallback else conf_path.parent.name
    model_match = re.search(r"_model_(\d+)\.json$", name)
    model = int(model_match.group(1)) if model_match else None
    return molport, model


def collect_rows(results_root: Path) -> pd.DataFrame:
    """Collect all confidence JSONs in results_root."""
    rows = []
    for conf_path in results_root.rglob("confidence_*.json"):
        try:
            data = json.loads(conf_path.read_text())
        except Exception:
            continue
        molport, model = extract_molport_and_model(conf_path)
        row = {
            "isoform": extract_isoform(conf_path),
            "molport_id": molport,
            "model_idx": model,
            "file": str(conf_path),
        }
        for metric in METRICS:
            value = data.get(metric)
            row[metric] = float(value) if isinstance(value, (int, float)) else None
        rows.append(row)
    return pd.DataFrame(rows)


def plot_isoform_distribution(df_iso: pd.DataFrame, out_png: Path) -> None:
    """Create boxplot with individual points for one isoform."""
    available = [m for m in METRICS if m in df_iso and df_iso[m].notna().any()]
    if not available:
        return

    melt = df_iso[available].melt(var_name="metric", value_name="value").dropna()
    fig, ax = plt.subplots(figsize=(12, 6))
    grouped = [melt.loc[melt["metric"] == metric, "value"].values for metric in available]
    ax.boxplot(grouped, tick_labels=available, showfliers=False)

    for idx, metric in enumerate(available, start=1):
        vals = melt.loc[melt["metric"] == metric, "value"].values
        jitter = pd.Series(range(len(vals))).mod(7).mul(0.03).sub(0.09).values
        x = idx + jitter
        ax.scatter(x, vals, s=12, alpha=0.35)

    ax.set_title(f"Boltz structure quality metrics - {df_iso['isoform'].iloc[0]}")
    ax.set_ylabel("Metric value")
    ax.set_ylim(0.0, 1.05)
    ax.grid(axis="y", alpha=0.2)
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_isoform_scatter(df_iso: pd.DataFrame, out_png: Path) -> None:
    """Plot pTM vs ipTM with confidence color-coding."""
    if "ptm" not in df_iso or "iptm" not in df_iso:
        return
    valid = df_iso[["ptm", "iptm", "confidence_score"]].dropna(subset=["ptm", "iptm"])
    if valid.empty:
        return

    fig, ax = plt.subplots(figsize=(7, 6))
    color = valid["confidence_score"] if "confidence_score" in valid else None
    sc = ax.scatter(valid["ptm"], valid["iptm"], c=color, cmap="viridis", s=20, alpha=0.75)
    if color is not None and not color.isna().all():
        cbar = fig.colorbar(sc, ax=ax)
        cbar.set_label("confidence_score")
    ax.set_title(f"pTM vs ipTM - {df_iso['isoform'].iloc[0]}")
    ax.set_xlabel("pTM")
    ax.set_ylabel("ipTM")
    ax.set_xlim(0.0, 1.02)
    ax.set_ylim(0.0, 1.02)
    ax.grid(alpha=0.25)
    plt.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_isoform_comparison(df: pd.DataFrame, out_png: Path) -> None:
    """Create per-metric isoform distributions (violin + box + jitter)."""
    available = [m for m in KEY_METRICS if m in df and df[m].notna().any()]
    if not available:
        return

    isoforms = sorted([iso for iso in df["isoform"].dropna().unique() if iso != "unknown"])
    if not isoforms:
        return

    fig, axes = plt.subplots(1, len(available), figsize=(4.0 * len(available), 4.7), sharey=True)
    if len(available) == 1:
        axes = [axes]

    for ax, metric in zip(axes, available):
        iso_vals = [(iso, df.loc[df["isoform"] == iso, metric].dropna().values) for iso in isoforms]
        iso_vals = [(iso, vals) for iso, vals in iso_vals if len(vals) > 0]
        if not iso_vals:
            continue

        labels = [iso for iso, _ in iso_vals]
        grouped = [vals for _, vals in iso_vals]
        positions = list(range(1, len(labels) + 1))
        violin = ax.violinplot(grouped, positions=positions, widths=0.85, showmeans=False, showmedians=False, showextrema=False)
        for body in violin["bodies"]:
            body.set_alpha(0.32)

        ax.boxplot(
            grouped,
            positions=positions,
            widths=0.26,
            showfliers=False,
            patch_artist=True,
            boxprops={"facecolor": "white", "edgecolor": "#1f1f1f", "alpha": 0.9},
            whiskerprops={"color": "#1f1f1f"},
            capprops={"color": "#1f1f1f"},
            medianprops={"color": "#cc2f2f", "linewidth": 1.6},
        )

        ax.set_title(metric)
        ax.set_xlabel("isoform")
        ax.set_xticks(positions)
        ax.set_xticklabels(labels)
        ax.set_ylim(0.0, 1.05)
        ax.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("metric value")
    fig.suptitle("Boltz structure quality distributions per isoform")
    plt.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def write_summary(df: pd.DataFrame, out_csv: Path) -> None:
    """Write compact summary stats per isoform."""
    agg = {}
    for metric in METRICS:
        if metric not in df:
            continue
        agg[metric] = ["count", "mean", "median", "std", "min", "max"]
    summary = df.groupby("isoform").agg(agg)
    summary.columns = [f"{metric}_{stat}" for metric, stat in summary.columns]
    summary = summary.reset_index()
    summary.to_csv(out_csv, index=False)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Plot Boltz structure quality metrics per isoform")
    ap.add_argument("--results-root", help="Root folder containing Boltz prediction outputs")
    ap.add_argument("--input-csv", help="Optional precomputed rows CSV (boltz_structure_quality_all_rows.csv)")
    ap.add_argument("--out-dir", required=True, help="Output folder for plots and CSVs")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.input_csv:
        df = pd.read_csv(args.input_csv)
        all_csv = out_dir / "boltz_structure_quality_all_rows.csv"
        if Path(args.input_csv).resolve() != all_csv.resolve():
            df.to_csv(all_csv, index=False)
    elif args.results_root:
        results_root = Path(args.results_root)
        df = collect_rows(results_root)
        all_csv = out_dir / "boltz_structure_quality_all_rows.csv"
        df.to_csv(all_csv, index=False)
    else:
        raise SystemExit("Either --results-root or --input-csv is required")

    if df.empty:
        print("No rows found for plotting.")
        return

    df = df.sort_values(["isoform", "molport_id", "model_idx"], na_position="last").reset_index(drop=True)

    summary_csv = out_dir / "boltz_structure_quality_summary_by_isoform.csv"
    write_summary(df, summary_csv)

    for isoform, df_iso in df.groupby("isoform", sort=True):
        if isoform == "unknown":
            continue
        plot_isoform_distribution(df_iso, out_dir / f"{isoform}_structure_quality_distribution.png")
        plot_isoform_scatter(df_iso, out_dir / f"{isoform}_ptm_vs_iptm.png")

    plot_isoform_comparison(df, out_dir / "isoform_structure_quality_comparison.png")

    print(f"Rows analyzed: {len(df)}")
    print(f"Isoforms found: {', '.join(sorted(df['isoform'].dropna().unique()))}")
    print(f"Saved: {all_csv}")
    print(f"Saved: {summary_csv}")
    print(f"Plots written to: {out_dir}")


if __name__ == "__main__":
    main()
