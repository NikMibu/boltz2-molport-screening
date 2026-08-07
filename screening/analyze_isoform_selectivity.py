#!/usr/bin/env python3
"""Analyze isoform selectivity and overlap between CA2, CA4, and CA7."""
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, Tuple


def load_isoform_data(ca2_path: str, ca4_path: str, ca7_path: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load ranked ligand data for all three isoforms."""
    print("Loading isoform data...")
    ca2 = pd.read_csv(ca2_path)
    ca4 = pd.read_csv(ca4_path)
    ca7 = pd.read_csv(ca7_path)
    
    print(f"  CA2: {len(ca2)} compounds")
    print(f"  CA4: {len(ca4)} compounds")
    print(f"  CA7: {len(ca7)} compounds")
    
    return ca2, ca4, ca7


def merge_isoform_data(ca2: pd.DataFrame, ca4: pd.DataFrame, ca7: pd.DataFrame) -> pd.DataFrame:
    """Merge data from all isoforms on molport_number."""
    # Select key columns and rename for each isoform
    ca2_sel = ca2[['molport_number', 'boltz_prob_mean', 'boltz_pred_value_mean', 'combined_score', 'rank']].copy()
    ca2_sel.columns = ['molport_number', 'ca2_prob', 'ca2_pred', 'ca2_score', 'ca2_rank']
    
    ca4_sel = ca4[['molport_number', 'boltz_prob_mean', 'boltz_pred_value_mean', 'combined_score', 'rank']].copy()
    ca4_sel.columns = ['molport_number', 'ca4_prob', 'ca4_pred', 'ca4_score', 'ca4_rank']
    
    ca7_sel = ca7[['molport_number', 'boltz_prob_mean', 'boltz_pred_value_mean', 'combined_score', 'rank']].copy()
    ca7_sel.columns = ['molport_number', 'ca7_prob', 'ca7_pred', 'ca7_score', 'ca7_rank']
    
    # Merge all isoforms
    merged = ca2_sel.merge(ca4_sel, on='molport_number', how='outer')
    merged = merged.merge(ca7_sel, on='molport_number', how='outer')
    
    print(f"\nMerged data: {len(merged)} unique compounds")
    print(f"  Present in CA2: {merged['ca2_score'].notna().sum()}")
    print(f"  Present in CA4: {merged['ca4_score'].notna().sum()}")
    print(f"  Present in CA7: {merged['ca7_score'].notna().sum()}")
    
    return merged


def calculate_selectivity(merged: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate isoform selectivity metrics.
    
    Note: With new scoring system, LOWER score = BETTER binding.
    Selectivity is calculated as typical : 
    - Best isoform has LOWEST score (strongest binding)
    - Selectivity ratio = second_best_score / best_score (higher = more selective)
    """
    df = merged.copy()
    
    # Filter out flagged compounds (score > 900) and replace with NaN
    scores = df[['ca2_score', 'ca4_score', 'ca7_score']].copy()
    scores = scores.mask(scores > 900, np.nan)
    
    # Calculate basic score metrics
    df['min_score'] = scores.min(axis=1)  # Best (lowest) score
    df['max_score'] = scores.max(axis=1)  # Worst (highest) score
    df['score_range'] = df['max_score'] - df['min_score']
    df['score_mean'] = scores.mean(axis=1)
    df['score_std'] = scores.std(axis=1)
    
    # Identify best isoform (LOWEST score = best)
    df['best_isoform'] = scores.idxmin(axis=1).str.replace('_score', '').str.upper()
    df['best_score'] = df['min_score']
    
    # Calculate selectivity ratio (second_best / best)
    # Higher ratio = more selective
    def calc_selectivity_ratio(row):
        valid_scores = row[row.notna()].values
        if len(valid_scores) < 2:
            return np.nan
        sorted_scores = np.sort(valid_scores)  # Ascending: [best, ..., worst]
        best = sorted_scores[0]
        second_best = sorted_scores[1]
        if best <= 0 or second_best <= 0:
            return np.nan
        # Selectivity ratio: how much worse is second best compared to best?
        return second_best / best
    
    df['selectivity_ratio'] = scores.apply(calc_selectivity_ratio, axis=1)
    
    # Calculate fold-change (alternative selectivity metric)
    # Difference between best and second-best
    def calc_fold_change(row):
        valid_scores = row[row.notna()].values
        if len(valid_scores) < 2:
            return np.nan
        sorted_scores = np.sort(valid_scores)
        return sorted_scores[1] - sorted_scores[0]  # Second_best - best
    
    df['selectivity_fold_change'] = scores.apply(calc_fold_change, axis=1)
    
    # Count valid isoforms (exclude flagged)
    df['isoform_count'] = scores.notna().sum(axis=1)
    
    # Calculate worst isoform (for comparison)
    df['worst_isoform'] = scores.idxmax(axis=1).str.replace('_score', '').str.upper()
    df['worst_score'] = df['max_score']
    
    return df


def identify_overlap(merged: pd.DataFrame, top_n: int = 100) -> Dict:
    """Identify overlap between top N compounds of each isoform."""
    # Get top N for each isoform
    top_ca2 = set(merged.nsmallest(top_n, 'ca2_rank')['molport_number'].dropna())
    top_ca4 = set(merged.nsmallest(top_n, 'ca4_rank')['molport_number'].dropna())
    top_ca7 = set(merged.nsmallest(top_n, 'ca7_rank')['molport_number'].dropna())
    
    overlap = {
        'top_n': top_n,
        'ca2_only': len(top_ca2 - top_ca4 - top_ca7),
        'ca4_only': len(top_ca4 - top_ca2 - top_ca7),
        'ca7_only': len(top_ca7 - top_ca2 - top_ca4),
        'ca2_ca4': len(top_ca2 & top_ca4 - top_ca7),
        'ca2_ca7': len(top_ca2 & top_ca7 - top_ca4),
        'ca4_ca7': len(top_ca4 & top_ca7 - top_ca2),
        'all_three': len(top_ca2 & top_ca4 & top_ca7),
        'total_unique': len(top_ca2 | top_ca4 | top_ca7),
    }
    
    return overlap


def plot_selectivity_analysis(df: pd.DataFrame, out_dir: str):
    """Create visualization plots for selectivity analysis."""
    # Filter compounds present in at least one isoform (and not flagged)
    df_plot = df[df[['ca2_score', 'ca4_score', 'ca7_score']].notna().any(axis=1)].copy()
    df_plot = df_plot[(df_plot[['ca2_score', 'ca4_score', 'ca7_score']] < 900).any(axis=1)]
    
    fig = plt.figure(figsize=(18, 12))
    
    # 1. Scatter: CA2 vs CA4 scores (lower = better)
    ax1 = plt.subplot(2, 3, 1)
    mask = (df_plot[['ca2_score', 'ca4_score']] < 900).all(axis=1)
    scatter1 = ax1.scatter(df_plot.loc[mask, 'ca2_score'], df_plot.loc[mask, 'ca4_score'],
                          alpha=0.4, s=25, c=df_plot.loc[mask, 'selectivity_ratio'],
                          cmap='RdYlGn', edgecolors='black', linewidth=0.3)
    # Add diagonal line
    lims = [min(df_plot.loc[mask, 'ca2_score'].min(), df_plot.loc[mask, 'ca4_score'].min()),
            max(df_plot.loc[mask, 'ca2_score'].max(), df_plot.loc[mask, 'ca4_score'].max())]
    ax1.plot(lims, lims, 'r--', alpha=0.5, linewidth=2, label='Equal affinity')
    ax1.set_xlabel('CA2 Score (lower = better)', fontsize=10)
    ax1.set_ylabel('CA4 Score (lower = better)', fontsize=10)
    ax1.set_title('CA2 vs CA4 Selectivity', fontsize=11, fontweight='bold')
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    plt.colorbar(scatter1, ax=ax1, label='Selectivity Ratio')
    
    # 2. Scatter: CA2 vs CA7 scores
    ax2 = plt.subplot(2, 3, 2)
    mask = (df_plot[['ca2_score', 'ca7_score']] < 900).all(axis=1)
    scatter2 = ax2.scatter(df_plot.loc[mask, 'ca2_score'], df_plot.loc[mask, 'ca7_score'],
                          alpha=0.4, s=25, c=df_plot.loc[mask, 'selectivity_ratio'],
                          cmap='RdYlGn', edgecolors='black', linewidth=0.3)
    lims = [min(df_plot.loc[mask, 'ca2_score'].min(), df_plot.loc[mask, 'ca7_score'].min()),
            max(df_plot.loc[mask, 'ca2_score'].max(), df_plot.loc[mask, 'ca7_score'].max())]
    ax2.plot(lims, lims, 'r--', alpha=0.5, linewidth=2, label='Equal affinity')
    ax2.set_xlabel('CA2 Score (lower = better)', fontsize=10)
    ax2.set_ylabel('CA7 Score (lower = better)', fontsize=10)
    ax2.set_title('CA2 vs CA7 Selectivity', fontsize=11, fontweight='bold')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    plt.colorbar(scatter2, ax=ax2, label='Selectivity Ratio')
    
    # 3. Scatter: CA4 vs CA7 scores
    ax3 = plt.subplot(2, 3, 3)
    mask = (df_plot[['ca4_score', 'ca7_score']] < 900).all(axis=1)
    scatter3 = ax3.scatter(df_plot.loc[mask, 'ca4_score'], df_plot.loc[mask, 'ca7_score'],
                          alpha=0.4, s=25, c=df_plot.loc[mask, 'selectivity_ratio'],
                          cmap='RdYlGn', edgecolors='black', linewidth=0.3)
    lims = [min(df_plot.loc[mask, 'ca4_score'].min(), df_plot.loc[mask, 'ca7_score'].min()),
            max(df_plot.loc[mask, 'ca4_score'].max(), df_plot.loc[mask, 'ca7_score'].max())]
    ax3.plot(lims, lims, 'r--', alpha=0.5, linewidth=2, label='Equal affinity')
    ax3.set_xlabel('CA4 Score (lower = better)', fontsize=10)
    ax3.set_ylabel('CA7 Score (lower = better)', fontsize=10)
    ax3.set_title('CA4 vs CA7 Selectivity', fontsize=11, fontweight='bold')
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)
    plt.colorbar(scatter3, ax=ax3, label='Selectivity Ratio')
    
    # 4. Selectivity ratio distribution 
    ax4 = plt.subplot(2, 3, 4)
    selectivity = df_plot['selectivity_ratio'].dropna()
    selectivity = selectivity[selectivity < 100]  # Remove extreme outliers for visualization
    ax4.hist(selectivity, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    ax4.set_xlabel('Selectivity Ratio (2nd_best / best)', fontsize=10)
    ax4.set_ylabel('Frequency')
    ax4.set_title('Selectivity Ratio Distribution\n(higher = more selective)', fontsize=11, fontweight='bold')
    ax4.axvline(selectivity.median(), color='red', linestyle='--', linewidth=2,
                label=f'Median: {selectivity.median():.2f}')
    # Add selectivity thresholds
    ax4.axvline(10, color='green', linestyle=':', linewidth=2, alpha=0.7, label='10-fold selective')
    ax4.axvline(100, color='darkgreen', linestyle=':', linewidth=2, alpha=0.7, label='100-fold selective')
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3)
    
    # 5. Best isoform distribution
    ax5 = plt.subplot(2, 3, 5)
    best_counts = df_plot['best_isoform'].value_counts()
    colors_iso = {'CA2': 'steelblue', 'CA4': 'coral', 'CA7': 'mediumseagreen'}
    bar_colors = [colors_iso.get(iso, 'gray') for iso in best_counts.index]
    bars = ax5.bar(best_counts.index, best_counts.values, color=bar_colors, alpha=0.7, edgecolor='black')
    ax5.set_xlabel('Best Isoform (lowest score)', fontsize=10)
    ax5.set_ylabel('Count')
    ax5.set_title('Compounds by Best Isoform', fontsize=11, fontweight='bold')
    ax5.grid(True, alpha=0.3, axis='y')
    # Add counts on bars
    for bar in bars:
        height = bar.get_height()
        ax5.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}', ha='center', va='bottom', fontsize=10)
    
    # 6. Fold-change distribution (alternative selectivity metric)
    ax6 = plt.subplot(2, 3, 6)
    fold_change = df_plot['selectivity_fold_change'].dropna()
    fold_change = fold_change[fold_change < 5]  # Remove extreme outliers
    ax6.hist(fold_change, bins=50, edgecolor='black', alpha=0.7, color='mediumpurple')
    ax6.set_xlabel('Selectivity Fold-Change\n(2nd_best - best)', fontsize=10)
    ax6.set_ylabel('Frequency')
    ax6.set_title('Fold-Change Distribution\n(higher = more selective)', fontsize=11, fontweight='bold')
    ax6.axvline(fold_change.median(), color='red', linestyle='--', linewidth=2,
                label=f'Median: {fold_change.median():.2f}')
    ax6.legend(fontsize=9)
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_path = os.path.join(out_dir, 'isoform_selectivity_analysis.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Saved selectivity plots to: {plot_path}")
    plt.close()


def plot_overlap_venn(overlap: Dict, out_dir: str):
    """Create Venn diagram for top N overlap."""
    try:
        from matplotlib_venn import venn3
    except ImportError:
        print("WARNING: matplotlib_venn not installed. Skipping Venn diagram.")
        print("  Install with: pip install matplotlib-venn")
        return
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    venn3(subsets=(
        overlap['ca2_only'],
        overlap['ca4_only'],
        overlap['ca2_ca4'],
        overlap['ca7_only'],
        overlap['ca2_ca7'],
        overlap['ca4_ca7'],
        overlap['all_three']
    ), set_labels=('CA2', 'CA4', 'CA7'), ax=ax)
    
    ax.set_title(f'Top {overlap["top_n"]} Compounds Overlap')
    plt.tight_layout()
    
    plot_path = os.path.join(out_dir, 'isoform_overlap_venn.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Saved Venn diagram to: {plot_path}")
    plt.close()


def plot_overlap_bar(overlap: Dict, out_dir: str):
    """Create bar chart for overlap statistics."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    categories = [
        'CA2 only', 'CA4 only', 'CA7 only',
        'CA2+CA4', 'CA2+CA7', 'CA4+CA7',
        'All three'
    ]
    values = [
        overlap['ca2_only'], overlap['ca4_only'], overlap['ca7_only'],
        overlap['ca2_ca4'], overlap['ca2_ca7'], overlap['ca4_ca7'],
        overlap['all_three']
    ]
    colors = ['blue', 'green', 'orange', 'cyan', 'magenta', 'yellow', 'red']
    
    bars = ax.bar(categories, values, color=colors, alpha=0.7, edgecolor='black')
    ax.set_ylabel('Number of Compounds')
    ax.set_title(f'Top {overlap["top_n"]} Compounds Overlap Analysis')
    ax.grid(True, alpha=0.3, axis='y')
    ax.tick_params(axis='x', rotation=45)
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}', ha='center', va='bottom')
    
    plt.tight_layout()
    plot_path = os.path.join(out_dir, 'isoform_overlap_bar.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Saved overlap bar chart to: {plot_path}")
    plt.close()


def main():
    ap = argparse.ArgumentParser(
        description="Analyze isoform selectivity and overlap between CA2, CA4, CA7",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Selectivity Metrics:
  - Lower score = better binding (stronger affinity)
  - Best isoform = isoform with LOWEST score
  - Selectivity ratio = second_best_score / best_score (higher = more selective)
  - Fold-change = second_best_score - best_score (higher = more selective)
  
  Typical selectivity:
    - >10-fold: Moderately selective
    - >100-fold: Highly selective
    - >1000-fold: Extremely selective

Examples:
  # Basic analysis
  python analyze_isoform_selectivity.py \\
    --ca2 results/analyze/analysis_ca2/ligand_analysis_filtered/ligands_ranked_full.csv \\
    --ca4 results/analyze/analysis_ca4/ligand_analysis_ca4_filtered/ligands_ranked_full.csv \\
    --ca7 results/analyze/analysis_ca7/ligand_analysis_ca7_filtered/ligands_ranked_full.csv \\
    --out-dir results/analyze/isoform_selectivity

  # With custom top N for overlap
  python analyze_isoform_selectivity.py \\
    --ca2 data_ca2.csv --ca4 data_ca4.csv --ca7 data_ca7.csv \\
    --out-dir results/ --top-n 50
        """
    )
    
    ap.add_argument("--ca2", required=True,
                    help="Path to CA2 ligands_ranked_full.csv")
    ap.add_argument("--ca4", required=True,
                    help="Path to CA4 ligands_ranked_full.csv")
    ap.add_argument("--ca7", required=True,
                    help="Path to CA7 ligands_ranked_full.csv")
    ap.add_argument("--out-dir", required=True,
                    help="Output directory for results")
    ap.add_argument("--top-n", type=int, default=100,
                    help="Top N compounds for overlap analysis (default: 100)")
    
    args = ap.parse_args()
    
    # Validate input files
    for isoform, path in [("CA2", args.ca2), ("CA4", args.ca4), ("CA7", args.ca7)]:
        if not os.path.isfile(path):
            print(f"ERROR: {isoform} file not found: {path}")
            return
    
    # Create output directory
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Load data
    ca2, ca4, ca7 = load_isoform_data(args.ca2, args.ca4, args.ca7)
    
    # Merge data
    merged = merge_isoform_data(ca2, ca4, ca7)
    
    # Calculate selectivity
    print("\nCalculating selectivity metrics...")
    df_selectivity = calculate_selectivity(merged)
    
    # Identify overlap
    print(f"\nAnalyzing overlap in top {args.top_n} compounds...")
    overlap = identify_overlap(df_selectivity, top_n=args.top_n)
    
    # Print overlap statistics
    print("\n" + "="*60)
    print("Overlap Statistics")
    print("="*60)
    print(f"Top {overlap['top_n']} compounds:")
    print(f"  CA2 only:        {overlap['ca2_only']}")
    print(f"  CA4 only:        {overlap['ca4_only']}")
    print(f"  CA7 only:        {overlap['ca7_only']}")
    print(f"  CA2 + CA4:       {overlap['ca2_ca4']}")
    print(f"  CA2 + CA7:       {overlap['ca2_ca7']}")
    print(f"  CA4 + CA7:       {overlap['ca4_ca7']}")
    print(f"  All three:       {overlap['all_three']}")
    print(f"  Total unique:    {overlap['total_unique']}")
    print(f"  Overlap rate:    {overlap['all_three'] / overlap['total_unique'] * 100:.1f}%")
    
    # Print selectivity statistics
    print("\n" + "="*60)
    print("Selectivity Statistics")
    print("="*60)
    print(f"Compounds with selectivity ratio: {df_selectivity['selectivity_ratio'].notna().sum()}")
    if df_selectivity['selectivity_ratio'].notna().sum() > 0:
        sel = df_selectivity['selectivity_ratio'].dropna()
        print(f"\nSelectivity Ratio (2nd_best / best):")
        print(f"  Mean:   {sel.mean():.3f}-fold")
        print(f"  Median: {sel.median():.3f}-fold")
        print(f"  Max:    {sel.max():.3f}-fold")
        
        # Count highly selective compounds
        sel_10 = (sel >= 10).sum()
        sel_100 = (sel >= 100).sum()
        sel_1000 = (sel >= 1000).sum()
        print(f"\nSelectivity Classification:")
        print(f"  ≥10-fold selective:   {sel_10} ({sel_10/len(sel)*100:.1f}%)")
        print(f"  ≥100-fold selective:  {sel_100} ({sel_100/len(sel)*100:.1f}%)")
        print(f"  ≥1000-fold selective: {sel_1000} ({sel_1000/len(sel)*100:.1f}%)")
    
    if 'selectivity_fold_change' in df_selectivity.columns:
        fc = df_selectivity['selectivity_fold_change'].dropna()
        if len(fc) > 0:
            print(f"\nFold-Change (2nd_best - best):")
            print(f"  Mean:   {fc.mean():.3f}")
            print(f"  Median: {fc.median():.3f}")
    
    print(f"\nBest Isoform Distribution (lowest score = best):")
    best_dist = df_selectivity['best_isoform'].value_counts()
    for isoform, count in best_dist.items():
        print(f"  {isoform}: {count} ({count/len(df_selectivity)*100:.1f}%)")
    
    # Save results
    print("\nSaving results...")
    
    # Full merged data with selectivity
    merged_path = os.path.join(args.out_dir, 'isoform_selectivity_full.csv')
    df_selectivity.to_csv(merged_path, index=False)
    print(f"  Saved full data to: {merged_path}")
    
    # Top selective compounds (high selectivity ratio)
    top_selective = df_selectivity.nlargest(100, 'selectivity_ratio')
    top_selective_path = os.path.join(args.out_dir, 'top_selective_compounds.csv')
    top_selective.to_csv(top_selective_path, index=False)
    print(f"  Saved top selective compounds to: {top_selective_path}")
    
    # Compounds in all three top N
    all_three = df_selectivity[
        (df_selectivity['ca2_rank'] <= args.top_n) &
        (df_selectivity['ca4_rank'] <= args.top_n) &
        (df_selectivity['ca7_rank'] <= args.top_n)
    ].sort_values('score_mean', ascending=False)
    all_three_path = os.path.join(args.out_dir, f'compounds_in_all_top{args.top_n}.csv')
    all_three.to_csv(all_three_path, index=False)
    print(f"  Saved compounds in all top {args.top_n} to: {all_three_path}")
    
    # Save summary
    summary_path = os.path.join(args.out_dir, 'selectivity_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("="*60 + "\n")
        f.write("Isoform Selectivity Analysis Summary\n")
        f.write("="*60 + "\n\n")
        f.write(f"Input files:\n")
        f.write(f"  CA2: {args.ca2}\n")
        f.write(f"  CA4: {args.ca4}\n")
        f.write(f"  CA7: {args.ca7}\n\n")
        f.write(f"Total unique compounds: {len(df_selectivity)}\n")
        f.write(f"  Present in CA2: {df_selectivity['ca2_score'].notna().sum()}\n")
        f.write(f"  Present in CA4: {df_selectivity['ca4_score'].notna().sum()}\n")
        f.write(f"  Present in CA7: {df_selectivity['ca7_score'].notna().sum()}\n\n")
        f.write(f"Overlap in top {args.top_n}:\n")
        f.write(f"  CA2 only: {overlap['ca2_only']}\n")
        f.write(f"  CA4 only: {overlap['ca4_only']}\n")
        f.write(f"  CA7 only: {overlap['ca7_only']}\n")
        f.write(f"  CA2 + CA4: {overlap['ca2_ca4']}\n")
        f.write(f"  CA2 + CA7: {overlap['ca2_ca7']}\n")
        f.write(f"  CA4 + CA7: {overlap['ca4_ca7']}\n")
        f.write(f"  All three: {overlap['all_three']}\n")
        f.write(f"  Total unique: {overlap['total_unique']}\n")
        f.write(f"  Overlap rate: {overlap['all_three'] / overlap['total_unique'] * 100:.1f}%\n")
    print(f"  Saved summary to: {summary_path}")
    
    # Create plots
    print("\nGenerating plots...")
    plot_selectivity_analysis(df_selectivity, args.out_dir)
    plot_overlap_bar(overlap, args.out_dir)
    plot_overlap_venn(overlap, args.out_dir)
    
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
