#!/usr/bin/env python3
"""Analyze ligand predictions from Boltz and create rankings with filtering options."""
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def calculate_combined_score(df: pd.DataFrame, 
                             confidence_threshold: float = 0.5) -> pd.DataFrame:
    """
    Two-stage approach:
    1. Filter by confidence (QC) - only consider predictions with sufficient confidence
    2. Rank by predicted pKd - lower pred_value = stronger binding
    
    Args:
        confidence_threshold: Minimum boltz_prob_mean to consider (0.5 = median)
    
    Returns:
        DataFrame with 'passes_confidence' and 'combined_score' columns.
        Combined score is LOWER for better ligands (because lower pred_value = better).
        Low-confidence predictions get score of 999 (worst).
    """
    df_scored = df.copy()
    
    # Stage 1: Confidence QC
    df_scored['passes_confidence'] = (
        df_scored['boltz_prob_mean'].fillna(0) >= confidence_threshold
    )
    
    # Stage 2: For those that pass, use pred_value as score (lower = better)
    # Add small penalty for lower confidence as tiebreaker
    df_scored['combined_score'] = np.where(
        df_scored['passes_confidence'],
        df_scored['boltz_pred_value_mean'] - 0.1 * df_scored['boltz_prob_mean'],
        999  # Flag low-confidence predictions with worst score
    )
    
    return df_scored


def create_ranking(df: pd.DataFrame, sort_by: str = "combined_score", ascending: bool = True) -> pd.DataFrame:
    """
    Create ranking based on specified column.
    
    Note: With new combined_score, LOWER is BETTER (ascending=True by default).
    """
    df_ranked = df.copy()
    df_ranked = df_ranked.sort_values(by=sort_by, ascending=ascending, na_position='last')
    df_ranked['rank'] = range(1, len(df_ranked) + 1)
    return df_ranked


def plot_distributions(df: pd.DataFrame, out_dir: str):
    """Create distribution plots for Boltz predictions."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    
    # 1. Histogram of boltz_prob_mean (Confidence)
    ax1 = axes[0, 0]
    prob_data = df['boltz_prob_mean'].dropna()
    prob_data.hist(bins=50, ax=ax1, edgecolor='black', alpha=0.7, color='steelblue')
    ax1.set_xlabel('Confidence (boltz_prob_mean)', fontsize=10)
    ax1.set_ylabel('Frequency')
    ax1.set_title('Distribution of Confidence\n(higher = more confident)', fontsize=11)
    ax1.axvline(prob_data.median(), color='red', linestyle='--', linewidth=2,
                label=f'Median: {prob_data.median():.3f}')
    # Add confidence threshold line if passes_confidence exists
    if 'passes_confidence' in df.columns:
        threshold = df[df['passes_confidence']]['boltz_prob_mean'].min()
        if pd.notna(threshold):
            ax1.axvline(threshold, color='green', linestyle=':', linewidth=2,
                       label=f'QC threshold: {threshold:.3f}')
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # 2. Histogram of boltz_pred_value_mean (pKd-like)
    ax2 = axes[0, 1]
    pred_data = df['boltz_pred_value_mean'].dropna()
    pred_data.hist(bins=50, ax=ax2, edgecolor='black', alpha=0.7, color='coral')
    ax2.set_xlabel('Predicted Affinity (boltz_pred_value_mean)', fontsize=10)
    ax2.set_ylabel('Frequency')
    ax2.set_title('Distribution of Predicted Affinity\n(lower/negative = stronger binding)', fontsize=11)
    ax2.axvline(pred_data.median(), color='red', linestyle='--', linewidth=2,
                label=f'Median: {pred_data.median():.3f}')
    ax2.axvline(0, color='black', linestyle=':', linewidth=1.5, alpha=0.5,
                label='Zero line')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3, axis='y')
    
    # 3. Histogram of combined_score
    ax3 = axes[0, 2]
    if 'combined_score' in df.columns:
        # Filter out flagged compounds (score=999) for visualization
        valid_scores = df[df['combined_score'] < 900]['combined_score'].dropna()
        if len(valid_scores) > 0:
            valid_scores.hist(bins=50, ax=ax3, edgecolor='black', alpha=0.7, color='mediumpurple')
            ax3.set_xlabel('Combined Score', fontsize=10)
            ax3.set_ylabel('Frequency')
            ax3.set_title('Distribution of Combined Score (Lower = better)', fontsize=11, fontweight='bold')
            ax3.axvline(valid_scores.median(), color='red', linestyle='--', linewidth=2,
                        label=f'Median: {valid_scores.median():.3f}')
            ax3.legend(fontsize=9)
            ax3.grid(True, alpha=0.3, axis='y')
    
    # 4. Scatter plot: prob vs pred_value with confidence regions
    ax4 = axes[1, 0]
    # Separate high and low confidence
    if 'passes_confidence' in df.columns:
        high_conf = df[df['passes_confidence']]
        low_conf = df[~df['passes_confidence']]
        ax4.scatter(low_conf['boltz_prob_mean'], low_conf['boltz_pred_value_mean'],
                   alpha=0.3, s=15, c='lightgray', label='Low confidence', edgecolors='none')
        scatter = ax4.scatter(high_conf['boltz_prob_mean'], high_conf['boltz_pred_value_mean'],
                             alpha=0.6, s=20, c=high_conf['heavy_atoms'], cmap='viridis',
                             label='High confidence')
        plt.colorbar(scatter, ax=ax4, label='Heavy Atoms')
    else:
        scatter = ax4.scatter(df['boltz_prob_mean'], df['boltz_pred_value_mean'],
                             alpha=0.8, s=20, c=df['heavy_atoms'], cmap='viridis')
        plt.colorbar(scatter, ax=ax4, label='Heavy Atoms')
    ax4.set_xlabel('Confidence (boltz_prob_mean)', fontsize=10)
    ax4.set_ylabel('Predicted Affinity (boltz_pred_value_mean)', fontsize=10)
    ax4.set_title('Confidence vs Affinity', fontsize=11)
    ax4.axhline(0, color='red', linestyle='--', linewidth=1, alpha=0.5, label='y=0')
    ax4.legend(fontsize=8, loc='best')
    ax4.grid(True, alpha=0.3)
    
    # 5. Scatter plot with combined score (only high confidence)
    ax5 = axes[1, 1]
    if 'combined_score' in df.columns:
        # Only show compounds that pass QC
        mask = df['combined_score'] < 900
        if mask.sum() > 0:
            scatter2 = ax5.scatter(df.loc[mask, 'boltz_prob_mean'], 
                                  df.loc[mask, 'boltz_pred_value_mean'],
                                  alpha=0.6, s=25, c=df.loc[mask, 'combined_score'], 
                                  cmap='RdYlGn_r', edgecolors='black', linewidth=0.5)
            ax5.set_xlabel('Confidence (boltz_prob_mean)', fontsize=10)
            ax5.set_ylabel('Predicted Affinity (boltz_pred_value_mean)', fontsize=10)
            ax5.set_title('Compounds colored by Combined Score\n(green=low=best, red=high=worst)', fontsize=11)
            ax5.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.3)
            cbar = plt.colorbar(scatter2, ax=ax5, label='Combined Score')
            cbar.ax.invert_yaxis()  # Invert so best (low) is at top
            ax5.grid(True, alpha=0.3)
    
    # 6. Confidence QC statistics
    ax6 = axes[1, 2]
    if 'passes_confidence' in df.columns:
        # Pie chart of confidence QC
        pass_count = df['passes_confidence'].sum()
        fail_count = (~df['passes_confidence']).sum()
        colors = ['#90EE90', '#FFB6C6']  # Light green, light red
        wedges, texts, autotexts = ax6.pie([pass_count, fail_count], 
                                            labels=['Pass QC', 'Fail QC'],
                                            autopct='%1.1f%%',
                                            colors=colors,
                                            startangle=90,
                                            textprops={'fontsize': 10})
        ax6.set_title(f'Confidence QC Results\n(n={len(df)} \n Treshhold = {threshold:.3f} Probability)', fontsize=11)
        # Add counts
        for i, (count, autotext) in enumerate(zip([pass_count, fail_count], autotexts)):
            autotext.set_text(f'{autotext.get_text()}\n({count:,})')
    
    plt.tight_layout()
    plot_path = os.path.join(out_dir, 'boltz_predictions_distributions.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Saved distribution plots to: {plot_path}")
    plt.close()


def plot_property_correlations(df: pd.DataFrame, out_dir: str):
    """Create correlation plots between molecular properties and predictions."""
    properties = ['mw', 'heavy_atoms', 'logp', 'hbd', 'hba', 'rotatable_bonds']
    
    # Only plot for compounds that pass confidence QC
    if 'passes_confidence' in df.columns:
        df_plot = df[df['passes_confidence']].copy()
        title_suffix = ' (High Confidence Only)'
    else:
        df_plot = df.copy()
        title_suffix = ''
    
    fig, axes = plt.subplots(len(properties), 3, figsize=(15, 3.5 * len(properties)))
    
    for i, prop in enumerate(properties):
        # Column 1: Property vs Confidence
        ax1 = axes[i, 0]
        ax1.scatter(df_plot[prop], df_plot['boltz_prob_mean'], 
                   alpha=0.4, s=15, color='steelblue', edgecolors='none')
        ax1.set_xlabel(prop.replace('_', ' ').title(), fontsize=9)
        ax1.set_ylabel('Confidence', fontsize=9)
        ax1.set_title(f'{prop.upper()} vs Confidence{title_suffix}', fontsize=10)
        ax1.grid(True, alpha=0.3)
        # Add correlation
        valid = df_plot[[prop, 'boltz_prob_mean']].dropna()
        if len(valid) > 1:
            corr = valid.corr().iloc[0, 1]
            color = 'green' if abs(corr) < 0.3 else 'orange' if abs(corr) < 0.5 else 'red'
            ax1.text(0.05, 0.95, f'r = {corr:.3f}', 
                    transform=ax1.transAxes, va='top', fontsize=9,
                    bbox=dict(boxstyle='round', facecolor=color, alpha=0.3))
        
        # Column 2: Property vs Predicted Affinity
        ax2 = axes[i, 1]
        # Color by combined_score if available
        if 'combined_score' in df_plot.columns:
            mask = df_plot['combined_score'] < 900
            scatter = ax2.scatter(df_plot.loc[mask, prop], 
                                 df_plot.loc[mask, 'boltz_pred_value_mean'],
                                 alpha=0.5, s=20, 
                                 c=df_plot.loc[mask, 'combined_score'],
                                 cmap='RdYlGn_r', edgecolors='black', linewidth=0.3)
            if i == 0:  # Only add colorbar to top plot
                plt.colorbar(scatter, ax=ax2, label='Combined Score')
        else:
            ax2.scatter(df_plot[prop], df_plot['boltz_pred_value_mean'],
                       alpha=0.4, s=15, color='coral', edgecolors='none')
        ax2.set_xlabel(prop.replace('_', ' ').title(), fontsize=9)
        ax2.set_ylabel('Predicted Affinity', fontsize=9)
        ax2.set_title(f'{prop.upper()} vs Affinity{title_suffix}', fontsize=10)
        ax2.axhline(0, color='red', linestyle='--', linewidth=1, alpha=0.3)
        ax2.grid(True, alpha=0.3)
        # Add correlation
        valid = df_plot[[prop, 'boltz_pred_value_mean']].dropna()
        if len(valid) > 1:
            corr = valid.corr().iloc[0, 1]
            color = 'green' if abs(corr) < 0.3 else 'orange' if abs(corr) < 0.5 else 'red'
            ax2.text(0.05, 0.95, f'r = {corr:.3f}', 
                    transform=ax2.transAxes, va='top', fontsize=9,
                    bbox=dict(boxstyle='round', facecolor=color, alpha=0.3))
        
        # Column 3: Property vs Combined Score
        ax3 = axes[i, 2]
        if 'combined_score' in df_plot.columns:
            mask = df_plot['combined_score'] < 900
            if mask.sum() > 0:
                ax3.scatter(df_plot.loc[mask, prop], 
                           df_plot.loc[mask, 'combined_score'],
                           alpha=0.4, s=15, color='mediumpurple', edgecolors='none')
                ax3.set_xlabel(prop.replace('_', ' ').title(), fontsize=9)
                ax3.set_ylabel('Combined Score', fontsize=9)
                ax3.set_title(f'{prop.upper()} vs Combined Score{title_suffix}', fontsize=10)
                ax3.grid(True, alpha=0.3)
                # Add correlation
                valid = df_plot.loc[mask, [prop, 'combined_score']].dropna()
                if len(valid) > 1:
                    corr = valid.corr().iloc[0, 1]
                    color = 'green' if abs(corr) < 0.3 else 'orange' if abs(corr) < 0.5 else 'red'
                    ax3.text(0.05, 0.95, f'r = {corr:.3f}', 
                            transform=ax3.transAxes, va='top', fontsize=9,
                            bbox=dict(boxstyle='round', facecolor=color, alpha=0.3))
                # Add arrow showing better direction
                if i == 0:
                    ax3.annotate('Better', xy=(0.95, 0.05), xycoords='axes fraction',
                                fontsize=9, color='green', ha='right', va='bottom',
                                arrowprops=dict(arrowstyle='<-', color='green', lw=1.5))
        else:
            ax3.text(0.5, 0.5, 'No Combined Score', 
                    transform=ax3.transAxes, ha='center', va='center')
            ax3.set_xlabel(prop.replace('_', ' ').title(), fontsize=9)
    
    plt.tight_layout()
    plot_path = os.path.join(out_dir, 'property_correlations.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Saved correlation plots to: {plot_path}")
    plt.close()


def apply_filters(df: pd.DataFrame, filters: dict) -> pd.DataFrame:
    """Apply filters to dataframe based on molecular properties."""
    df_filtered = df.copy()
    
    for prop, (min_val, max_val) in filters.items():
        if prop not in df_filtered.columns:
            print(f"WARNING: Property '{prop}' not found in dataframe")
            continue
        
        if min_val is not None:
            df_filtered = df_filtered[df_filtered[prop] >= min_val]
        if max_val is not None:
            df_filtered = df_filtered[df_filtered[prop] <= max_val]
    
    return df_filtered


def print_statistics(df: pd.DataFrame, title: str = "Statistics"):
    """Print summary statistics."""
    print(f"\n{'='*60}")
    print(f"{title}")
    print(f"{'='*60}")
    print(f"Total compounds: {len(df)}")
    print(f"Compounds with predictions: {df['boltz_prob_mean'].notna().sum()}")
    print()
    
    if df['boltz_prob_mean'].notna().sum() > 0:
        print("Boltz Probability Mean:")
        print(f"  Mean:   {df['boltz_prob_mean'].mean():.4f}")
        print(f"  Median: {df['boltz_prob_mean'].median():.4f}")
        print(f"  Std:    {df['boltz_prob_mean'].std():.4f}")
        print(f"  Min:    {df['boltz_prob_mean'].min():.4f}")
        print(f"  Max:    {df['boltz_prob_mean'].max():.4f}")
        print()
        
        print("Boltz Pred Value Mean:")
        print(f"  Mean:   {df['boltz_pred_value_mean'].mean():.4f}")
        print(f"  Median: {df['boltz_pred_value_mean'].median():.4f}")
        print(f"  Std:    {df['boltz_pred_value_mean'].std():.4f}")
        print(f"  Min:    {df['boltz_pred_value_mean'].min():.4f}")
        print(f"  Max:    {df['boltz_pred_value_mean'].max():.4f}")
        print()
        
        if 'combined_score' in df.columns:
            # Filter out flagged compounds (score=999) for stats
            valid_scores = df[df['combined_score'] < 900]['combined_score']
            if len(valid_scores) > 0:
                print("Combined Score (lower = better, based on pKd + confidence):")
                print(f"  Mean:   {valid_scores.mean():.4f}")
                print(f"  Median: {valid_scores.median():.4f}")
                print(f"  Std:    {valid_scores.std():.4f}")
                print(f"  Min:    {valid_scores.min():.4f} (best)")
                print(f"  Max:    {valid_scores.max():.4f}")
            if 'passes_confidence' in df.columns:
                print(f"  Compounds passing confidence QC: {df['passes_confidence'].sum()}/{len(df)}")
            print()
    
    # Molecular properties statistics
    print("Molecular Properties:")
    for prop in ['mw', 'heavy_atoms', 'logp', 'hbd', 'hba', 'rotatable_bonds']:
        if prop in df.columns:
            print(f"  {prop}: {df[prop].mean():.2f} ± {df[prop].std():.2f} "
                  f"(range: {df[prop].min():.2f} - {df[prop].max():.2f})")


def main():
    ap = argparse.ArgumentParser(
        description="Analyze ligand predictions and create rankings",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis (ranks by combined_score: confidence QC + pKd)
  python analyze_ligand_predictions.py --input data.csv --out-dir results/

  # With filters
  python analyze_ligand_predictions.py --input data.csv --out-dir results/ \\
    --min-heavy-atoms 15 --max-heavy-atoms 30 --min-prob 0.5

  # Adjust confidence threshold (default 0.5)
  python analyze_ligand_predictions.py --input data.csv --out-dir results/ \\
    --confidence-threshold 0.6

  # Rank by single metric instead of combined score
  python analyze_ligand_predictions.py --input data.csv --out-dir results/ \\
    --rank-by boltz_pred_value_mean
        """
    )
    ap.add_argument("--input", required=True, 
                    help="Input CSV with Boltz predictions")
    ap.add_argument("--out-dir", required=True, 
                    help="Output directory for analysis results")
    
    # Ranking options
    ap.add_argument("--rank-by", default="combined_score",
                    choices=["combined_score", "boltz_prob_mean", "boltz_pred_value_mean"],
                    help="Column to rank by (default: combined_score)")
    ap.add_argument("--ascending", action="store_true",
                    help="Rank in ascending order (default: True for combined_score, False for others)")
    ap.add_argument("--top-n", type=int, default=100,
                    help="Number of top compounds to save separately (default: 100)")
    ap.add_argument("--confidence-threshold", type=float, default=0.5,
                    help="Minimum confidence (boltz_prob_mean) to pass QC (default: 0.5)")
    
    # Filter options
    ap.add_argument("--min-heavy-atoms", type=int,
                    help="Minimum number of heavy atoms")
    ap.add_argument("--max-heavy-atoms", type=int,
                    help="Maximum number of heavy atoms")
    ap.add_argument("--min-mw", type=float,
                    help="Minimum molecular weight")
    ap.add_argument("--max-mw", type=float,
                    help="Maximum molecular weight")
    ap.add_argument("--min-logp", type=float,
                    help="Minimum logP")
    ap.add_argument("--max-logp", type=float,
                    help="Maximum logP")
    ap.add_argument("--min-hbd", type=int,
                    help="Minimum H-bond donors")
    ap.add_argument("--max-hbd", type=int,
                    help="Maximum H-bond donors")
    ap.add_argument("--min-hba", type=int,
                    help="Minimum H-bond acceptors")
    ap.add_argument("--max-hba", type=int,
                    help="Maximum H-bond acceptors")
    ap.add_argument("--min-rotatable", type=int,
                    help="Minimum rotatable bonds")
    ap.add_argument("--max-rotatable", type=int,
                    help="Maximum rotatable bonds")
    ap.add_argument("--min-prob", type=float,
                    help="Minimum Boltz probability mean")
    ap.add_argument("--max-prob", type=float,
                    help="Maximum Boltz probability mean")
    ap.add_argument("--min-pred-value", type=float,
                    help="Minimum Boltz pred value mean")
    ap.add_argument("--max-pred-value", type=float,
                    help="Maximum Boltz pred value mean")
    ap.add_argument("--no-pains", action="store_true",
                    help="Filter out PAINS compounds")
    
    args = ap.parse_args()
    
    # Load data
    print(f"Loading data from: {args.input}")
    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} compounds")
    
    # Print initial statistics
    print_statistics(df, "Initial Dataset Statistics")
    
    # Calculate combined score
    print("\nCalculating combined score (two-stage: confidence QC + pKd ranking)...")
    print(f"  Confidence threshold: {args.confidence_threshold}")
    df = calculate_combined_score(df, confidence_threshold=args.confidence_threshold)
    print(f"  Compounds passing confidence QC: {df['passes_confidence'].sum()}/{len(df)}")
    
    # Apply filters
    filters = {}
    if args.min_heavy_atoms is not None or args.max_heavy_atoms is not None:
        filters['heavy_atoms'] = (args.min_heavy_atoms, args.max_heavy_atoms)
    if args.min_mw is not None or args.max_mw is not None:
        filters['mw'] = (args.min_mw, args.max_mw)
    if args.min_logp is not None or args.max_logp is not None:
        filters['logp'] = (args.min_logp, args.max_logp)
    if args.min_hbd is not None or args.max_hbd is not None:
        filters['hbd'] = (args.min_hbd, args.max_hbd)
    if args.min_hba is not None or args.max_hba is not None:
        filters['hba'] = (args.min_hba, args.max_hba)
    if args.min_rotatable is not None or args.max_rotatable is not None:
        filters['rotatable_bonds'] = (args.min_rotatable, args.max_rotatable)
    if args.min_prob is not None or args.max_prob is not None:
        filters['boltz_prob_mean'] = (args.min_prob, args.max_prob)
    if args.min_pred_value is not None or args.max_pred_value is not None:
        filters['boltz_pred_value_mean'] = (args.min_pred_value, args.max_pred_value)
    
    if filters:
        print(f"\nApplying filters: {filters}")
        df = apply_filters(df, filters)
        print(f"After filtering: {len(df)} compounds")
    
    # Filter PAINS if requested
    if args.no_pains and 'pains' in df.columns:
        before = len(df)
        df = df[df['pains'] == False]
        print(f"Filtered out PAINS: {before} -> {len(df)} compounds")
    
    # Print filtered statistics
    if filters or args.no_pains:
        print_statistics(df, "Filtered Dataset Statistics")
    
    # Create output directory
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Determine sort order (lower combined_score is better, higher prob is better)
    if args.rank_by == "combined_score":
        sort_ascending = True  # Lower score = better
    else:
        sort_ascending = args.ascending if args.ascending else False
    
    # Create ranking
    print(f"\nCreating ranking by {args.rank_by} ({'ascending' if sort_ascending else 'descending'})...")
    df_ranked = create_ranking(df, sort_by=args.rank_by, ascending=sort_ascending)
    
    # Define column order for CSV output (most important columns first)
    column_order = [
        'rank',                          # Ranking position
        'molport_number',                # Compound identifier
        'smiles',                        # SMILES string
        'combined_score',                # Combined score
        'boltz_prob_mean',               # Confidence
        'boltz_pred_value_mean',         # Predicted affinity
        'passes_confidence',             # QC flag
        'mw',                            # Molecular weight
        'heavy_atoms',                   # Heavy atom count
        'num_atoms',                     # Total atom count
        'logp',                          # LogP
        'hbd',                           # H-bond donors
        'hba',                           # H-bond acceptors
        'rotatable_bonds',                # Rotatable bonds
        'pains',                         # PAINS flag
        'compound_id',                   # Compound ID
        'name',                          # Compound name
    ]
    
    # Reorder columns: first the defined order, then any remaining columns
    existing_cols = [c for c in column_order if c in df_ranked.columns]
    remaining_cols = [c for c in df_ranked.columns if c not in column_order]
    df_ranked = df_ranked[existing_cols + remaining_cols]
    
    # Save full ranked list
    ranked_csv = os.path.join(args.out_dir, "ligands_ranked_full.csv")
    df_ranked.to_csv(ranked_csv, index=False)
    print(f"Saved full ranked list to: {ranked_csv}")
    
    # Save top N
    if len(df_ranked) >= args.top_n:
        top_n_csv = os.path.join(args.out_dir, f"ligands_ranked_top{args.top_n}.csv")
        df_ranked.head(args.top_n).to_csv(top_n_csv, index=False)
        print(f"Saved top {args.top_n} compounds to: {top_n_csv}")
        
        # Print top 10
        print(f"\n{'='*60}")
        print(f"Top 10 Compounds (ranked by {args.rank_by})")
        print(f"{'='*60}")
        display_cols = ['rank', 'molport_number', 'combined_score', 'boltz_prob_mean', 
                       'boltz_pred_value_mean', 'heavy_atoms', 'mw', 'logp']
        display_cols = [c for c in display_cols if c in df_ranked.columns]
        print(df_ranked[display_cols].head(10).to_string(index=False))
    
    # Create plots
    print("\nGenerating plots...")
    plot_distributions(df_ranked, args.out_dir)
    plot_property_correlations(df_ranked, args.out_dir)
    
    # Save summary statistics
    summary_path = os.path.join(args.out_dir, "analysis_summary.txt")
    with open(summary_path, 'w') as f:
        f.write("="*60 + "\n")
        f.write("Ligand Prediction Analysis Summary\n")
        f.write("="*60 + "\n\n")
        f.write(f"Input file: {args.input}\n")
        f.write(f"Ranking by: {args.rank_by} ({'ascending' if args.ascending else 'descending'})\n")
        f.write(f"Filters applied: {filters if filters else 'None'}\n")
        f.write(f"PAINS filtered: {args.no_pains}\n\n")
        f.write(f"Total compounds: {len(df_ranked)}\n")
        f.write(f"Compounds with predictions: {df_ranked['boltz_prob_mean'].notna().sum()}\n\n")
        
        if df_ranked['boltz_prob_mean'].notna().sum() > 0:
            f.write("Boltz Probability Mean Statistics:\n")
            f.write(f"  Mean:   {df_ranked['boltz_prob_mean'].mean():.4f}\n")
            f.write(f"  Median: {df_ranked['boltz_prob_mean'].median():.4f}\n")
            f.write(f"  Std:    {df_ranked['boltz_prob_mean'].std():.4f}\n")
            f.write(f"  Range:  [{df_ranked['boltz_prob_mean'].min():.4f}, {df_ranked['boltz_prob_mean'].max():.4f}]\n\n")
            
            f.write("Boltz Pred Value Mean Statistics:\n")
            f.write(f"  Mean:   {df_ranked['boltz_pred_value_mean'].mean():.4f}\n")
            f.write(f"  Median: {df_ranked['boltz_pred_value_mean'].median():.4f}\n")
            f.write(f"  Std:    {df_ranked['boltz_pred_value_mean'].std():.4f}\n")
            f.write(f"  Range:  [{df_ranked['boltz_pred_value_mean'].min():.4f}, {df_ranked['boltz_pred_value_mean'].max():.4f}]\n\n")
            
            if 'combined_score' in df_ranked.columns:
                f.write("Combined Score Statistics (high prob + negative pred is good):\n")
                f.write(f"  Mean:   {df_ranked['combined_score'].mean():.4f}\n")
                f.write(f"  Median: {df_ranked['combined_score'].median():.4f}\n")
                f.write(f"  Std:    {df_ranked['combined_score'].std():.4f}\n")
                f.write(f"  Range:  [{df_ranked['combined_score'].min():.4f}, {df_ranked['combined_score'].max():.4f}]\n")
    
    print(f"Saved analysis summary to: {summary_path}")
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
