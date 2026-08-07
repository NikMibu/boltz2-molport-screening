#!/usr/bin/env python3
"""Analyze Boltz predictions with constraints and compare rankings."""
import argparse
import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml


def safe_mean_from_prefix(d: dict, prefix: str) -> float:
    vals = [float(v) for k, v in d.items() if k.startswith(prefix) and isinstance(v, (int, float))]
    return float(np.mean(vals)) if vals else float("nan")


def find_affinity_json(base_dir: str, molport_id: str, constraint_type: str):
    if not os.path.exists(base_dir):
        return None
    
    for entry in os.scandir(base_dir):
        if not entry.is_dir():
            continue
        
        predictions_dir = os.path.join(entry.path, "predictions")
        if not os.path.exists(predictions_dir):
            predictions_dir = entry.path
        
        for pred_entry in os.scandir(predictions_dir):
            if not pred_entry.is_dir():
                continue
            
            if molport_id in pred_entry.name and constraint_type in pred_entry.name:
                aff_files = [f for f in os.listdir(pred_entry.path) 
                           if f.startswith("affinity_") and f.endswith(".json")]
                if aff_files:
                    try:
                        with open(os.path.join(pred_entry.path, aff_files[0])) as f:
                            return json.load(f)
                    except:
                        pass
    return None


def collect_predictions(results_dir: str, molport_ids: list, constraint_types: list):
    rows = []
    for molport_id in molport_ids:
        row = {'molport_number': molport_id}
        for ctype in constraint_types:
            affinity_dict = find_affinity_json(results_dir, molport_id, ctype)
            if affinity_dict:
                row[f'{ctype}_prob'] = safe_mean_from_prefix(affinity_dict, 'affinity_probability_binary')
                row[f'{ctype}_pred'] = safe_mean_from_prefix(affinity_dict, 'affinity_pred_value')
            else:
                row[f'{ctype}_prob'] = np.nan
                row[f'{ctype}_pred'] = np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def calc_combined_score(df: pd.DataFrame, prefix: str, conf_threshold: float = 0.5):
    prob_col = f'{prefix}_prob'
    pred_col = f'{prefix}_pred'
    
    if prob_col not in df.columns or pred_col not in df.columns:
        return df
    
    df[f'{prefix}_pass_qc'] = df[prob_col].fillna(0) >= conf_threshold
    df[f'{prefix}_score'] = np.where(
        df[f'{prefix}_pass_qc'],
        df[pred_col] - 0.1 * df[prob_col],
        999
    )
    return df


def plot_comparison(df: pd.DataFrame, constraint_types: list, out_dir: str):
    """Compact comparison plots."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. Violin plots for all constraint types
    ax1 = axes[0, 0]
    prob_data = []
    labels = []
    for ctype in constraint_types:
        if f'{ctype}_prob' in df.columns:
            data = df[f'{ctype}_prob'].dropna()
            if len(data) > 0:
                prob_data.append(data)
                labels.append(ctype.capitalize())
    
    if prob_data:
        ax1.violinplot(prob_data, showmeans=True, showmedians=True)
        ax1.set_xticks(range(1, len(labels) + 1))
        ax1.set_xticklabels(labels)
        ax1.set_ylabel('Binding Probability')
        ax1.set_title('Binding Probability by Constraint Type', fontweight='bold')
        ax1.grid(True, alpha=0.3, axis='y')
    
    # 2. Predicted affinity distributions
    ax2 = axes[0, 1]
    pred_data = []
    for ctype in constraint_types:
        if f'{ctype}_pred' in df.columns:
            data = df[f'{ctype}_pred'].dropna()
            if len(data) > 0:
                pred_data.append(data)
    
    if pred_data:
        ax2.violinplot(pred_data, showmeans=True, showmedians=True)
        ax2.set_xticks(range(1, len(labels) + 1))
        ax2.set_xticklabels(labels)
        ax2.set_ylabel('Predicted Affinity')
        ax2.set_title('Affinity by Constraint Type', fontweight='bold')
        ax2.axhline(0, color='red', linestyle='--', linewidth=1, alpha=0.5)
        ax2.grid(True, alpha=0.3, axis='y')
    
    # 3. Default vs Pocket correlation
    ax3 = axes[1, 0]
    if 'default_pred' in df.columns and 'pocket_pred' in df.columns:
        valid = df[['default_pred', 'pocket_pred']].dropna()
        if len(valid) > 0:
            ax3.scatter(valid['default_pred'], valid['pocket_pred'], 
                       alpha=0.5, s=20, edgecolors='black', linewidth=0.5)
            ax3.set_xlabel('Default Affinity')
            ax3.set_ylabel('Pocket Affinity')
            ax3.set_title('Default vs Pocket Correlation', fontweight='bold')
            
            lims = [min(ax3.get_xlim()[0], ax3.get_ylim()[0]),
                    max(ax3.get_xlim()[1], ax3.get_ylim()[1])]
            ax3.plot(lims, lims, 'r--', alpha=0.5)
            ax3.grid(True, alpha=0.3)
            
            if len(valid) > 1:
                corr = valid.corr().iloc[0, 1]
                ax3.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax3.transAxes, 
                        va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # 4. Rank shift default vs pocket
    ax4 = axes[1, 1]
    if 'default_score' in df.columns and 'pocket_score' in df.columns:
        df_both = df[(df['default_score'] < 900) & (df['pocket_score'] < 900)].copy()
        if len(df_both) > 0:
            df_both = df_both.sort_values('default_score')
            df_both['default_rank'] = range(1, len(df_both) + 1)
            df_both = df_both.sort_values('pocket_score')
            df_both['pocket_rank'] = range(1, len(df_both) + 1)
            df_both['shift'] = df_both['default_rank'] - df_both['pocket_rank']
            
            ax4.hist(df_both['shift'], bins=30, edgecolor='black', alpha=0.7, color='steelblue')
            ax4.axvline(0, color='red', linestyle='--', linewidth=2)
            ax4.set_xlabel('Rank Shift (Default - Pocket)')
            ax4.set_ylabel('Frequency')
            ax4.set_title('Ranking Shifts', fontweight='bold')
            ax4.grid(True, alpha=0.3, axis='y')
            
            mean_shift = df_both['shift'].mean()
            ax4.text(0.05, 0.95, f'Mean: {mean_shift:.1f}', transform=ax4.transAxes, 
                    va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'constraint_comparison.png'), dpi=300, bbox_inches='tight')
    print(f"Saved comparison plot")
    plt.close()


def plot_old_vs_new(df: pd.DataFrame, old_col: str, new_col: str, out_dir: str):
    """Compact old vs new ranking comparison."""
    df_both = df[[old_col, new_col, 'molport_number']].dropna()
    if len(df_both) == 0:
        return
    
    df_both = df_both.copy()
    df_both['shift'] = df_both[old_col] - df_both[new_col]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # 1. Scatter old vs new
    ax1 = axes[0]
    scatter = ax1.scatter(df_both[old_col], df_both[new_col],
                         alpha=0.6, s=40, c=df_both['shift'], 
                         cmap='RdYlGn', edgecolors='black', linewidth=0.5,
                         vmin=-df_both['shift'].abs().max(),
                         vmax=df_both['shift'].abs().max())
    ax1.set_xlabel('Old Rank (without constraints)', fontweight='bold')
    ax1.set_ylabel('New Rank (with constraints)', fontweight='bold')
    ax1.set_title('Old vs New Ranking', fontweight='bold')
    
    max_rank = max(df_both[old_col].max(), df_both[new_col].max())
    ax1.plot([1, max_rank], [1, max_rank], 'r--', alpha=0.5, linewidth=2)
    ax1.grid(True, alpha=0.3)
    
    plt.colorbar(scatter, ax=ax1, label='Rank Shift (positive = improved)')
    
    if len(df_both) > 1:
        corr = df_both[[old_col, new_col]].corr().iloc[0, 1]
        ax1.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax1.transAxes, 
                va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    
    # 2. Histogram
    ax2 = axes[1]
    ax2.hist(df_both['shift'], bins=30, edgecolor='black', alpha=0.7, color='steelblue')
    ax2.axvline(0, color='red', linestyle='--', linewidth=2)
    ax2.set_xlabel('Rank Shift (Old - New)', fontweight='bold')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Distribution of Rank Shifts', fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    mean_shift = df_both['shift'].mean()
    median_shift = df_both['shift'].median()
    improved = (df_both['shift'] > 0).sum()
    worsened = (df_both['shift'] < 0).sum()
    
    stats_text = (f'Mean: {mean_shift:.1f}\n'
                 f'Median: {median_shift:.1f}\n'
                 f'Improved: {improved} ({100*improved/len(df_both):.0f}%)\n'
                 f'Worsened: {worsened} ({100*worsened/len(df_both):.0f}%)')
    ax2.text(0.95, 0.95, stats_text, transform=ax2.transAxes, 
            va='top', ha='right', bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
    
    # 3. Top 20 overlap
    ax3 = axes[2]
    top_n = min(20, len(df_both))
    old_top = set(df_both.nsmallest(top_n, old_col)['molport_number'].tolist())
    new_top = set(df_both.nsmallest(top_n, new_col)['molport_number'].tolist())
    
    overlap = len(old_top & new_top)
    only_old = len(old_top - new_top)
    only_new = len(new_top - old_top)
    
    bars = ax3.bar(['In both', 'Only old', 'Only new'], 
                   [overlap, only_old, only_new],
                   color=['#90EE90', '#FFB6C6', '#87CEEB'],
                   edgecolor='black', linewidth=1.5, alpha=0.8)
    ax3.set_ylabel('Number of Compounds')
    ax3.set_title(f'Top {top_n} Overlap', fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='y')
    
    for bar, count in zip(bars, [overlap, only_old, only_new]):
        ax3.text(bar.get_x() + bar.get_width()/2., bar.get_height(),
                f'{count}', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'old_vs_new_comparison.png'), dpi=300, bbox_inches='tight')
    print(f"Saved old vs new comparison")
    plt.close()
    
    # Print stats
    print(f"\n{'='*60}")
    print(f"OLD vs NEW RANKING")
    print(f"{'='*60}")
    print(f"Total compared: {len(df_both)}")
    print(f"Mean shift: {mean_shift:.1f} (positive = improved)")
    print(f"Improved: {improved} ({100*improved/len(df_both):.0f}%)")
    print(f"Worsened: {worsened} ({100*worsened/len(df_both):.0f}%)")
    print(f"Top {top_n} overlap: {overlap}/{top_n}")
    if len(df_both) > 1:
        print(f"Correlation: r = {corr:.3f}")


def print_stats(df: pd.DataFrame, prefix: str):
    prob_col = f'{prefix}_prob'
    pred_col = f'{prefix}_pred'
    score_col = f'{prefix}_score'
    
    if prob_col not in df.columns:
        return
    
    n = df[prob_col].notna().sum()
    if n == 0:
        return
    
    print(f"\n{prefix.upper()}: n={n}")
    print(f"  Binding Probability: {df[prob_col].mean():.3f} ± {df[prob_col].std():.3f}")
    print(f"  Affinity:   {df[pred_col].mean():.3f} ± {df[pred_col].std():.3f}")
    
    if score_col in df.columns:
        valid = df[df[score_col] < 900][score_col]
        if len(valid) > 0:
            print(f"  Score:      {valid.mean():.3f} ± {valid.std():.3f} (lower=better)")
            print(f"  Pass QC:    {len(valid)}/{n}")


def load_config(config_path: str) -> dict:
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser(description="Analyze Boltz predictions with constraints")
    parser.add_argument("--config", help="Config YAML")
    parser.add_argument("--input", help="Input CSV with molport_number")
    parser.add_argument("--results-dir", help="Boltz results directory")
    parser.add_argument("--out-dir", help="Output directory")
    parser.add_argument("--constraint-types", nargs='+', 
                       default=['default', 'pocket', 'contact'],
                       help="Constraint types")
    parser.add_argument("--rank-by", default="pocket", help="Rank by this constraint type")
    parser.add_argument("--confidence-threshold", type=float, default=0.5)
    parser.add_argument("--top-n", type=int, default=100)
    
    args = parser.parse_args()
    
    # Load from config or args
    if args.config:
        config = load_config(args.config)
        input_csv = args.input or config['input']['csv_path']
        results_dir = args.results_dir or os.path.join(config['output']['base_dir'], 'boltz_output')
        out_dir = args.out_dir or os.path.join(config['output']['base_dir'], 'analysis')
    else:
        if not args.input or not args.results_dir or not args.out_dir:
            parser.error("Need --config or --input, --results-dir, --out-dir")
        input_csv = args.input
        results_dir = args.results_dir
        out_dir = args.out_dir
    
    # Load input
    print(f"Loading: {input_csv}")
    df_input = pd.read_csv(input_csv)
    
    if 'molport_number' not in df_input.columns:
        print(f"ERROR: molport_number not found")
        sys.exit(1)
    
    print(f"Found {len(df_input)} compounds")
    
    # Collect predictions
    print(f"\nCollecting predictions from: {results_dir}")
    df_pred = collect_predictions(results_dir, df_input['molport_number'].tolist(), args.constraint_types)
    df = df_input.merge(df_pred, on='molport_number', how='left')
    
    # Calculate scores
    print(f"\nCalculating scores (threshold: {args.confidence_threshold})...")
    for ctype in args.constraint_types:
        df = calc_combined_score(df, ctype, args.confidence_threshold)
    
    os.makedirs(out_dir, exist_ok=True)
    
    # Print stats
    for ctype in args.constraint_types:
        print_stats(df, ctype)
    
    # Create plots
    print("\nGenerating plots...")
    plot_comparison(df, args.constraint_types, out_dir)
    
    # Create ranking
    score_col = f'{args.rank_by}_score'
    if score_col not in df.columns:
        print(f"ERROR: {score_col} not found")
        sys.exit(1)
    
    print(f"\nRanking by {args.rank_by}...")
    
    # Check for old ranking
    old_rank_col = None
    if 'rank' in df.columns:
        if 'combined_score' in df.columns or 'boltz_prob_mean' in df.columns:
            old_rank_col = 'rank'
            print(f"Found old ranking")
    
    # Sort and rank
    df = df.sort_values(score_col, na_position='last')
    if old_rank_col:
        df = df.rename(columns={'rank': 'rank_old'})
        old_rank_col = 'rank_old'
    df['rank'] = range(1, len(df) + 1)
    
    # Reorder columns
    cols = ['rank', 'molport_number', 'smiles']
    if old_rank_col:
        cols.append(old_rank_col)
    
    for ctype in args.constraint_types:
        cols.extend([f'{ctype}_score', f'{ctype}_prob', f'{ctype}_pred'])
    
    cols = [c for c in cols if c in df.columns]
    remaining = [c for c in df.columns if c not in cols]
    df = df[cols + remaining]
    
    # Save
    out_csv = os.path.join(out_dir, f"ranked_by_{args.rank_by}_full.csv")
    df.to_csv(out_csv, index=False)
    print(f"Saved: {out_csv}")
    
    if len(df) >= args.top_n:
        top_csv = os.path.join(out_dir, f"ranked_by_{args.rank_by}_top{args.top_n}.csv")
        df.head(args.top_n).to_csv(top_csv, index=False)
        print(f"Saved top {args.top_n}: {top_csv}")
    
    # Print top 10
    print(f"\n{'='*60}")
    print(f"Top 10 ({args.rank_by})")
    print(f"{'='*60}")
    display_cols = ['rank', 'molport_number', score_col, f'{args.rank_by}_prob', f'{args.rank_by}_pred']
    if old_rank_col:
        display_cols.insert(2, old_rank_col)
    display_cols = [c for c in display_cols if c in df.columns]
    print(df[display_cols].head(10).to_string(index=False))
    
    # Compare old vs new
    if old_rank_col:
        print(f"\nComparing old vs new ranking...")
        plot_old_vs_new(df, old_rank_col, 'rank', out_dir)
    
    print("\nDone!")


if __name__ == "__main__":
    main()
