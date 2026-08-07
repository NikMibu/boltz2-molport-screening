#!/usr/bin/env python3
"""
Quick Heavy Atom Comparison: Your hits vs. known CA inhibitors from ChEMBL
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import requests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

# ============================================================================
# 1. GET KNOWN CA2 INHIBITORS FROM CHEMBL
# ============================================================================

def get_ca2_inhibitors_from_chembl(max_compounds=500, activity_threshold_nm=10000, activity_types=['IC50', 'Ki']):
    """
    Fetch CA2 inhibitors from ChEMBL API
    
    CA2 Target ID: CHEMBL205
    Accepts both IC50 and Ki values
    """
    print(f"Fetching CA2 inhibitors from ChEMBL...")
    print(f"Parameters: Activity < {activity_threshold_nm} nM, types: {activity_types}, max {max_compounds} compounds")
    
    # ChEMBL API endpoint
    base_url = "https://www.ebi.ac.uk/chembl/api/data"
    
    # Collect compounds from all activity types (IC50, Ki, etc.)
    all_compounds = []
    
    for activity_type in activity_types:
        print(f"  Querying {activity_type}...")
        
        params = {
            'target_chembl_id': 'CHEMBL205',  # Carbonic anhydrase II
            'standard_type': activity_type,
            'standard_units': 'nM',
            'limit': max_compounds,
            'format': 'json'
        }
        
        url = f"{base_url}/activity.json"
        
        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            activities = data.get('activities', [])
            print(f"    Found {len(activities)} {activity_type} activities")
            
            # Extract molecule ChEMBL IDs and activity values
            for act in activities:
                try:
                    activity_value = float(act.get('standard_value', 999999))
                    relation = act.get('standard_relation', '=')
                    
                    # Filter by activity threshold
                    if activity_value <= activity_threshold_nm:
                        all_compounds.append({
                            'molecule_chembl_id': act.get('molecule_chembl_id'),
                            'activity_value_nm': activity_value,
                            'activity_type': activity_type,
                            'relation': relation
                        })
                except (ValueError, TypeError):
                    continue
        
        except Exception as e:
            print(f"    Error fetching {activity_type}: {e}")
            continue
    
    print(f"\nTotal activities collected: {len(all_compounds)}")
    
    # Get unique compounds (keep best activity per compound)
    df = pd.DataFrame(all_compounds)
    if len(df) > 0:
        df = df.sort_values('activity_value_nm').drop_duplicates('molecule_chembl_id', keep='first')
        print(f"Unique compounds (keeping best activity): {len(df)}")
    
    return df


def get_smiles_from_chembl(chembl_ids):
    """Get SMILES for ChEMBL IDs"""
    print(f"Fetching SMILES for {len(chembl_ids)} compounds...")
    
    base_url = "https://www.ebi.ac.uk/chembl/api/data"
    
    smiles_dict = {}
    batch_size = 50

    with requests.Session() as session:
        for i in range(0, len(chembl_ids), batch_size):
            batch = [cid for cid in chembl_ids[i:i+batch_size] if cid]
            if not batch:
                continue

            url = f"{base_url}/molecule.json"
            params = {
                "molecule_chembl_id__in": ",".join(batch),
                "limit": len(batch),
                "format": "json",
            }

            try:
                response = session.get(url, params=params, timeout=30)
                response.raise_for_status()
                data = response.json()

                for mol in data.get("molecules", []):
                    chembl_id = mol.get("molecule_chembl_id")
                    smiles = (mol.get("molecule_structures") or {}).get("canonical_smiles")
                    if chembl_id and smiles:
                        smiles_dict[chembl_id] = smiles

            except Exception as e:
                print(f"Error in batch {i//batch_size}: {e}")
                continue
    
    print(f"Retrieved SMILES for {len(smiles_dict)} compounds")
    return smiles_dict


def calculate_heavy_atoms(smiles):
    """Calculate heavy atoms from SMILES"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Descriptors.HeavyAtomCount(mol)
    except:
        return None


# ============================================================================
# 2. LOAD YOUR HITS
# ============================================================================

def load_your_hits(csv_path):
    """
    Load your Top-100 hits
    Expected columns: SMILES, combined_score (or similar)
    """
    print(f"\nLoading your hits from {csv_path}...")
    
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} compounds")

    smiles_col = None
    if "SMILES" in df.columns:
        smiles_col = "SMILES"
    elif "smiles" in df.columns:
        smiles_col = "smiles"
    else:
        raise KeyError("Input CSV must contain a 'SMILES' (or 'smiles') column")

    df['heavy_atoms'] = df[smiles_col].apply(calculate_heavy_atoms)
    df = df.dropna(subset=['heavy_atoms'])
    
    print(f"Valid compounds after heavy atom calculation: {len(df)}")
    
    return df


# ============================================================================
# 3. COMPARISON & VISUALIZATION
# ============================================================================

def compare_and_plot(your_hits_df, known_inhibitors_df, library_df=None, output_prefix='heavy_atom_comparison'):
    """
    Compare heavy atom distributions and create plots
    """
    print("\n" + "="*60)
    print("HEAVY ATOM COMPARISON")
    print("="*60)

    datasets = [
        ("Known CA2 Inhibitors", known_inhibitors_df['heavy_atoms'].dropna(), 'steelblue'),
        ("Boltz-2 Hits", your_hits_df['heavy_atoms'].dropna(), 'coral'),
    ]
    if library_df is not None:
        datasets.append(("Full Library", library_df['heavy_atoms'].dropna(), '#9FC5A8'))

    for label, values, _ in datasets:
        print(f"\n{label} (n={len(values)}):")
        print(f"  Mean ± Std: {values.mean():.1f} ± {values.std():.1f}")
        print(f"  Median: {values.median():.1f}")
        print(f"  Range: {values.min():.0f} - {values.max():.0f}")
        print(f"  25th-75th percentile: {values.quantile(0.25):.0f} - {values.quantile(0.75):.0f}")

    your_ha = datasets[1][1]
    for label, values, _ in datasets:
        if label == "Boltz-2 Hits":
            continue
        overlap_min = max(your_ha.min(), values.min())
        overlap_max = min(your_ha.max(), values.max())
        if overlap_min <= overlap_max:
            in_overlap = your_ha[(your_ha >= overlap_min) & (your_ha <= overlap_max)]
            print(f"\nOverlap with {label}: {overlap_min:.0f} - {overlap_max:.0f} heavy atoms")
            print(f"Hits in overlap: {len(in_overlap)}/{len(your_ha)} ({100*len(in_overlap)/len(your_ha):.1f}%)")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    ax1 = axes[0]
    all_values = pd.concat([values for _, values, _ in datasets], ignore_index=True)
    bins = range(int(all_values.min()), int(all_values.max()) + 2, 2)

    hist_plot_order = sorted(datasets, key=lambda item: 0 if item[0] == "Full Library" else 1)
    for label, values, color in hist_plot_order:
        is_library = label == "Full Library"
        ax1.hist(
            values,
            bins=bins,
            alpha=0.35 if is_library else 0.55,
            label=label,
            color=color,
            edgecolor='#6E9B79' if is_library else 'black',
            linewidth=0.6 if is_library else 0.8,
            zorder=1 if is_library else 2,
        )
        ax1.axvline(
            values.median(),
            color='#5E8A69' if is_library else color,
            linestyle=':' if is_library else '--',
            linewidth=1.8 if is_library else 2,
            alpha=1.0,
            label=f'{label} median: {values.median():.0f}',
            zorder=3,
        )

    ax1.set_xlabel('Heavy Atoms', fontsize=11)
    ax1.set_ylabel('Count', fontsize=11)
    ax1.set_title('Heavy Atom Distribution Comparison', fontsize=12, fontweight='bold')
    ax1.legend()
    ax1.grid(alpha=0.3)

    ax2 = axes[1]
    boxplot_order = [item for item in datasets if item[0] != "Full Library"]
    if any(label == "Full Library" for label, _, _ in datasets):
        boxplot_order.append(next(item for item in datasets if item[0] == "Full Library"))

    data_to_plot = [values for _, values, _ in boxplot_order]
    positions = list(range(1, len(boxplot_order) + 1))

    bp = ax2.boxplot(data_to_plot, positions=positions, widths=0.6,
                     patch_artist=True, showmeans=True,
                     meanprops=dict(marker='D', markerfacecolor='red', markersize=8))

    for patch, (label, _, color) in zip(bp['boxes'], boxplot_order):
        patch.set_facecolor(color)
        patch.set_alpha(0.5 if label == "Full Library" else 0.7)

    ax2.set_xticks(positions)
    ax2.set_xticklabels([label.replace(' ', '\n') for label, _, _ in boxplot_order])
    ax2.set_ylabel('Heavy Atoms', fontsize=11)
    ax2.set_title('Heavy Atom Distribution (Box Plot)', fontsize=12, fontweight='bold')
    ax2.grid(alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(f'{output_prefix}.png', dpi=300, bbox_inches='tight')
    print(f"\nPlot saved: {output_prefix}.png")

    summary = pd.DataFrame([
        {
            'Dataset': label,
            'N': len(values),
            'Mean': values.mean(),
            'Std': values.std(),
            'Median': values.median(),
            'Min': values.min(),
            'Max': values.max(),
            'Q25': values.quantile(0.25),
            'Q75': values.quantile(0.75),
        }
        for label, values, _ in datasets
    ])

    summary.to_csv(f'{output_prefix}_summary.csv', index=False)
    print(f"Summary saved: {output_prefix}_summary.csv")

    return summary


# ============================================================================
# MAIN
# ============================================================================

DEFAULT_KNOWN = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "reference", "known_ca2_inhibitors.csv",
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare heavy atom counts of a hit list against known CA2 inhibitors",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--hits", "--your-hits-csv", dest="hits", required=True,
                        help="CSV with the hit list; needs a 'SMILES' or 'smiles' column")
    parser.add_argument("--library-csv", default=None,
                        help="Optional full-library CSV for a third distribution")
    parser.add_argument("--known-csv", default=DEFAULT_KNOWN,
                        help="Reference set of known CA2 inhibitors. Pass --refetch to "
                             "rebuild it from ChEMBL instead of reading this file")
    parser.add_argument("--refetch", action="store_true",
                        help="Query the ChEMBL API instead of using --known-csv. "
                             "Needs network access and gives a different set as ChEMBL grows")
    parser.add_argument("--out-dir", default=".", help="Directory for the plot and the summary CSV")
    # 100000 nM = 100 uM. The threshold is deliberately permissive: it defines
    # "known binder", not "potent binder", so the novelty comparison is not
    # flattered by an artificially small reference set.
    parser.add_argument("--activity-threshold-nm", type=float, default=100000,
                        help="Keep ChEMBL activities at or below this value, in nM (only with --refetch)")
    parser.add_argument("--max-compounds", type=int, default=1000,
                        help="ChEMBL result limit per activity type (only with --refetch)")
    args = parser.parse_args()

    ACTIVITY_TYPES = ['IC50', 'Ki']
    os.makedirs(args.out_dir, exist_ok=True)

    print("="*60)
    print("HEAVY ATOM COMPARISON: hit list vs known CA2 inhibitors")
    print("="*60)

    # Step 1: reference set - from file (default) or freshly from ChEMBL
    if args.refetch:
        ca2_activities = get_ca2_inhibitors_from_chembl(
            max_compounds=args.max_compounds,
            activity_threshold_nm=args.activity_threshold_nm,
            activity_types=ACTIVITY_TYPES,
        )
        if len(ca2_activities) == 0:
            print("No CA2 inhibitors returned by ChEMBL. Exiting.")
            exit(1)

        chembl_ids = ca2_activities['molecule_chembl_id'].tolist()
        smiles_dict = get_smiles_from_chembl(chembl_ids)
        ca2_activities['smiles'] = ca2_activities['molecule_chembl_id'].map(smiles_dict)
        ca2_activities = ca2_activities.dropna(subset=['smiles'])
        ca2_activities['heavy_atoms'] = ca2_activities['smiles'].apply(calculate_heavy_atoms)
        ca2_activities = ca2_activities.dropna(subset=['heavy_atoms'])

        refreshed = os.path.join(args.out_dir, 'known_ca2_inhibitors.csv')
        ca2_activities.to_csv(refreshed, index=False)
        print(f"Saved refreshed reference set: {refreshed}")
    else:
        print(f"\nReading reference set from {args.known_csv}")
        ca2_activities = pd.read_csv(args.known_csv)
        if 'heavy_atoms' not in ca2_activities.columns:
            ca2_activities['heavy_atoms'] = ca2_activities['smiles'].apply(calculate_heavy_atoms)
        ca2_activities = ca2_activities.dropna(subset=['heavy_atoms'])

    print(f"Known CA2 inhibitors: {len(ca2_activities)} compounds")

    # Step 2: load the hit list
    try:
        your_hits = load_your_hits(args.hits)
    except FileNotFoundError:
        print(f"\nERROR: could not find {args.hits}")
        print(f"Working directory: {os.getcwd()}")
        print("The hit list is not shipped with this repository - pass your own CSV")
        print("with at least a 'SMILES' column via --hits.")
        exit(1)

    library_hits = None
    if args.library_csv:
        try:
            library_hits = load_your_hits(args.library_csv)
        except FileNotFoundError:
            print(f"\nWARNING: could not find library CSV: {args.library_csv}")

    # Step 3: compare and visualise
    prefix = os.path.join(args.out_dir, 'heavy_atom_comparison')
    summary = compare_and_plot(your_hits, ca2_activities,
                               library_df=library_hits, output_prefix=prefix)

    print("\nDone. Output files:")
    print(f"  {prefix}.png")
    print(f"  {prefix}_summary.csv")