#!/usr/bin/env python3
"""
Convert SMILES from top-ranked ligands to 3D SDF files for DiffDock
"""
import argparse
import os
import sys
import pandas as pd
import yaml


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    if not os.path.exists(config_path):
        print(f"[ERROR] Config file not found: {config_path}")
        sys.exit(1)
    
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def smiles_to_sdf(smiles: str, output_path: str) -> bool:
    """Convert SMILES to 3D SDF file using RDKit."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        print("[ERROR] RDKit not installed. Install with: pip install rdkit")
        sys.exit(1)
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        mol = Chem.AddHs(mol)
        
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if result != 0:
            return False
        
        try:
            AllChem.UFFOptimizeMolecule(mol)
        except AttributeError:
            AllChem.UFFOptimize(mol)
        
        writer = Chem.SDWriter(output_path)
        writer.write(mol)
        writer.close()
        
        return True
        
    except Exception as e:
        print(f"    Error: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Convert SMILES to 3D SDF files for DiffDock"
    )
    parser.add_argument(
        "--config",
        help="Path to config YAML file (optional, overrides other args)"
    )
    parser.add_argument(
        "--input",
        help="Input CSV file with smiles column"
    )
    parser.add_argument(
        "--output-dir",
        help="Output directory for SDF files"
    )
    parser.add_argument(
        "--id-column",
        help="Column name for ligand ID (default from config or 'molport_number')"
    )
    
    args = parser.parse_args()
    
    # Load config if provided
    if args.config:
        config = load_config(args.config)
        input_csv = args.input or config['input']['csv_path']
        output_dir = args.output_dir or os.path.join(config['output']['base_dir'], 'sdf')
        id_column = args.id_column or config['input']['id_column']
    else:
        if not args.input or not args.output_dir:
            parser.error("Either --config or both --input and --output-dir are required")
        input_csv = args.input
        output_dir = args.output_dir
        id_column = args.id_column or "molport_number"
    
    if not os.path.exists(input_csv):
        print(f"[ERROR] Input file not found: {input_csv}")
        sys.exit(1)
    
    df = pd.read_csv(input_csv)
    
    if 'smiles' not in df.columns:
        print(f"[ERROR] 'smiles' column not found in input CSV")
        sys.exit(1)
    
    if id_column not in df.columns:
        print(f"[ERROR] '{id_column}' column not found in input CSV")
        sys.exit(1)
    
    ensure_dir(output_dir)
    
    print(f"\n{'='*70}")
    print(f"  LIGAND PREPARATION (SMILES → SDF)")
    print(f"{'='*70}\n")
    print(f"Input:  {input_csv}")
    print(f"Output: {output_dir}")
    print(f"ID Col: {id_column}")
    print(f"Ligands: {len(df)}\n")
    
    sdf_paths = []
    successful = 0
    failed = 0
    
    for idx, row in df.iterrows():
        ligand_id = row[id_column]
        smiles = row['smiles']
        
        if pd.isna(smiles) or not smiles:
            print(f"  [SKIP] {ligand_id}: No SMILES")
            sdf_paths.append(None)
            failed += 1
            continue
        
        sdf_path = os.path.join(output_dir, f"{ligand_id}.sdf")
        
        if smiles_to_sdf(smiles, sdf_path):
            print(f"  ✓ {ligand_id}")
            sdf_paths.append(sdf_path)
            successful += 1
        else:
            print(f"  ✗ {ligand_id}: Failed")
            sdf_paths.append(None)
            failed += 1
    
    df['sdf_path'] = sdf_paths
    
    output_csv = input_csv.replace('.csv', '_with_sdf.csv')
    df.to_csv(output_csv, index=False)
    
    print(f"\n{'='*70}")
    print(f"  SUMMARY")
    print(f"{'='*70}")
    print(f"Successful: {successful}/{len(df)}")
    print(f"Failed:     {failed}/{len(df)}")
    print(f"\nUpdated CSV saved: {output_csv}\n")


if __name__ == "__main__":
    main()
