#!/usr/bin/env python3
"""
Run DiffDock molecular docking for prepared ligands
"""
import argparse
import os
import sys
import subprocess
import shutil
import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import require_diffdock_home, require_micromamba


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    if not os.path.exists(config_path):
        print(f"[ERROR] Config file not found: {config_path}")
        sys.exit(1)
    
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def prepare_apo_protein(protein_pdb: str, output_path: str) -> bool:
    """Remove any existing ligands from protein structure."""
    try:
        from Bio import PDB
    except ImportError:
        print("[ERROR] Biopython not installed. Install with: pip install biopython")
        sys.exit(1)
    
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', protein_pdb)
        
        for model in structure:
            for chain in model:
                residues_to_remove = []
                for residue in chain:
                    if residue.id[0] != ' ':
                        residues_to_remove.append(residue.id)
                
                for res_id in residues_to_remove:
                    chain.detach_child(res_id)
        
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(output_path)
        
        return True
        
    except Exception as e:
        print(f"    Error preparing apo: {e}")
        return False


def run_diffdock(protein_path: str, ligand_sdf: str, output_dir: str,
                 diffdock_repo: str, micromamba_exe: str, micromamba_env: str,
                 samples: int = 40, timeout: int = 600,
                 inference_steps: int = 20, actual_steps: int = 18,
                 batch_size: int = 10) -> str:
    """Run DiffDock docking.

    inference_steps, actual_steps and batch_size used to be written into the
    command as literals while the config carried keys of the same name that
    nothing read. The values are unchanged (20/18/10, as in the thesis run);
    they now actually come from the config.
    """

    ensure_dir(output_dir)

    protein_abs = os.path.abspath(protein_path)
    ligand_abs = os.path.abspath(ligand_sdf)
    output_abs = os.path.abspath(output_dir)

    cmd = [
        micromamba_exe, "run", "-n", micromamba_env,
        "python", "-m", "inference",
        "--protein_path", protein_abs,
        "--ligand", ligand_abs,
        "--out_dir", output_abs,
        "--inference_steps", str(inference_steps),
        "--samples_per_complex", str(samples),
        "--batch_size", str(batch_size),
        "--actual_steps", str(actual_steps),
    ]
    
    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            cwd=diffdock_repo
        )
        
        complex_dir = os.path.join(output_dir, "complex_0")
        
        if os.path.exists(complex_dir):
            rank1_path = os.path.join(complex_dir, "rank1.sdf")
            if os.path.exists(rank1_path):
                return rank1_path
            
            sdf_files = [f for f in os.listdir(complex_dir) 
                       if f.startswith('rank1') and f.endswith('.sdf')]
            if sdf_files:
                return os.path.join(complex_dir, sdf_files[0])
        
        return None
        
    except subprocess.TimeoutExpired:
        print("    Timeout")
        return None
    except subprocess.CalledProcessError as e:
        stderr = e.stderr.decode() if e.stderr else ""
        print(f"    DiffDock error: {stderr[:100]}")
        return None
    except Exception as e:
        print(f"    Error: {e}")
        return None


def create_complex(protein_path: str, pose_sdf: str, output_path: str) -> bool:
    """Combine protein and docked pose into complex PDB."""
    try:
        from rdkit import Chem
    except ImportError:
        print("[ERROR] RDKit not installed")
        return False
    
    try:
        with open(protein_path, 'r') as f:
            protein_lines = f.readlines()
        
        suppl = Chem.SDMolSupplier(pose_sdf, removeHs=False)
        mol = next(suppl)
        
        if mol is None:
            return False
        
        ligand_pdb = Chem.MolToPDBBlock(mol)
        
        with open(output_path, 'w') as f:
            for line in protein_lines:
                if not line.startswith('END'):
                    f.write(line)
            
            for line in ligand_pdb.split('\n'):
                if line.startswith(('HETATM', 'ATOM')):
                    line = line[:21] + 'L' + line[22:]
                    f.write(line + '\n')
            
            f.write('END\n')
        
        return True
        
    except Exception as e:
        print(f"    Error creating complex: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Run DiffDock molecular docking"
    )
    parser.add_argument(
        "--config",
        help="Path to config YAML file (optional, overrides other args)"
    )
    parser.add_argument(
        "--input",
        help="Input CSV with sdf_path column (from script 05)"
    )
    parser.add_argument(
        "--protein",
        help="Protein PDB file (apo structure)"
    )
    parser.add_argument(
        "--output-dir",
        help="Output directory for DiffDock results"
    )
    parser.add_argument(
        "--diffdock-repo",
        help="Path to the DiffDock repository (default: $DIFFDOCK_HOME, then the config)"
    )
    parser.add_argument(
        "--micromamba",
        help="Micromamba executable"
    )
    parser.add_argument(
        "--env",
        help="Micromamba environment name"
    )
    parser.add_argument(
        "--samples",
        type=int,
        help="Samples per complex"
    )
    parser.add_argument(
        "--timeout",
        type=int,
        help="Timeout per ligand in seconds"
    )
    
    args = parser.parse_args()
    
    # Load config if provided
    if args.config:
        config = load_config(args.config)
        input_csv = args.input or config['input']['csv_path'].replace('.csv', '_with_sdf.csv')
        protein_pdb = args.protein or config['target']['protein_structure']
        output_dir = args.output_dir or os.path.join(config['output']['base_dir'], 'docking')
        diffdock_repo = args.diffdock_repo or config['diffdock'].get('repo_path')
        micromamba_exe = args.micromamba or config.get('micromamba', {}).get('executable', 'micromamba')
        micromamba_env = args.env or config.get('micromamba', {}).get('diffdock_env', 'diffdock')
        samples = args.samples or config['diffdock']['samples_per_complex']
        timeout = args.timeout or config['diffdock']['timeout']
        inference_steps = config['diffdock'].get('inference_steps', 20)
        actual_steps = config['diffdock'].get('actual_steps', 18)
        batch_size = config['diffdock'].get('batch_size', 10)
    else:
        if not args.input or not args.protein or not args.output_dir:
            parser.error("Either --config or --input, --protein, and --output-dir are required")
        input_csv = args.input
        protein_pdb = args.protein
        output_dir = args.output_dir
        diffdock_repo = args.diffdock_repo
        micromamba_exe = args.micromamba or "micromamba"
        micromamba_env = args.env or "diffdock"
        samples = args.samples or 40
        timeout = args.timeout or 600
        inference_steps = 20
        actual_steps = 18
        batch_size = 10

    # $DIFFDOCK_HOME wins over the config, and micromamba is looked up the same
    # way everywhere — see constraints/utils.py for why neither is a plain path.
    diffdock_repo = require_diffdock_home(diffdock_repo)
    micromamba_exe = require_micromamba(micromamba_exe)

    if not os.path.exists(input_csv):
        print(f"[ERROR] Input file not found: {input_csv}")
        sys.exit(1)
    
    if not os.path.exists(protein_pdb):
        print(f"[ERROR] Protein file not found: {protein_pdb}")
        sys.exit(1)
    
    df = pd.read_csv(input_csv)
    
    if 'sdf_path' not in df.columns:
        print(f"[ERROR] 'sdf_path' column not found. Run script 05 first.")
        sys.exit(1)
    
    poses_dir = os.path.join(output_dir, "poses")
    complex_dir = os.path.join(output_dir, "complexes")
    apo_dir = os.path.join(output_dir, "apo")
    
    ensure_dir(poses_dir)
    ensure_dir(complex_dir)
    ensure_dir(apo_dir)
    
    apo_path = os.path.join(apo_dir, "protein_apo.pdb")
    if not os.path.exists(apo_path):
        print("Preparing apo protein...")
        if not prepare_apo_protein(protein_pdb, apo_path):
            print("[ERROR] Failed to prepare apo protein")
            sys.exit(1)
        print(f"  ✓ Apo protein saved: {apo_path}\n")
    else:
        print(f"Using existing apo protein: {apo_path}\n")
    
    print(f"\n{'='*70}")
    print(f"  DIFFDOCK MOLECULAR DOCKING")
    print(f"{'='*70}\n")
    print(f"Input:    {input_csv}")
    print(f"Protein:  {apo_path}")
    print(f"Output:   {output_dir}")
    print(f"Samples:  {samples}")
    print(f"Timeout:  {timeout}s")
    print(f"Ligands:  {len(df)}\n")
    
    pose_paths = []
    complex_paths = []
    successful = 0
    failed = 0
    
    for idx, row in df.iterrows():
        ligand_id = row.get('molport_number', idx)
        sdf_path = row.get('sdf_path')
        
        if pd.isna(sdf_path) or not os.path.exists(str(sdf_path)):
            print(f"  [SKIP] {ligand_id}: No SDF file")
            pose_paths.append(None)
            complex_paths.append(None)
            failed += 1
            continue
        
        print(f"  [{idx+1}/{len(df)}] {ligand_id}...", end=" ")
        
        ligand_output = os.path.join(poses_dir, str(ligand_id))
        
        pose_path = run_diffdock(
            apo_path, sdf_path, ligand_output,
            diffdock_repo, micromamba_exe, micromamba_env,
            samples, timeout, inference_steps, actual_steps, batch_size
        )
        
        if pose_path and os.path.exists(pose_path):
            complex_path = os.path.join(complex_dir, f"{ligand_id}_complex.pdb")
            
            if create_complex(apo_path, pose_path, complex_path):
                print(f"✓")
                pose_paths.append(pose_path)
                complex_paths.append(complex_path)
                successful += 1
            else:
                print(f"✗ (complex creation failed)")
                pose_paths.append(pose_path)
                complex_paths.append(None)
                failed += 1
        else:
            print(f"✗ (docking failed)")
            pose_paths.append(None)
            complex_paths.append(None)
            failed += 1
    
    df['pose_path'] = pose_paths
    df['complex_path'] = complex_paths
    
    output_csv = input_csv.replace('.csv', '_with_diffdock.csv')
    df.to_csv(output_csv, index=False)
    
    print(f"\n{'='*70}")
    print(f"  SUMMARY")
    print(f"{'='*70}")
    print(f"Successful: {successful}/{len(df)}")
    print(f"Failed:     {failed}/{len(df)}")
    print(f"\nUpdated CSV saved: {output_csv}\n")


if __name__ == "__main__":
    main()
