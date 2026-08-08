#!/usr/bin/env python3
"""
Analyze DiffDock complexes with PLIP and generate Boltz-2 YAMLs with constraints
"""
import argparse
import os
import sys
import subprocess
import xml.etree.ElementTree as ET
import re
import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import require_micromamba


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def load_config(config_path: str) -> dict:
    """Load configuration from YAML file."""
    if not os.path.exists(config_path):
        print(f"[ERROR] Config file not found: {config_path}")
        sys.exit(1)
    
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def yaml_quote(s: str) -> str:
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def run_plip(complex_pdb: str, output_dir: str, micromamba_exe: str, 
             micromamba_env: str, timeout: int = 120) -> str:
    """Run PLIP analysis on complex."""
    
    ensure_dir(output_dir)
    
    cmd = [
        micromamba_exe, "run", "-n", micromamba_env,
        "plip",
        "-f", os.path.abspath(complex_pdb),
        "-x",
        "-o", os.path.abspath(output_dir)
    ]
    
    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout
        )
        
        xml_path1 = os.path.join(output_dir, "report.xml")
        xml_path2 = os.path.join(output_dir, f"{os.path.basename(complex_pdb).replace('.pdb', '')}_report.xml")
        
        if os.path.exists(xml_path1):
            return xml_path1
        elif os.path.exists(xml_path2):
            return xml_path2
        else:
            return None
            
    except subprocess.TimeoutExpired:
        return None
    except subprocess.CalledProcessError:
        return None
    except Exception:
        return None


def parse_residue_id(residue_text: str):
    """Parse residue identifier like '199A' or '199:A'."""
    if not residue_text:
        return None, None
    
    residue_text = residue_text.strip()
    
    if ":" in residue_text:
        parts = residue_text.split(":")
        try:
            return parts[1], int(parts[0])
        except (ValueError, IndexError):
            return None, None
    
    match = re.match(r"(-?\d+)([A-Za-z]+)$", residue_text)
    if match:
        resnr = int(match.group(1))
        chain = match.group(2)
        return chain, resnr
    
    return None, None


def extract_plip_contacts(xml_path: str, ignore_metal: bool = True):
    """Extract pocket contacts from PLIP XML."""
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
    except Exception:
        return []
    
    all_contacts = []
    
    for site in root.findall(".//bindingsite"):
        hetid_elem = site.find(".//hetid")
        if hetid_elem is None:
            continue
        
        hetid = hetid_elem.text.strip() if hetid_elem.text else ""
        
        if hetid not in ['UNL', 'LIG']:
            continue
        
        if ignore_metal:
            ligtype_elem = site.find(".//ligtype")
            ligtype = ligtype_elem.text.strip().upper() if ligtype_elem is not None and ligtype_elem.text else ""
            
            if ligtype in {"ION", "METAL", "COFACTOR"}:
                continue
        
        bs_residues = site.find("bs_residues")
        if bs_residues is not None:
            for bs_res in bs_residues.findall("bs_residue"):
                contact_flag = bs_res.get("contact", "False").lower() == "true"
                if not contact_flag:
                    continue
                
                chain, resnr = parse_residue_id(bs_res.text)
                if chain and resnr:
                    all_contacts.append([chain, str(resnr)])
        
        interaction_types = [
            ".//hydrogen_bond",
            ".//hydrophobic_interaction",
            ".//pi_stacking",
            ".//salt_bridge",
            ".//halogen_bond"
        ]
        
        for interaction_type in interaction_types:
            for interaction in site.findall(interaction_type):
                chain_elem = interaction.find("reschain")
                resnr_elem = interaction.find("resnr")
                
                if chain_elem is not None and resnr_elem is not None:
                    chain = chain_elem.text.strip()
                    resnr = resnr_elem.text.strip()
                    all_contacts.append([chain, resnr])
    
    unique_contacts = []
    seen = set()
    for contact in all_contacts:
        key = tuple(contact)
        if key not in seen:
            seen.add(key)
            unique_contacts.append(contact)
    
    return unique_contacts


def generate_yaml(sequence: str, ligand_id: str, smiles: str, msa_path: str,
                  contacts: list, max_distance: float = 10.0, 
                  constraint_type: str = 'pocket') -> str:
    """Generate Boltz-2 YAML with constraints."""
    
    lines = []
    lines.append("version: 1")
    lines.append("")
    lines.append("sequences:")
    lines.append("  - protein:")
    lines.append('      id: ["A"]')
    lines.append(f"      sequence: {yaml_quote(sequence)}")
    
    if msa_path and os.path.exists(msa_path):
        lines.append(f"      msa: {yaml_quote(msa_path)}")
    else:
        lines.append("      msa: empty")
    
    lines.append("")
    lines.append("  - ligand:")
    lines.append(f"      id: [{yaml_quote(ligand_id)}]")
    lines.append(f"      smiles: {yaml_quote(smiles)}")
    lines.append("")
    lines.append("properties:")
    lines.append("  - affinity:")
    lines.append(f"      binder: {yaml_quote(ligand_id)}")
    
    if contacts and constraint_type == 'pocket':
        lines.append("")
        lines.append("constraints:")
        lines.append("  - pocket:")
        lines.append(f"      binder: {yaml_quote(ligand_id)}")
        
        contacts_str = "["
        for i, contact in enumerate(contacts):
            if i > 0:
                contacts_str += ", "
            chain = contact[0]
            residue = contact[1]
            try:
                residue_int = int(str(residue))
                contacts_str += f'["{chain}", {residue_int}]'
            except ValueError:
                contacts_str += f'["{chain}", "{residue}"]'
        contacts_str += "]"
        
        lines.append(f"      contacts: {contacts_str}")
        lines.append(f"      max_distance: {max_distance}")
        lines.append("      force: false")
    
    elif contacts and constraint_type == 'contact':
        lines.append("")
        lines.append("constraints:")
        
        for contact in contacts:
            chain = contact[0]
            residue = contact[1]
            
            lines.append("  - contact:")
            lines.append(f'      token1: [{yaml_quote(ligand_id)}]')
            
            try:
                residue_int = int(str(residue))
                lines.append(f'      token2: ["{chain}", {residue_int}]')
            except ValueError:
                lines.append(f'      token2: ["{chain}", "{residue}"]')
            
            lines.append(f"      max_distance: {max_distance}")
            lines.append("      force: false")
    
    lines.append("")
    return "\n".join(lines)


def read_fasta_sequence(fasta_path: str) -> str:
    if not os.path.isfile(fasta_path):
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")
    seq = []
    with open(fasta_path, "r", encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith(">"):
                continue
            seq.append(s)
    if not seq:
        raise ValueError(f"No sequence in FASTA: {fasta_path}")
    return "".join(seq).upper()


def main():
    parser = argparse.ArgumentParser(
        description="Generate Boltz-2 YAMLs with PLIP-derived pocket constraints",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Two distances are involved and they are easy to confuse:

  PLIP binding-site radius   which residues count as contacts. This script
                             passes no distance argument, so PLIP uses its own
                             default (BS_DIST, 7.5 A in 2.3.x).
  --max-distance             how far Boltz-2 may place the binder from those
                             contacts. Written into the YAML. Default 6.0.

Only the second one is set here.
"""
    )
    parser.add_argument(
        "--config",
        help="Path to config YAML file (optional, overrides other args)"
    )
    parser.add_argument(
        "--input",
        help="Input CSV with complex_path column (from script 06)"
    )
    parser.add_argument(
        "--fasta",
        help="Protein sequence FASTA file"
    )
    parser.add_argument(
        "--output-dir",
        help="Output directory for YAML files"
    )
    parser.add_argument(
        "--plip-output",
        help="Output directory for PLIP reports"
    )
    parser.add_argument(
        "--msa",
        help="MSA file (optional)"
    )
    parser.add_argument(
        "--micromamba",
        help="Micromamba executable"
    )
    parser.add_argument(
        "--env",
        help="Micromamba environment for PLIP"
    )
    parser.add_argument(
        "--max-distance",
        type=float,
        help="max_distance of the Boltz-2 pocket constraint, in Angstrom. "
             "NOT a PLIP setting: PLIP is invoked without a distance argument "
             "and detects contacts at its own default radius. (default: 6.0)"
    )
    parser.add_argument(
        "--constraint-type",
        choices=['pocket', 'contact', 'both'],
        help="Constraint type"
    )
    parser.add_argument(
        "--ligand-id",
        help="Ligand ID for YAML"
    )
    
    args = parser.parse_args()
    
    # Load config if provided
    if args.config:
        config = load_config(args.config)
        input_csv = args.input or config['input']['csv_path'].replace('.csv', '_with_sdf_with_diffdock.csv')
        fasta_path = args.fasta or config['target']['fasta_path']
        output_dir = args.output_dir or os.path.join(config['output']['base_dir'], 'yamls')
        plip_output = args.plip_output or os.path.join(config['output']['base_dir'], 'plip')
        msa_path = args.msa or config['target'].get('msa_path')
        micromamba_exe = args.micromamba or config.get('micromamba', {}).get('executable', 'micromamba')
        micromamba_env = args.env or config.get('micromamba', {}).get('plip_env', 'plip')
        max_distance = args.max_distance or config['constraints']['max_distance']
        constraint_type = args.constraint_type or ('both' if 'contact' in config['constraints']['constraint_types'] else 'pocket')
        ligand_id = args.ligand_id or config['constraints']['ligand_id']
    else:
        if not args.input or not args.fasta or not args.output_dir or not args.plip_output:
            parser.error("Either --config or --input, --fasta, --output-dir, and --plip-output are required")
        input_csv = args.input
        fasta_path = args.fasta
        output_dir = args.output_dir
        plip_output = args.plip_output
        msa_path = args.msa
        micromamba_exe = args.micromamba or "micromamba"
        micromamba_env = args.env or "plip"
        max_distance = args.max_distance or 6.0
        constraint_type = args.constraint_type or 'pocket'
        ligand_id = args.ligand_id or "LIG"

    # PLIP is invoked without a distance flag; micromamba is resolved the same
    # way in every script — see constraints/utils.py.
    micromamba_exe = require_micromamba(micromamba_exe)
    
    if not os.path.exists(input_csv):
        print(f"[ERROR] Input file not found: {input_csv}")
        sys.exit(1)
    
    if not os.path.exists(fasta_path):
        print(f"[ERROR] FASTA file not found: {fasta_path}")
        sys.exit(1)
    
    df = pd.read_csv(input_csv)
    
    if 'complex_path' not in df.columns:
        print(f"[ERROR] 'complex_path' column not found. Run script 06 first.")
        sys.exit(1)
    
    sequence = read_fasta_sequence(fasta_path)
    
    ensure_dir(output_dir)
    ensure_dir(plip_output)
    
    print(f"\n{'='*70}")
    print(f"  PLIP ANALYSIS & CONSTRAINT YAML GENERATION")
    print(f"{'='*70}\n")
    print(f"Input:            {input_csv}")
    print(f"FASTA:            {fasta_path}")
    print(f"YAMLs:            {output_dir}")
    print(f"PLIP reports:     {plip_output}")
    print(f"Ligands:          {len(df)}")
    print(f"Constraint type:  {constraint_type}")
    print()
    print(f"PLIP detection:   default binding-site radius (BS_DIST, 7.5 A in PLIP 2.3.x)")
    print(f"                  no distance argument is passed to plip")
    print(f"Boltz max_distance: {max_distance} A")
    print(f"                  written into the pocket constraint; how far the")
    print(f"                  binder may sit from the contacts PLIP found\n")
    
    yaml_default_paths = []
    yaml_pocket_paths = []
    yaml_contact_paths = []
    num_contacts_list = []
    
    successful = 0
    no_constraints = 0
    failed = 0
    
    for idx, row in df.iterrows():
        molport_number = row.get('molport_number', idx)
        smiles = row.get('smiles', '')
        complex_path = row.get('complex_path')
        
        # Generate short ligand ID for Boltz-2 (max 5 chars)
        if 'rank' in df.columns:
            rank = row.get('rank', idx + 1)
            short_ligand_id = f"R{int(rank)}"
        else:
            short_ligand_id = ligand_id  # Use config default (usually "LIG")
        
        if pd.isna(complex_path) or not os.path.exists(str(complex_path)):
            print(f"  [SKIP] {molport_number}: No complex file")
            yaml_default_paths.append(None)
            yaml_pocket_paths.append(None)
            yaml_contact_paths.append(None)
            num_contacts_list.append(0)
            failed += 1
            continue
        
        print(f"  [{idx+1}/{len(df)}] {molport_number}...", end=" ")
        
        plip_ligand_output = os.path.join(plip_output, str(molport_number))
        xml_path = run_plip(complex_path, plip_ligand_output, micromamba_exe, micromamba_env)
        
        if not xml_path or not os.path.exists(xml_path):
            print(f"✗ (PLIP failed)")
            yaml_default_paths.append(None)
            yaml_pocket_paths.append(None)
            yaml_contact_paths.append(None)
            num_contacts_list.append(0)
            failed += 1
            continue
        
        contacts = extract_plip_contacts(xml_path)
        num_contacts = len(contacts)
        num_contacts_list.append(num_contacts)
        
        # Use molport_number for filename, short_ligand_id for YAML content
        yaml_default = os.path.join(output_dir, f"{molport_number}_default.yaml")
        with open(yaml_default, 'w') as f:
            f.write(generate_yaml(sequence, short_ligand_id, smiles, msa_path, 
                                 None, max_distance, 'pocket'))
        yaml_default_paths.append(yaml_default)
        
        if num_contacts > 0:
            yaml_pocket = os.path.join(output_dir, f"{molport_number}_pocket.yaml")
            with open(yaml_pocket, 'w') as f:
                f.write(generate_yaml(sequence, short_ligand_id, smiles, msa_path,
                                     contacts, max_distance, 'pocket'))
            yaml_pocket_paths.append(yaml_pocket)
            
            if constraint_type in ['contact', 'both']:
                yaml_contact = os.path.join(output_dir, f"{molport_number}_contact.yaml")
                with open(yaml_contact, 'w') as f:
                    f.write(generate_yaml(sequence, short_ligand_id, smiles, msa_path,
                                         contacts, max_distance, 'contact'))
                yaml_contact_paths.append(yaml_contact)
            else:
                yaml_contact_paths.append(None)
            
            print(f"✓ ({num_contacts} contacts)")
            successful += 1
        else:
            yaml_pocket_paths.append(None)
            yaml_contact_paths.append(None)
            print(f"⚠ (no contacts)")
            no_constraints += 1
    
    df['yaml_default'] = yaml_default_paths
    df['yaml_pocket'] = yaml_pocket_paths
    df['yaml_contact'] = yaml_contact_paths
    df['num_contacts'] = num_contacts_list
    
    output_csv = input_csv.replace('.csv', '_with_yamls.csv')
    df.to_csv(output_csv, index=False)
    
    print(f"\n{'='*70}")
    print(f"  SUMMARY")
    print(f"{'='*70}")
    print(f"With constraints:    {successful}/{len(df)}")
    print(f"Without constraints: {no_constraints}/{len(df)}")
    print(f"Failed:              {failed}/{len(df)}")
    
    if num_contacts_list:
        valid_contacts = [n for n in num_contacts_list if n > 0]
        if valid_contacts:
            print(f"\nContact statistics:")
            print(f"  Min:  {min(valid_contacts)}")
            print(f"  Max:  {max(valid_contacts)}")
            print(f"  Mean: {sum(valid_contacts)/len(valid_contacts):.1f}")
    
    print(f"\nUpdated CSV saved: {output_csv}\n")


if __name__ == "__main__":
    main()
