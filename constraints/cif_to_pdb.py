#!/usr/bin/env python3
"""Convert a Boltz-2 mmCIF prediction to PDB.

Boltz-2 writes predicted complexes as mmCIF. PDB is what the tools downstream of
it expect — a viewer such as ChimeraX or PyMOL to look at a pose, or PLIP to
profile its interactions.

The catch is the chain identifier. PDB allots it a single character, and Boltz
names the ligand chain after the ligand id in the input YAML, which
01_generate_yamls.py writes as L1, L2, ... — or L4635 in a large run. Those
names do not fit, so a naive conversion drops the ligand and leaves a bare
protein that looks fine until you notice nothing is bound.

Every chain is therefore renamed to one character first, keeping the leading
character where it is free. The ligand keeps chain L, which is where PLIP
reports it elsewhere in this pipeline.

Not part of run_constraints.sh: that pipeline runs PLIP on DiffDock poses, which
are already PDB. This is for looking at, or profiling, what Boltz-2 itself
predicted.
"""
import argparse
import os
import sys

FALLBACK_CHARS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"


def choose_single_char_chain_ids(structure) -> dict:
    """Rename chains in place to one character each. Returns the mapping."""
    model = structure[0]
    used = set()
    fallback = iter(FALLBACK_CHARS)
    mapping = {}

    for chain in model:
        old = chain.name
        if len(old) == 1 and old not in used:
            new = old
        elif old and old[0].isalnum() and old[0].upper() not in used:
            # Keep the first character. Boltz derives the chain name from the
            # ligand id in the YAML, which 01_generate_yamls.py writes as
            # L1, L2, ... — so the ligand keeps chain L, which is what PLIP
            # reports it under everywhere else in this pipeline.
            new = old[0].upper()
        elif old.startswith("R") and "L" not in used:
            # Older runs named the ligand chain R<something>.
            new = "L"
        else:
            new = None
            for candidate in fallback:
                if candidate not in used:
                    new = candidate
                    break
            if new is None:
                raise SystemExit(
                    f"ERROR: {len(model)} chains is more than the {len(FALLBACK_CHARS)} "
                    "single-character names PDB allows."
                )
        chain.name = new
        used.add(new)
        mapping[old] = new
    return mapping


def convert(cif_path, pdb_path, verbose=True):
    import gemmi

    structure = gemmi.read_structure(str(cif_path))
    if len(structure) == 0:
        raise SystemExit(f"ERROR: no model in {cif_path}")
    mapping = choose_single_char_chain_ids(structure)
    structure.setup_entities()
    os.makedirs(os.path.dirname(os.path.abspath(pdb_path)) or ".", exist_ok=True)
    structure.write_pdb(str(pdb_path))
    if verbose:
        renamed = {o: n for o, n in mapping.items() if o != n}
        note = f"  chains: {', '.join(f'{o}->{n}' for o, n in renamed.items())}" if renamed else ""
        print(f"{os.path.basename(cif_path)} -> {pdb_path}{note}")
    return mapping


def main():
    ap = argparse.ArgumentParser(
        description="Convert Boltz-2 mmCIF complexes to PLIP-readable PDB",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  # one pose, to open in ChimeraX or PyMOL
  python constraints/cif_to_pdb.py --input pred.cif --output pred.pdb

  # every prediction of a Boltz run
  python constraints/cif_to_pdb.py --input results/.../boltz_output --out-dir pdb/

PDB allows one character per chain. Boltz chain names are longer, so they are
renamed; the ligand chain becomes L. Without that the ligand is dropped and the
file opens as a bare protein.
""",
    )
    ap.add_argument("--input", required=True,
                    help="An .cif file, or a folder to search recursively")
    ap.add_argument("--output", help="Output file (only with a single --input file)")
    ap.add_argument("--out-dir", help="Output folder (required when --input is a folder)")
    args = ap.parse_args()

    try:
        import gemmi  # noqa: F401
    except ImportError:
        raise SystemExit("ERROR: gemmi is not installed. pip install gemmi")

    if os.path.isfile(args.input):
        out = args.output or os.path.splitext(args.input)[0] + ".pdb"
        convert(args.input, out)
        return 0

    if not os.path.isdir(args.input):
        raise SystemExit(f"ERROR: not found: {args.input}")
    if not args.out_dir:
        raise SystemExit("ERROR: --out-dir is required when --input is a folder")

    cifs = []
    for root, _dirs, files in os.walk(args.input):
        for name in sorted(files):
            if name.endswith(".cif"):
                cifs.append(os.path.join(root, name))

    if not cifs:
        raise SystemExit(f"ERROR: no .cif files under {args.input}")

    failed = 0
    for cif in cifs:
        stem = os.path.splitext(os.path.basename(cif))[0]
        try:
            convert(cif, os.path.join(args.out_dir, stem + ".pdb"))
        except SystemExit:
            raise
        except Exception as exc:
            print(f"  FAILED {cif}: {exc}", file=sys.stderr)
            failed += 1

    converted = len(cifs) - failed
    print(f"\nConverted {converted}/{len(cifs)} structure(s) into {args.out_dir}")
    if converted == 0:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
