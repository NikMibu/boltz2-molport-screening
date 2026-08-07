#!/usr/bin/env python3
"""Point Boltz-2 input YAMLs at a pre-computed MSA.

Handles three cases:

1. ``msa: empty``      -> replaced with the given path
2. no ``msa:`` field   -> a line is inserted after the protein sequence
3. an ``msa:`` path is already present -> skipped, unless ``--replace``

Case 3 matters more than it looks. YAMLs generated on another machine carry an
absolute MSA path that does not exist here; skipping them leaves the whole set
pointing at nothing and Boltz-2 only finds out at run time. Those files are
therefore reported explicitly, with the path they currently hold.
"""
import argparse
import re
import sys
from pathlib import Path

MSA_LINE = re.compile(r"^\s*msa:\s*(?P<value>.*)$", re.MULTILINE)


def update_yaml_with_msa(content: str, msa_path: str) -> str:
    """Insert or replace the MSA path in a Boltz YAML."""

    # Case 1: replace the "empty" placeholder
    if "msa: empty" in content:
        return content.replace("msa: empty", f"msa: {msa_path}")

    # Case 2: no MSA field — insert after the protein sequence line
    lines = content.split("\n")
    new_lines = []
    in_protein_block = False
    msa_added = False

    for line in lines:
        new_lines.append(line)

        if "- protein:" in line or ("protein:" in line and "- " in line):
            in_protein_block = True

        if in_protein_block and ("- ligand:" in line or "- rna:" in line or "- dna:" in line):
            in_protein_block = False

        if in_protein_block and "sequence:" in line and not msa_added:
            indent = len(line) - len(line.lstrip())
            new_lines.append(" " * indent + f"msa: {msa_path}")
            msa_added = True
            in_protein_block = False  # add the MSA only once

    return "\n".join(new_lines)


def replace_existing_msa(content: str, msa_path: str) -> str:
    """Overwrite an MSA path that is already present, keeping the indentation."""

    def sub(match: re.Match) -> str:
        line = match.group(0)
        indent = len(line) - len(line.lstrip())
        return " " * indent + f"msa: {msa_path}"

    return MSA_LINE.sub(sub, content, count=1)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Insert or replace the MSA path in a folder of Boltz-2 YAMLs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  # add the MSA to freshly generated YAMLs
  python 02_update_msa.py output/ca2/yamls data/msa/P00918.a3m

  # repoint YAMLs that already carry a path from another machine
  python 02_update_msa.py output/ca2/yamls data/msa/P00918.a3m --replace
""",
    )
    ap.add_argument("yaml_dir", help="Directory containing the .yaml files")
    ap.add_argument("msa_path", help="Path to the .a3m alignment")
    ap.add_argument(
        "--replace",
        action="store_true",
        help="Overwrite an MSA path that is already present (default: leave it alone)",
    )
    ap.add_argument(
        "--relative",
        action="store_true",
        help="Write the path as given instead of resolving it to an absolute path",
    )
    args = ap.parse_args()

    yaml_dir = Path(args.yaml_dir).expanduser().resolve()
    msa_arg = Path(args.msa_path).expanduser()
    msa_path = msa_arg if args.relative else msa_arg.resolve()

    if not yaml_dir.is_dir():
        print(f"ERROR: directory not found: {yaml_dir}", file=sys.stderr)
        return 2

    if not msa_arg.expanduser().exists():
        print(f"ERROR: MSA file not found: {msa_arg}", file=sys.stderr)
        return 2

    print(f"Updating YAMLs in: {yaml_dir}")
    print(f"With MSA path: {msa_path}")

    count = 0
    updated = 0
    replaced = 0
    unchanged = 0
    stale = []

    for yaml_file in sorted(yaml_dir.glob("*.yaml")):
        count += 1
        content = yaml_file.read_text()
        existing = MSA_LINE.search(content)

        if existing and existing.group("value").strip() not in ("", "empty"):
            current = existing.group("value").strip()
            if current == str(msa_path):
                unchanged += 1
                continue
            if not args.replace:
                stale.append((yaml_file.name, current))
                continue
            new_content = replace_existing_msa(content, str(msa_path))
            if new_content != content:
                yaml_file.write_text(new_content)
                replaced += 1
            continue

        new_content = update_yaml_with_msa(content, str(msa_path))
        if new_content != content:
            yaml_file.write_text(new_content)
            updated += 1

    if count == 0:
        print(f"ERROR: no .yaml files found in {yaml_dir}", file=sys.stderr)
        return 1

    print(f"Processed {count} file(s)")
    print(f"  MSA added:      {updated}")
    if replaced:
        print(f"  MSA replaced:   {replaced}")
    if unchanged:
        print(f"  Already correct: {unchanged}")

    if stale:
        shown = stale[:3]
        print(
            f"\nWARNING: {len(stale)} file(s) already carry a different MSA path and were "
            f"left untouched:",
            file=sys.stderr,
        )
        for name, current in shown:
            print(f"  {name}: {current}", file=sys.stderr)
        if len(stale) > len(shown):
            print(f"  ... and {len(stale) - len(shown)} more", file=sys.stderr)
        print(
            "Those paths are not checked for existence. If they come from another "
            "machine, Boltz-2 will fail at run time — re-run with --replace.",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
