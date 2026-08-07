#!/usr/bin/env python3
"""Setup validation.

Checks what the pipelines actually call, before a run starts, and reads the
configs so that environment names and paths are the ones the code will use
rather than assumptions baked in here. A validator that checks its own idea of
the setup is worse than none: it reports green where things are broken and red
where they work.

Groups:
  core        - needed for anything
  screening   - run_screening.sh: sequences, alignments, Boltz-2
  constraints - run_constraints.sh: micromamba, DiffDock, PLIP

Constraints failures are warnings: the screening pipeline stands on its own and
does not need DiffDock or PLIP.

Exit codes: 0 all required checks passed, 1 a required check failed.
"""
import argparse
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "constraints"))
try:
    # Share the pipeline's own resolvers so the two cannot drift apart.
    from utils import resolve_micromamba, resolve_diffdock_home
except Exception:  # missing deps are reported by the core checks below
    resolve_micromamba = None
    resolve_diffdock_home = None

ROOT = os.path.dirname(os.path.abspath(__file__))

OK = "[ ok ]"
FAIL = "[fail]"
WARN = "[warn]"


class Report:
    def __init__(self):
        self.required = []
        self.optional = []

    def add(self, name, passed, detail="", optional=False, hint=""):
        """detail is shown always, hint only when the check did not pass."""
        mark = OK if passed else (WARN if optional else FAIL)
        suffix = detail if passed else (hint or detail)
        print(f"  {mark} {name}" + (f" - {suffix}" if suffix else ""))
        (self.optional if optional else self.required).append((name, passed))
        return passed


def load_config(path):
    try:
        import yaml
    except ImportError:
        return None, "PyYAML is not installed, cannot read the config"
    if not os.path.exists(path):
        return None, f"{path} not found"
    try:
        with open(path) as f:
            return yaml.safe_load(f), ""
    except Exception as e:
        return None, str(e)


def rel(path):
    """Resolve a config path against the repository root."""
    return path if os.path.isabs(path) else os.path.join(ROOT, path)


def fasta_sequence(path):
    seq = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line)
    return "".join(seq).upper()


def a3m_query_sequence(path):
    """First sequence entry of an a3m — by construction the target sequence.

    Insertions relative to the query are written lowercase, so the query line
    itself carries none; it is read verbatim.
    """
    with open(path, encoding="utf-8") as fh:
        seen_header = False
        seq = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                if seen_header:
                    break
                seen_header = True
                continue
            if seen_header:
                seq.append(line.strip())
    return "".join(seq).upper()


def check_core(rep):
    print("\nCore")
    v = sys.version_info
    rep.add("Python >= 3.10", v >= (3, 10), f"{v.major}.{v.minor}.{v.micro}")

    for label, module in [("pandas", "pandas"), ("numpy", "numpy"),
                          ("PyYAML", "yaml"), ("RDKit", "rdkit"),
                          ("matplotlib", "matplotlib"), ("seaborn", "seaborn"),
                          ("requests", "requests")]:
        try:
            __import__(module)
            rep.add(label, True)
        except ImportError:
            rep.add(label, False, hint="pip install -r requirements.txt")

    try:
        __import__("Bio")
        rep.add("BioPython", True)
    except ImportError:
        rep.add("BioPython", False, optional=True,
                hint="only needed for the constraint pipeline (06 strips HETATM records)")


def resolve_boltz(explicit=None):
    """Locate the boltz executable the way the wrappers do."""
    import shutil
    candidate = explicit or os.environ.get("BOLTZ_EXE")
    if candidate:
        return candidate if os.path.isfile(candidate) and os.access(candidate, os.X_OK) else None
    return shutil.which("boltz")


def boltz_version(exe):
    """Read the version from the package metadata.

    Boltz has no --version, and --help imports torch, which costs 45 seconds.
    The interpreter next to the executable knows the answer in under a second.
    """
    python = os.path.join(os.path.dirname(exe), "python")
    if not os.path.isfile(python):
        python = sys.executable
    try:
        result = subprocess.run(
            [python, "-c",
             "import importlib.metadata as m; print(m.version('boltz'))"],
            capture_output=True, timeout=60, text=True)
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    except Exception:
        pass
    return None


def check_screening(rep, cfg):
    print("\nScreening workflow")
    if not cfg:
        rep.add("screening config", False, hint="could not be read, skipping the rest")
        return

    targets = cfg.get("targets", {})
    rep.add("targets defined", bool(targets), ", ".join(targets) or "none")

    for key, spec in targets.items():
        fasta = rel(spec.get("fasta", ""))
        msa = rel(spec.get("msa", ""))

        if not os.path.isfile(fasta):
            rep.add(f"{key}: FASTA", False, hint=f"not found: {fasta}")
            continue
        if not os.path.isfile(msa):
            rep.add(f"{key}: MSA", False, hint=f"not found: {msa}")
            continue

        # The pairing check. An a3m built from a different sequence is accepted
        # by Boltz-2 without complaint and quietly degrades the prediction —
        # the trap that the trimmed and untrimmed hCA IV entries set, since
        # they carry the same UniProt header.
        try:
            seq = fasta_sequence(fasta)
            query = a3m_query_sequence(msa)
        except Exception as e:
            rep.add(f"{key}: FASTA/MSA pair", False, hint=f"could not read: {e}")
            continue

        if seq == query:
            rep.add(f"{key}: FASTA/MSA pair", True, f"{len(seq)} residues")
        elif len(seq) == len(query):
            rep.add(f"{key}: FASTA/MSA pair", False,
                    hint=f"same length ({len(seq)}) but different sequence — "
                         f"{os.path.basename(msa)} was built from another protein")
        else:
            rep.add(f"{key}: FASTA/MSA pair", False,
                    hint=f"{os.path.basename(fasta)} has {len(seq)} residues, "
                         f"{os.path.basename(msa)} was built from {len(query)} — "
                         f"trimmed and untrimmed sequences are easy to mix up")

    exe = resolve_boltz()
    if not exe:
        rep.add("Boltz-2", False,
                hint="not on the PATH and $BOLTZ_EXE is not set. Boltz commonly "
                     "lives in its own venv: export BOLTZ_EXE=/path/to/bin/boltz")
    else:
        version = boltz_version(exe)
        rep.add("Boltz-2", True, f"{version or 'version unknown'} at {exe}")
        if version and version != "2.2.1":
            rep.add("Boltz-2 version 2.2.1", False, optional=True,
                    hint=f"found {version}; the published results were produced with 2.2.1")


def check_constraints(rep, cfg):
    print("\nConstraint workflow (optional)")
    if resolve_micromamba is None or resolve_diffdock_home is None:
        rep.add("constraints/utils.py", False, optional=True,
                hint="could not be imported, skipping the rest")
        return
    if not cfg:
        rep.add("constraints config", False, optional=True,
                hint="constraints/config.yaml not found. "
                     "cp constraints/config.example.yaml constraints/config.yaml")
        return

    mm = cfg.get("micromamba", {})
    exe = resolve_micromamba(mm.get("executable", "micromamba"))
    if not exe:
        rep.add("micromamba", False, optional=True,
                hint="not found via the config, $MAMBA_EXE, the PATH or the usual "
                     "locations. The default install leaves it as a shell function; "
                     "export MAMBA_EXE=/path/to/bin/micromamba")
        return
    rep.add("micromamba", True, exe)

    for env_key, default, tool in [("diffdock_env", "diffdock", None),
                                   ("plip_env", "plip", "plip")]:
        env = mm.get(env_key, default)
        try:
            result = subprocess.run([exe, "env", "list"],
                                    capture_output=True, timeout=30, text=True)
            found = any(line.split() and line.split()[0] == env
                        for line in result.stdout.splitlines())
        except Exception as e:
            rep.add(f"environment '{env}'", False, optional=True, hint=str(e))
            continue
        rep.add(f"environment '{env}'", found, optional=True,
                hint=f"micromamba knows no environment '{env}'")

        if found and tool:
            try:
                result = subprocess.run([exe, "run", "-n", env, tool, "--help"],
                                        capture_output=True, timeout=120)
                runs = result.returncode in (0, 1, 2)
            except Exception:
                runs = False
            rep.add(f"{tool} runs in '{env}'", runs, optional=True,
                    hint=f"'{tool}' is not callable inside the environment")

    repo = resolve_diffdock_home(cfg.get("diffdock", {}).get("repo_path"))
    if not repo:
        rep.add("DiffDock repository", False, optional=True,
                hint="$DIFFDOCK_HOME is not set and the config path does not exist")
    elif not os.path.isfile(os.path.join(repo, "inference.py")):
        rep.add("DiffDock repository", False, optional=True,
                hint=f"{repo} holds no inference.py")
    else:
        rep.add("DiffDock repository", True, repo)

    receptor = rel(cfg.get("target", {}).get("protein_structure", ""))
    rep.add("apo receptor PDB", os.path.isfile(receptor), receptor,
            optional=True, hint=f"not found: {receptor}")


def check_network(rep):
    print("\nNetwork (optional)")
    try:
        import requests
    except ImportError:
        rep.add("requests", False, optional=True, hint="pip install requests")
        return
    for label, url in [("UniProt", "https://rest.uniprot.org/uniprotkb/P00918.fasta"),
                       ("RCSB", "https://files.rcsb.org/download/1CNC.pdb")]:
        try:
            r = requests.get(url, timeout=20)
            rep.add(label, r.status_code == 200, f"HTTP {r.status_code}",
                    optional=True)
        except Exception as e:
            rep.add(label, False, optional=True, hint=str(e))


def main():
    p = argparse.ArgumentParser(
        description="Check that the pipelines' external dependencies are in place",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  python3 validate_setup.py
  python3 validate_setup.py --skip-network
  python3 validate_setup.py --constraints-config constraints/my.yaml

Constraint-workflow failures are warnings: the screening pipeline needs neither
DiffDock nor PLIP.
""",
    )
    p.add_argument("--screening-config",
                   default=os.path.join(ROOT, "screening", "screening_config.yaml"),
                   help="screening config to validate against")
    p.add_argument("--constraints-config",
                   default=os.path.join(ROOT, "constraints", "config.yaml"),
                   help="constraint config to validate against")
    p.add_argument("--skip-network", action="store_true",
                   help="do not contact UniProt or RCSB")
    args = p.parse_args()

    print("=" * 62)
    print("  boltz2-molport-screening — setup validation")
    print("=" * 62)

    rep = Report()
    check_core(rep)

    scr_cfg, scr_err = load_config(args.screening_config)
    if scr_err:
        print(f"\n  note: {scr_err}")
    check_screening(rep, scr_cfg)

    con_cfg, _ = load_config(args.constraints_config)
    check_constraints(rep, con_cfg)

    if not args.skip_network:
        check_network(rep)

    failed = [n for n, ok in rep.required if not ok]
    warned = [n for n, ok in rep.optional if not ok]

    print("\n" + "=" * 62)
    print(f"  required: {len(rep.required) - len(failed)}/{len(rep.required)} passed")
    print(f"  optional: {len(rep.optional) - len(warned)}/{len(rep.optional)} passed")
    if warned:
        print(f"  warnings: {', '.join(warned)}")
    if failed:
        print(f"  FAILED:   {', '.join(failed)}")
        print("=" * 62)
        return 1
    print("  ready")
    print("=" * 62)
    return 0


if __name__ == "__main__":
    sys.exit(main())
