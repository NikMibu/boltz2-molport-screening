#!/usr/bin/env python3
"""Shared resolution of the external tools the constraint pipeline shells out to.

Both micromamba and DiffDock are invoked through ``subprocess``, which does not
go through a shell. That makes them harder to find than they look from an
interactive prompt, and the failure modes are quiet. The lookups live here so
that every script — and the setup validator — asks the same question the same
way. A validator that checks something the code never calls is worse than no
validator: it reports green where things are broken and red where they work.
"""
import os
import shutil
from typing import Optional

_MICROMAMBA_CANDIDATES = (
    "~/.local/bin/micromamba",
    "~/micromamba/bin/micromamba",
    "~/bin/micromamba",
    "/usr/local/bin/micromamba",
    "/opt/micromamba/bin/micromamba",
)


def resolve_micromamba(configured: str = "micromamba") -> Optional[str]:
    """Locate the micromamba executable, or return None.

    ``shutil.which`` is not enough. The standard micromamba install defines a
    *shell function* named ``micromamba`` that wraps the binary, and puts the
    binary somewhere that is not on PATH. So ``micromamba run ...`` works when
    typed interactively while ``which micromamba`` finds nothing — and
    subprocess, which does not go through the shell, cannot call it at all.

    The shell hook exports the real path as ``$MAMBA_EXE``, which is the
    reliable way to find it from Python.

    Kept deliberately identical to ``plip_pipeline/utils.py`` in the
    plip_constraints_pipeline repository.
    """
    # 1. An explicit path in the config always wins.
    if configured and os.path.sep in configured:
        candidate = os.path.expanduser(os.path.expandvars(configured))
        if os.path.exists(candidate):
            return candidate

    # 2. Exported by the micromamba shell hook.
    mamba_exe = os.environ.get("MAMBA_EXE")
    if mamba_exe and os.path.exists(mamba_exe):
        return mamba_exe

    # 3. Plain PATH lookup.
    found = shutil.which(configured or "micromamba")
    if found:
        return found

    # 4. Usual install locations.
    for candidate in _MICROMAMBA_CANDIDATES:
        candidate = os.path.expanduser(candidate)
        if os.path.exists(candidate):
            return candidate

    return None


def require_micromamba(configured: str = "micromamba") -> str:
    """Resolve micromamba or exit with an explanation."""
    resolved = resolve_micromamba(configured)
    if resolved:
        return resolved
    raise SystemExit(
        "ERROR: micromamba not found.\n"
        f"       Tried: the config value '{configured}', $MAMBA_EXE, the PATH, "
        "and the usual install locations.\n"
        "       The standard install leaves micromamba as a shell function, so it\n"
        "       works when you type it but is invisible to subprocess. Export the\n"
        "       binary explicitly:\n"
        "         export MAMBA_EXE=/path/to/bin/micromamba\n"
        "       or set micromamba.executable to a full path in the config."
    )


def resolve_diffdock_home(configured: Optional[str] = None) -> Optional[str]:
    """Locate the DiffDock repository.

    ``$DIFFDOCK_HOME`` takes precedence over the config, so that one checkout
    can be pointed at a different DiffDock without editing the config. Returns
    None if neither resolves to a directory holding ``inference.py``.
    """
    for candidate in (os.environ.get("DIFFDOCK_HOME"), configured):
        if not candidate:
            continue
        path = os.path.expanduser(os.path.expandvars(candidate))
        if "$" in path:
            # An unexpanded ${DIFFDOCK_HOME} from the example config.
            continue
        if os.path.isdir(path):
            return path
    return None


def require_diffdock_home(configured: Optional[str] = None) -> str:
    """Resolve the DiffDock repository or exit with an explanation."""
    resolved = resolve_diffdock_home(configured)
    if resolved:
        if not os.path.isfile(os.path.join(resolved, "inference.py")):
            raise SystemExit(
                f"ERROR: {resolved} exists but holds no inference.py.\n"
                "       That is not a DiffDock checkout."
            )
        return resolved
    raise SystemExit(
        "ERROR: DiffDock repository not found.\n"
        f"       $DIFFDOCK_HOME is {os.environ.get('DIFFDOCK_HOME') or 'not set'}, "
        f"config value is {configured or 'unset'}.\n"
        "       Point at your checkout:\n"
        "         export DIFFDOCK_HOME=/path/to/DiffDock"
    )
