#!/usr/bin/env bash
# Boltz-2 inference over a folder of YAMLs.
#   run_boltz.sh CONFIG TARGET
# Normally invoked by run_screening.sh, not directly.
set -euo pipefail

usage() {
    cat <<EOF
Usage: bash screening/run_boltz.sh CONFIG TARGET [--smoke]

Runs 'boltz predict' over the YAML folder for one target, using the inference
parameters from CONFIG. Normally invoked by run_screening.sh.

  CONFIG    Pipeline config (e.g. screening/screening_config.yaml)
  TARGET    Target key from the config (ca2, ca4, ca7)
  --smoke   Use the reduced smoke_* inference parameters

The Boltz executable is taken from \$BOLTZ_EXE if set, otherwise 'boltz' on the
PATH. Boltz-2 often lives in its own venv rather than on the PATH:

  export BOLTZ_EXE=/path/to/exp_env/bin/boltz
EOF
}

if [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
    usage; exit 0
fi

CONFIG="${1:?CONFIG required — see --help}"
TARGET="${2:?TARGET required — see --help}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# Prefix for the inference-parameter keys: "" (full run) or "smoke_"
P=""
[ "${3:-}" = "--smoke" ] && P="smoke_"

read_cfg() {
    python3 - "$CONFIG" "$1" "${2:-}" <<'PY'
import sys, yaml
with open(sys.argv[1]) as f:
    cfg = yaml.safe_load(f)
v = cfg
for k in sys.argv[2].split("."):
    v = v[k]
if len(sys.argv) > 3 and sys.argv[3]:
    v = str(v).replace("{target}", sys.argv[3])
print(v)
PY
}

YAML_DIR="$ROOT/$(read_cfg yaml_dir "$TARGET")"
RESULTS_DIR="$ROOT/$(read_cfg results_dir "$TARGET")"
RECYCLING=$(read_cfg "${P}recycling_steps")
SAMPLING=$(read_cfg "${P}sampling_steps")
DIFF_SAMPLES_AFF=$(read_cfg "${P}diffusion_samples_affinity")
SAMPLING_AFF=$(read_cfg "${P}sampling_steps_affinity")
NO_KERNELS=$(read_cfg no_kernels)

# Boltz has no --version; the executable is resolved here so that a missing
# install fails before the YAMLs are handed over, not halfway through a run.
BOLTZ="${BOLTZ_EXE:-}"
if [ -z "$BOLTZ" ]; then
    if command -v boltz >/dev/null 2>&1; then
        BOLTZ="$(command -v boltz)"
    else
        echo "ERROR: 'boltz' not found on the PATH and \$BOLTZ_EXE is not set." >&2
        echo "       Boltz-2 is commonly installed in its own venv:" >&2
        echo "       export BOLTZ_EXE=/path/to/exp_env/bin/boltz" >&2
        exit 1
    fi
fi
if [ ! -x "$BOLTZ" ]; then
    echo "ERROR: \$BOLTZ_EXE is not an executable file: $BOLTZ" >&2
    exit 1
fi

if [ ! -d "$YAML_DIR" ]; then
    echo "ERROR: YAML folder not found: $YAML_DIR" >&2
    echo "       Run 01_generate_yamls.py first." >&2
    exit 1
fi
YAML_COUNT=$(find "$YAML_DIR" -name '*.yaml' 2>/dev/null | wc -l || true)
if [ "$YAML_COUNT" -eq 0 ]; then
    echo "ERROR: no .yaml files in $YAML_DIR" >&2
    exit 1
fi

EXTRA=()
if [ "$NO_KERNELS" = "True" ] || [ "$NO_KERNELS" = "true" ]; then
    EXTRA+=(--no_kernels)
elif grep -qi microsoft /proc/version 2>/dev/null; then
    # WSL's NVML has no nvmlDeviceGetNumGpuCores; the cuEquivariance kernels
    # call it during import and the run dies with NVMLError_NotSupported after
    # the model is already loaded. Cheaper to say so now than 45 seconds in.
    echo "WARNING: this looks like WSL and no_kernels is false." >&2
    echo "         Boltz's cuEquivariance kernels crash there with" >&2
    echo "         'NVMLError_NotSupported' once prediction starts." >&2
    echo "         Set 'no_kernels: true' in $(basename "$CONFIG")." >&2
    echo "" >&2
fi

mkdir -p "$RESULTS_DIR"

echo "Boltz:      $BOLTZ"
echo "Input:      $YAML_DIR ($YAML_COUNT YAMLs)"
echo "Output:     $RESULTS_DIR"
echo "Parameters: ${RECYCLING}/${SAMPLING}/${DIFF_SAMPLES_AFF}/${SAMPLING_AFF}${P:+ (smoke)}"
if [ ${#EXTRA[@]} -gt 0 ]; then
    echo "cuEquivariance kernels disabled (no_kernels: true)"
fi

"$BOLTZ" predict "$YAML_DIR" \
    --out_dir "$RESULTS_DIR" \
    --recycling_steps "$RECYCLING" \
    --sampling_steps "$SAMPLING" \
    --diffusion_samples_affinity "$DIFF_SAMPLES_AFF" \
    --sampling_steps_affinity "$SAMPLING_AFF" \
    "${EXTRA[@]}"
