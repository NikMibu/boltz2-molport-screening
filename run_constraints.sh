#!/usr/bin/env bash
# Constraint-based re-ranking with DiffDock and PLIP (Chapter 6.6).
#   bash run_constraints.sh --config constraints/config.yaml
#   bash run_constraints.sh --config constraints/config.yaml --smoke
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
CONFIG="${ROOT}/constraints/config.yaml"
CONSTRAINTS="${ROOT}/constraints"

usage() {
    cat <<EOF
Usage: bash run_constraints.sh [OPTIONS]

Re-ranks a hit list under interaction constraints: 3D ligands -> DiffDock poses
-> PLIP contacts -> two Boltz-2 YAMLs per compound (unconstrained and pocket)
-> comparison of the two rankings.

Options:
  --config PATH     Pipeline config (default: constraints/config.yaml)
  --smoke           Reduced run: few compounds, few DiffDock poses, reduced
                    inference parameters. Minutes instead of a day. Proves the
                    chain runs; reproduces nothing.
  --steps LIST      Comma-separated subset of: prep,dock,plip,boltz,analyze
                    (default: all five)
  --skip-boltz      Alias for --steps prep,dock,plip
  -h, --help        Show this message and exit

Examples:
  bash run_constraints.sh --config constraints/config.yaml --smoke
  bash run_constraints.sh --steps analyze

Requirements beyond the screening pipeline:
  \$DIFFDOCK_HOME   path to a DiffDock checkout
  \$MAMBA_EXE       path to the micromamba binary, if micromamba is installed
                   as a shell function (the default install is)
  \$BOLTZ_EXE       path to the boltz executable, if not on the PATH
Run 'python3 validate_setup.py' to check all of them before a long run.

Inference runs at 3/200/5/200, the Boltz-2 defaults. This differs from
run_screening.sh, which uses the reduced 2/50/1/50 — see the configs.
EOF
}

SMOKE=false
STEPS="prep,dock,plip,boltz,analyze"

while [ $# -gt 0 ]; do
    case "$1" in
        --config)      CONFIG="${2:?--config requires a path}"; shift 2 ;;
        --config=*)    CONFIG="${1#*=}"; shift ;;
        --smoke)       SMOKE=true; shift ;;
        --steps)       STEPS="${2:?--steps requires a list}"; shift 2 ;;
        --steps=*)     STEPS="${1#*=}"; shift ;;
        --skip-boltz)  STEPS="prep,dock,plip"; shift ;;
        -h|--help)     usage; exit 0 ;;
        *)
            echo "ERROR: unknown option '$1'" >&2
            echo "" >&2
            usage >&2
            exit 2
            ;;
    esac
done

for s in ${STEPS//,/ }; do
    case "$s" in
        prep|dock|plip|boltz|analyze) ;;
        *)
            echo "ERROR: unknown step '$s'. Known steps: prep, dock, plip, boltz, analyze" >&2
            exit 2
            ;;
    esac
done
has_step() { case ",${STEPS}," in *",$1,"*) return 0 ;; *) return 1 ;; esac; }

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: config not found: $CONFIG" >&2
    if [ -f "${ROOT}/constraints/config.example.yaml" ]; then
        echo "       Start from the example:" >&2
        echo "         cp constraints/config.example.yaml constraints/config.yaml" >&2
    fi
    exit 1
fi

read_cfg() {
    python3 - "$CONFIG" "$1" <<'PY'
import sys, yaml
with open(sys.argv[1]) as f:
    cfg = yaml.safe_load(f)
v = cfg
for k in sys.argv[2].split("."):
    v = v[k]
print(v)
PY
}

BASE_DIR="${ROOT}/$(read_cfg output.base_dir)"
INPUT_CSV=$(read_cfg input.csv_path)
case "$INPUT_CSV" in /*) ;; *) INPUT_CSV="${ROOT}/${INPUT_CSV}" ;; esac
ID_COL=$(read_cfg input.id_column)
CONF_THRESHOLD=$(read_cfg confidence_threshold)
RANK_BY=$(read_cfg rank_by)
TOP_N=$(read_cfg top_n)

# Prefix for the inference-parameter keys: "" (full run) or "smoke_"
P=""
if [ "$SMOKE" = true ]; then
    P="smoke_"
    SAMPLES=$(read_cfg smoke_samples_per_complex)
    N_LIGANDS=$(read_cfg smoke_n_ligands)
else
    SAMPLES=$(read_cfg diffdock.samples_per_complex)
    N_LIGANDS=""
fi
RECYCLING=$(read_cfg "${P}recycling_steps")
SAMPLING=$(read_cfg "${P}sampling_steps")
DIFF_SAMPLES_AFF=$(read_cfg "${P}diffusion_samples_affinity")
SAMPLING_AFF=$(read_cfg "${P}sampling_steps_affinity")
NO_KERNELS=$(read_cfg no_kernels)

LOG_DIR="${ROOT}/logging"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/constraints_$(date +%Y%m%d_%H%M%S).log"
log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG_FILE"; }

START_TIME=$(date +%s)
log "=== Constraint re-ranking started ==="
log "Config: $CONFIG"
log "Steps:  $STEPS"
[ "$SMOKE" = true ] && log "SMOKE mode: ${N_LIGANDS} ligands, ${SAMPLES} poses, ${RECYCLING}/${SAMPLING}/${DIFF_SAMPLES_AFF}/${SAMPLING_AFF}"

if [ ! -f "$INPUT_CSV" ]; then
    log "ERROR: hit list not found: $INPUT_CSV"
    log "       This repository ships no hit list. Point input.csv_path at yours."
    exit 1
fi

# In smoke mode, work on a truncated copy so the config's input stays untouched.
WORK_CSV="$INPUT_CSV"
if [ "$SMOKE" = true ]; then
    mkdir -p "$BASE_DIR"
    WORK_CSV="${BASE_DIR}/smoke_input.csv"
    python3 - "$INPUT_CSV" "$WORK_CSV" "$N_LIGANDS" <<'PY'
import sys, pandas as pd
df = pd.read_csv(sys.argv[1]).head(int(sys.argv[3]))
df.to_csv(sys.argv[2], index=False)
print(f"smoke input: {len(df)} ligand(s) -> {sys.argv[2]}")
PY
fi

# --- Step 1: SMILES -> 3D SDF ---
if has_step prep; then
    log "--- prep: SMILES to 3D SDF ---"
    python3 "${CONSTRAINTS}/05_prepare_ligands_for_diffdock.py" \
        --config "$CONFIG" \
        --input "$WORK_CSV" \
        --output-dir "${BASE_DIR}/sdf" 2>&1 | tee -a "$LOG_FILE"
fi

# --- Step 2: DiffDock ---
if has_step dock; then
    log "--- dock: DiffDock (${SAMPLES} poses/ligand) ---"
    python3 "${CONSTRAINTS}/06_run_diffdock.py" \
        --config "$CONFIG" \
        --input "${WORK_CSV%.csv}_with_sdf.csv" \
        --output-dir "${BASE_DIR}/docking" \
        --samples "$SAMPLES" 2>&1 | tee -a "$LOG_FILE"
fi

# --- Step 3: PLIP -> constraint YAMLs ---
if has_step plip; then
    log "--- plip: contacts and constraint YAMLs ---"
    python3 "${CONSTRAINTS}/07_generate_constraint_yamls.py" \
        --config "$CONFIG" \
        --input "${WORK_CSV%.csv}_with_sdf_with_diffdock.csv" \
        --output-dir "${BASE_DIR}/yamls" \
        --plip-output "${BASE_DIR}/plip" 2>&1 | tee -a "$LOG_FILE"
fi

# --- Step 4: Boltz-2 over both YAML variants ---
if has_step boltz; then
    log "--- boltz: inference at ${RECYCLING}/${SAMPLING}/${DIFF_SAMPLES_AFF}/${SAMPLING_AFF} ---"

    BOLTZ="${BOLTZ_EXE:-}"
    if [ -z "$BOLTZ" ]; then
        if command -v boltz >/dev/null 2>&1; then
            BOLTZ="$(command -v boltz)"
        else
            log "ERROR: 'boltz' not found on the PATH and \$BOLTZ_EXE is not set."
            log "       export BOLTZ_EXE=/path/to/exp_env/bin/boltz"
            exit 1
        fi
    fi
    if [ ! -x "$BOLTZ" ]; then
        log "ERROR: \$BOLTZ_EXE is not an executable file: $BOLTZ"
        exit 1
    fi

    if [ -d "${BASE_DIR}/yamls" ]; then
        YAML_COUNT=$(find "${BASE_DIR}/yamls" -name '*.yaml' 2>/dev/null | wc -l)
    else
        YAML_COUNT=0
    fi
    if [ "$YAML_COUNT" -eq 0 ]; then
        log "ERROR: no YAMLs in ${BASE_DIR}/yamls — run the plip step first."
        exit 1
    fi
    log "Boltz: $BOLTZ  ($YAML_COUNT YAMLs)"

    EXTRA=()
    if [ "$NO_KERNELS" = "True" ] || [ "$NO_KERNELS" = "true" ]; then
        EXTRA+=(--no_kernels)
        log "cuEquivariance kernels disabled (no_kernels: true)"
    fi

    "$BOLTZ" predict "${BASE_DIR}/yamls" \
        --out_dir "${BASE_DIR}/boltz_output" \
        --recycling_steps "$RECYCLING" \
        --sampling_steps "$SAMPLING" \
        --diffusion_samples_affinity "$DIFF_SAMPLES_AFF" \
        --sampling_steps_affinity "$SAMPLING_AFF" \
        "${EXTRA[@]}" 2>&1 | tee -a "$LOG_FILE"
fi

# --- Step 5: compare the two rankings ---
if has_step analyze; then
    log "--- analyze: constrained vs. unconstrained (QC >= ${CONF_THRESHOLD}) ---"
    python3 "${CONSTRAINTS}/08_analyze_constraint_predictions.py" \
        --config "$CONFIG" \
        --input "$WORK_CSV" \
        --results-dir "${BASE_DIR}/boltz_output" \
        --out-dir "${BASE_DIR}/analysis" \
        --confidence-threshold "$CONF_THRESHOLD" \
        --rank-by "$RANK_BY" \
        --top-n "$TOP_N" 2>&1 | tee -a "$LOG_FILE"
fi

END_TIME=$(date +%s)
log ""
log "=== Constraint re-ranking completed in $((END_TIME - START_TIME))s ==="
log "Output: $BASE_DIR"
log "Log:    $LOG_FILE"
if [ "$SMOKE" = true ]; then
    log ""
    log "This was a smoke run. Reduced sampling degrades pose quality badly;"
    log "do not read pose or ranking conclusions out of it."
fi
