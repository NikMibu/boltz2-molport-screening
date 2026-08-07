#!/usr/bin/env bash
# Prospective virtual screening with Boltz-2 (Chapter 6.4).
#   bash run_screening.sh --target ca2
#   bash run_screening.sh --target all --input data/demo/demo_100.csv
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
CONFIG="${ROOT}/screening/screening_config.yaml"
SCREENING="${ROOT}/screening"

usage() {
    cat <<EOF
Usage: bash run_screening.sh [OPTIONS]

Screens a compound table against one or more carbonic anhydrase isoforms:
YAML generation -> MSA injection -> Boltz-2 inference -> ranking.

Options:
  --target KEY      ca2, ca4, ca7, or all (default: ca2)
  --input CSV       Compound table; overrides input_csv from the config.
                    Needs a SMILES column and an ID column.
  --config PATH     Pipeline config (default: screening/screening_config.yaml)
  --limit N         Only generate the first N YAMLs. For trying the wiring.
  --smoke           Reduced run: the first smoke_n_compounds of the input at
                    reduced inference parameters (see the config). Minutes
                    instead of hours. Proves the chain runs; reproduces nothing.
  --skip-boltz      Generate YAMLs but stop before inference. Use when the GPU
                    run happens elsewhere.
  --analysis-only   Skip generation and inference, re-run ranking on existing
                    Boltz output.
  -h, --help        Show this message and exit

Examples:
  bash run_screening.sh --target ca2 --smoke
  bash run_screening.sh --target all
  bash run_screening.sh --target ca2 --analysis-only

Inference parameters come from the config and are 2/50/1/50, the reduced
"baseline" settings used for the prospective screen — not the Boltz-2 defaults.
The QC threshold is 0.6, not the 0.5 that 04_analyze_ligand_predictions.py
defaults to. Both are documented in the config; changing them changes the
rankings.

Runtime: the full library is 8,963 compounds per isoform. Budget days on a
single GPU, and start with --limit.

Boltz-2 is taken from \$BOLTZ_EXE if set, otherwise 'boltz' on the PATH.
EOF
}

TARGET="ca2"
INPUT_OVERRIDE=""
LIMIT=""
SMOKE=false
SKIP_BOLTZ=false
ANALYSIS_ONLY=false

while [ $# -gt 0 ]; do
    case "$1" in
        --target)         TARGET="${2:?--target requires a value}"; shift 2 ;;
        --target=*)       TARGET="${1#*=}"; shift ;;
        --input)          INPUT_OVERRIDE="${2:?--input requires a path}"; shift 2 ;;
        --input=*)        INPUT_OVERRIDE="${1#*=}"; shift ;;
        --config)         CONFIG="${2:?--config requires a path}"; shift 2 ;;
        --config=*)       CONFIG="${1#*=}"; shift ;;
        --limit)          LIMIT="${2:?--limit requires a number}"; shift 2 ;;
        --limit=*)        LIMIT="${1#*=}"; shift ;;
        --smoke)          SMOKE=true; shift ;;
        --skip-boltz)     SKIP_BOLTZ=true; shift ;;
        --analysis-only)  ANALYSIS_ONLY=true; shift ;;
        -h|--help)        usage; exit 0 ;;
        *)
            echo "ERROR: unknown option '$1'" >&2
            echo "" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: config not found: $CONFIG" >&2
    exit 1
fi

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

# --- Resolve targets ---
ALL_TARGETS=$(python3 - "$CONFIG" <<'PY'
import sys, yaml
with open(sys.argv[1]) as f:
    cfg = yaml.safe_load(f)
print(" ".join(cfg["targets"].keys()))
PY
)

if [ "$TARGET" = "all" ]; then
    TARGETS="$ALL_TARGETS"
else
    case " $ALL_TARGETS " in
        *" $TARGET "*) TARGETS="$TARGET" ;;
        *)
            echo "ERROR: unknown target '$TARGET'. Config defines: $ALL_TARGETS (or 'all')" >&2
            exit 2
            ;;
    esac
fi

INPUT_CSV="${INPUT_OVERRIDE:-$(read_cfg input_csv)}"
case "$INPUT_CSV" in
    /*) ;;
    *) INPUT_CSV="${ROOT}/${INPUT_CSV}" ;;
esac
if [ "$ANALYSIS_ONLY" = false ] && [ ! -f "$INPUT_CSV" ]; then
    echo "ERROR: compound table not found: $INPUT_CSV" >&2
    echo "       This repository ships no screening library. Pass your own with" >&2
    echo "       --input, or build the demo set first (see data/demo/README.md)." >&2
    exit 1
fi

# --smoke caps the compound count too; an explicit --limit still wins.
if [ "$SMOKE" = true ] && [ -z "$LIMIT" ]; then
    LIMIT=$(read_cfg smoke_n_compounds)
fi

SMILES_COL=$(read_cfg smiles_column)
ID_COL=$(read_cfg id_column)
CONF_THRESHOLD=$(read_cfg confidence_threshold)
TOP_N=$(read_cfg top_n)
MAX_HEAVY=$(read_cfg max_heavy_atoms)
MAX_MW=$(read_cfg max_mw)

LOG_DIR="${ROOT}/logging"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/screening_$(date +%Y%m%d_%H%M%S).log"
log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG_FILE"; }

START_TIME=$(date +%s)
log "=== Boltz-2 screening started ==="
log "Config:  $CONFIG"
log "Targets: $TARGETS"
[ "$ANALYSIS_ONLY" = false ] && log "Input:   $INPUT_CSV"

for T in $TARGETS; do
    log ""
    log "########## Target: $T ##########"

    FASTA="${ROOT}/$(read_cfg "targets.${T}.fasta")"
    MSA="${ROOT}/$(read_cfg "targets.${T}.msa")"
    YAML_DIR="${ROOT}/$(read_cfg yaml_dir "$T")"
    RESULTS_DIR="${ROOT}/$(read_cfg results_dir "$T")"
    ANALYSIS_DIR="${ROOT}/$(read_cfg analysis_dir "$T")"

    if [ "$ANALYSIS_ONLY" = false ]; then
        for f in "$FASTA" "$MSA"; do
            if [ ! -f "$f" ]; then
                log "ERROR: missing input file: $f"
                exit 1
            fi
        done

        # --- Step 1: Generate YAMLs ---
        log "--- Step 1: YAML generation ---"
        # Nothing clears the folder, and boltz predicts everything it finds.
        # YAMLs left over from an earlier input would be predicted alongside
        # the current one and then show up as unmatched in step 4.
        STALE_BEFORE=$(find "$YAML_DIR" -name '*.yaml' 2>/dev/null | wc -l)
        YAML_CMD=(python3 "${SCREENING}/01_generate_yamls.py"
                  --smiles-csv "$INPUT_CSV"
                  --smiles-column "$SMILES_COL"
                  --id-column "$ID_COL"
                  --fasta "$FASTA"
                  --msa-path "$MSA"
                  --outdir "$YAML_DIR")
        [ -n "$LIMIT" ] && YAML_CMD+=(--limit "$LIMIT")
        "${YAML_CMD[@]}" 2>&1 | tee -a "$LOG_FILE"

        # --- Step 2: MSA paths ---
        # The alignments are pre-computed locally rather than fetched per run
        # with --use_msa_server. YAMLs hold an absolute path to the .a3m, so
        # this step brings the whole folder up to date with wherever the
        # alignment currently lives — including YAMLs generated on another
        # machine, which is the normal case for a server run.
        log "--- Step 2: MSA paths ---"
        python3 "${SCREENING}/02_update_msa.py" "$YAML_DIR" "$MSA" 2>&1 | tee -a "$LOG_FILE"

        YAML_COUNT=$(find "$YAML_DIR" -name '*.yaml' | wc -l)
        log "$YAML_COUNT YAML(s) ready"

        EXPECTED=$(python3 -c "
import csv, sys
n = sum(1 for _ in csv.DictReader(open('$INPUT_CSV', newline='')))
print(min(n, ${LIMIT:-n}) if '${LIMIT}' else n)
" 2>/dev/null || echo "")
        if [ -n "$EXPECTED" ] && [ "$YAML_COUNT" -gt "$EXPECTED" ]; then
            log "WARNING: $YAML_DIR holds $YAML_COUNT YAMLs but the input has $EXPECTED"
            log "         compound(s). $((YAML_COUNT - EXPECTED)) file(s) are left over from an"
            log "         earlier run and will be predicted too. Clear the folder first."
        fi

        # --- Step 3: Boltz-2 inference ---
        if [ "$SKIP_BOLTZ" = true ]; then
            log "--- Step 3: SKIPPED (--skip-boltz) ---"
            continue
        fi
        log "--- Step 3: Boltz-2 inference ---"
        SMOKE_ARG=""
        [ "$SMOKE" = true ] && SMOKE_ARG="--smoke"
        bash "${SCREENING}/run_boltz.sh" "$CONFIG" "$T" $SMOKE_ARG 2>&1 | tee -a "$LOG_FILE"
    fi

    # --- Step 4: Collect predictions ---
    log "--- Step 4: Collecting predictions ---"
    mkdir -p "$ANALYSIS_DIR"
    python3 "${SCREENING}/03_make_boltz_predictions_summary.py" \
        --results-dir "$RESULTS_DIR" \
        --input "$INPUT_CSV" \
        --id-column "$ID_COL" \
        --out-dir "$ANALYSIS_DIR" \
        --out-name "predictions_${T}.csv" 2>&1 | tee -a "$LOG_FILE"

    # --- Step 5: Ranking ---
    log "--- Step 5: Ranking (QC >= ${CONF_THRESHOLD}, top ${TOP_N}) ---"
    python3 "${SCREENING}/04_analyze_ligand_predictions.py" \
        --input "${ANALYSIS_DIR}/predictions_${T}.csv" \
        --out-dir "${ANALYSIS_DIR}/ranking" \
        --confidence-threshold "$CONF_THRESHOLD" \
        --top-n "$TOP_N" \
        --max-heavy-atoms "$MAX_HEAVY" \
        --max-mw "$MAX_MW" 2>&1 | tee -a "$LOG_FILE"

    log "Target $T done -> ${ANALYSIS_DIR}/ranking"
done

END_TIME=$(date +%s)
log ""
log "=== Screening completed in $((END_TIME - START_TIME))s ==="
log "Log: $LOG_FILE"
if [ "$(echo "$TARGETS" | wc -w)" -gt 1 ]; then
    log ""
    log "All three isoforms ranked. For the selectivity comparison (Chapter 6.4.6):"
    log "  python3 screening/analyze_isoform_selectivity.py \\"
    log "      --ca2 results/ca2/ranking/ligands_ranked_full.csv \\"
    log "      --ca4 results/ca4/ranking/ligands_ranked_full.csv \\"
    log "      --ca7 results/ca7/ranking/ligands_ranked_full.csv \\"
    log "      --out-dir results/isoform_selectivity"
fi
