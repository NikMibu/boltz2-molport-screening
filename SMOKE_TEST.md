# Smoke test

**Run it from a fresh clone, not from a working copy.** The point is to test
what is actually on GitHub. A working copy quietly contains files that were
never committed — that is how three sets of results in this project turned out
to exist only locally, and how a `set -e` bug survived every test until the
first clone hit it.

## Four runs, two purposes

Each pipeline has a smoke variant and a full one. They differ by two orders of
magnitude in cost. Do the smoke variant before committing to the full one.

| | | Compounds | Parameters | Wall clock on an RTX 4060 Ti |
|---|---|---|---|---|
| Screening | `--smoke` | 5 per target | 1/10/1/10 | minutes |
| Screening | full | 100 per target | 2/50/1/50 | ~1.5 h per isoform, **~5 h for all three** |
| Constraints | `--smoke` | 3, 4 poses each | 1/10/1/10 | minutes, plus the first-run table build |
| Constraints | full | 10, 40 poses each | 3/200/5/200 | **2–3 h** |

The smoke runs answer "does the chain execute". Only the full runs produce
output worth putting in `results/example/`, because at 1/10 the numbers mean
nothing.

Note the constraint stage costs far more per compound: every compound is
predicted twice, at four times the sampling of the screening stage. Ten
compounds there are comparable to a hundred here.

## 0. Clone and environment

On WSL, clone to a **native path**. A repository under `/mnt/c` or `/mnt/e` goes
through the Windows filesystem layer, and Boltz-2's checkpoint loading and
DiffDock's dataset caching are both punished heavily for it.

```bash
cd ~ && git clone https://github.com/NikMibu/boltz2-molport-screening.git
cd boltz2-molport-screening

python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Point at the three external tools. None of them is reliably on the `PATH`:

```bash
export BOLTZ_EXE=/path/to/boltz-env/bin/boltz
export MAMBA_EXE=/path/to/bin/micromamba
export DIFFDOCK_HOME=/path/to/DiffDock
```

`micromamba` in particular: the standard install defines a shell *function* and
leaves the binary off the `PATH`, so it works when typed and is invisible to
`subprocess`. `$MAMBA_EXE` is what the shell hook exports and what the code
looks for.

**On WSL2, set `no_kernels: true`** in `screening/screening_config.yaml` and
`constraints/config.yaml`. Boltz's cuEquivariance kernels call
`nvmlDeviceGetNumGpuCores`, which WSL's NVML does not implement, and prediction
crashes. Slower, same results. Not needed on native Linux.

## 1. Validate

```bash
python3 validate_setup.py
```

Everything under **Core** and **Screening workflow** must pass. Constraint
warnings are fine if you only want the screening half.

Watch the `FASTA/MSA pair` lines: they check that each target's alignment was
built from the sequence it is paired with. hCA IV ships in a trimmed and an
untrimmed form under the same UniProt header, and Boltz-2 accepts a mismatched
pair without complaint.

## 2. Screening

Start small:

```bash
bash run_screening.sh --target ca2 --smoke
```

Five compounds at 1/10/1/10, a few minutes. If that produces a ranking, the
wiring is sound and the long run is worth starting.

```bash
bash run_screening.sh --target all      # ~5 h, see the table above
```

Expect five of the 100 compounds to be skipped per isoform with
`Failed to process … Skipping.` — four are the multi-fragment SMILES that
Boltz-2 2.2.1 cannot reduce, which `01_generate_yamls.py` names at the start of
the run; the fifth is a large polyketide that fails for another reason. This is
expected; see [`data/demo/`](data/demo/README.md).

If the run is interrupted, the ranking can be redone from whatever predictions
exist without re-running inference:

```bash
bash run_screening.sh --target ca2 --analysis-only
```

**Check, per isoform:**

- `results/<t>/predictions_<t>.csv` — 100 rows, `boltz_prob_mean` and
  `boltz_pred_value_mean` populated for 95 of them
- step 4 reports the skipped compounds as folders without `affinity_*.json`,
  and **no unmatched folders**
- `results/<t>/ranking/ligands_ranked_full.csv` — `passes_confidence` true only
  at or above 0.6, failures carry `combined_score` 999, sorted ascending
- `results/<t>/ranking/boltz_predictions_distributions.png` exists

Then the isoform comparison:

```bash
python3 screening/analyze_isoform_selectivity.py \
    --ca2 results/ca2/ranking/ligands_ranked_full.csv \
    --ca4 results/ca4/ranking/ligands_ranked_full.csv \
    --ca7 results/ca7/ranking/ligands_ranked_full.csv \
    --out-dir results/isoform_selectivity
```

## 3. Constraints

```bash
cp constraints/config.example.yaml constraints/config.yaml
sed -i 's/^no_kernels: false/no_kernels: true/' constraints/config.yaml   # WSL only
bash run_constraints.sh --smoke
```

One command, five steps: 3D ligands, DiffDock, PLIP, Boltz-2 over both YAML
variants, comparison. `--smoke` cuts it to 3 compounds, 4 poses each, inference
at 1/10/1/10 — minutes once the lookup tables are built.

Then the real demo run:

```bash
bash run_constraints.sh
```

10 compounds, 40 poses each, 20 YAMLs at 3/200/5/200. **Budget 2-3 hours.**
Per compound this stage costs an order of magnitude more than screening: twice
the predictions, at four times the sampling.

Individual steps can be re-run without repeating the expensive ones:

```bash
bash run_constraints.sh --steps analyze
bash run_constraints.sh --steps plip,boltz,analyze
```

**The first DiffDock run builds its SO(3) lookup tables. That takes 5–10 minutes
and produces no output — it is not hung.**

**Check:**

- `results/constraints/ca2/sdf/` — one SDF per compound
- `results/constraints/ca2/docking/*/complex_0/rank1.sdf` — a pose per compound
- `results/constraints/ca2/plip/*/report.xml` — a PLIP report per compound
- `results/constraints/ca2/yamls/` — **two** YAMLs per compound, `_default` and
  `_pocket`, and the pocket ones carry a non-empty `contacts:` list
- `results/constraints/ca2/analysis/` — comparison plots and a ranking

An empty `contacts:` list is not automatically a bug — a pose can genuinely make
no detected contact — but if *every* compound comes out empty, PLIP is failing
silently. Check one `report.xml` by hand.

## 4. What a smoke run does not tell you

At 1/10 sampling the structures are poor. In the sister repository, median pose
RMSD went from 0.48 Å at full parameters to roughly 18 Å at 1/10 — enough to
conclude the model cannot place a ligand, which is not true. Affinity rankings
degrade more gracefully but are still not comparable across parameter settings.

A smoke run answers one question: does the chain execute without dropping data
on the floor. Nothing else.
