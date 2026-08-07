# Smoke test

Exercises both pipelines end to end on the 100-compound demo set. Roughly an
hour of wall clock on a single consumer GPU, most of it DiffDock.

**Run it from a fresh clone, not from a working copy.** The point is to test
what is actually on GitHub. A working copy quietly contains files that were
never committed — that is how three sets of results in this project turned out
to exist only locally.

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
export BOLTZ_EXE=/mnt/e/master_thesis/exp_env/bin/boltz
export MAMBA_EXE=/mnt/e/master_thesis/bin/micromamba
export DIFFDOCK_HOME=/mnt/e/master_thesis/DiffDock
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

## 2. Screening — minutes

```bash
bash run_screening.sh --target ca2 --smoke
```

Five compounds at 1/10/1/10. Expect a few minutes on an RTX 4060 Ti.

Then all three isoforms at the real parameters:

```bash
bash run_screening.sh --target all
```

100 compounds × 3 isoforms at 2/50/1/50. Budget around 30–60 minutes per
isoform on a 16 GB consumer card.

**Check:**

- `results/ca2/predictions_ca2.csv` — 100 rows, `boltz_prob_mean` and
  `boltz_pred_value_mean` populated
- `results/ca2/ranking/ligands_ranked_full.csv` — `passes_confidence` true only
  above 0.6, failures carry `combined_score` 999
- step 4 reports no unmatched prediction folders
- `results/*/ranking/boltz_predictions_distributions.png` exists

Then the isoform comparison:

```bash
python3 screening/analyze_isoform_selectivity.py \
    --ca2 results/ca2/ranking/ligands_ranked_full.csv \
    --ca4 results/ca4/ranking/ligands_ranked_full.csv \
    --ca7 results/ca7/ranking/ligands_ranked_full.csv \
    --out-dir results/isoform_selectivity
```

## 3. Constraints — the long half

```bash
cp constraints/config.example.yaml constraints/config.yaml
# edit input.csv_path -> ./data/demo/demo_10.csv
bash run_constraints.sh --smoke
```

10 compounds, 4 DiffDock poses each, inference at 1/10/1/10.

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

## 5. Filling `results/example/`

Once the full demo run is through, copy the outputs a reviewer should see
without running anything:

```bash
cp results/ca2/ranking/ligands_ranked_full.csv   results/example/demo_ranked_full.csv
cp results/ca2/ranking/ligands_ranked_top100.csv results/example/demo_ranked_top100.csv
cp results/ca2/ranking/boltz_predictions_distributions.png results/example/
cp results/constraints/ca2/analysis/constraint_comparison.png results/example/
cp results/constraints/ca2/analysis/old_vs_new_comparison.png results/example/
```

The `.gitignore` allows `results/example/README.md`, `demo_*.csv`, and any
`.png` or `.yaml` beneath it — and nothing else. Run `git status` before
committing and check that no ranking of a real hit list slipped in.
