# Boltz-2 MolPort Screening

Prospective virtual screening of a commercial natural-product library against
human carbonic anhydrase isoforms (CA II, CA IV, CA VII) with **Boltz-2**,
followed by constraint-based re-ranking with DiffDock and PLIP.

> **Code-only repository.** The screening library and every result derived from
> it — predictions, rankings, hit lists — are proprietary and are not published
> here. The code runs on any hit list with a `SMILES` column; see
> [Input formats](#input-formats). Public reference data lives in
> [`data/reference/`](data/reference/).

## Status

| | |
|---|---|
| [`screening/`](screening/) | YAML generation, MSA injection, Boltz-2 inference, ranking, isoform selectivity |
| [`analysis/`](analysis/) | Novelty and property comparison of a hit list against known hCA II inhibitors |
| [`data/`](data/README.md) | Target sequences, alignments, apo receptor, ChEMBL reference set |
| `constraints/` | *to follow* — DiffDock poses, PLIP contacts, constrained re-scoring |

## Screening

```bash
# one isoform, five compounds, no GPU — checks the wiring
bash run_screening.sh --target ca2 --input my_compounds.csv --limit 5 --skip-boltz

# all three isoforms, end to end
bash run_screening.sh --target all --input my_compounds.csv
```

The wrapper runs YAML generation → MSA injection → `boltz predict` → prediction
collection → ranking, per target. It answers `--help`, takes `--config`, and
exits 2 on an unknown flag. Boltz-2 is taken from `$BOLTZ_EXE` if set, otherwise
`boltz` on the `PATH`.

**The input is yours to supply.** No compound library ships with this
repository. Any CSV with a SMILES column and an ID column works; the columns
are named in the config or auto-detected.

### Scoring

Two stages, both configurable in
[`screening/screening_config.yaml`](screening/screening_config.yaml):

1. **Quality control** — a compound must reach `boltz_prob_mean >= 0.6`.
   Compounds below it receive the sentinel score 999 and drop out of the
   ranking. A strong predicted affinity does not rescue a low-confidence
   prediction.
2. **Ranking** — `combined_score = boltz_pred_value_mean - 0.1 * boltz_prob_mean`,
   **lower is better**.

Before either stage, compounds above 56 heavy atoms or 700 Da are filtered out:
Boltz-2's affinity module is not recommended beyond that size.

> The argparse default inside `04_analyze_ligand_predictions.py` is `0.5`, not
> `0.6`. The wrapper always passes the configured value. Running the script
> directly without `--confidence-threshold` will not reproduce the published
> counts.

### Inference parameters

The prospective screen uses **2/50/1/50** — recycling steps, sampling steps,
affinity recycling steps, affinity sampling steps. These are *reduced* settings,
not the Boltz-2 defaults of 3/200/5/200: a library of ~9,000 compounds against
three isoforms is not affordable at the defaults. The constraint re-ranking of
the top 100 does use the defaults.

Reduced sampling costs pose quality. That is acceptable when the output is an
affinity ranking, and not acceptable when the output is an RMSD — see the
[sister repository](https://github.com/NikMibu/plip_constraints_pipeline) for
what happens to pose accuracy at low sampling.

### Targets

| Key | Isoform | Sequence |
|---|---|---|
| `ca2` | hCA II | P00918, 260 aa |
| `ca4` | hCA IV | P22748 **mature**, 266 aa — signal peptide and GPI-anchor signal trimmed |
| `ca7` | hCA VII | P43166, 264 aa |

Once all three are ranked, `screening/analyze_isoform_selectivity.py` compares
the rank lists; `screening/14_plot_boltz_structure_quality.py` compares the
structural confidence metrics across targets. Both take file paths and ship no
data of their own.

## Analyses

All three take the hit list as an argument. None of them ship one.

```bash
# Structural novelty: max Tanimoto against known hCA II inhibitors
python analysis/novelty_check.py --hits my_hits.csv --out-dir results/

# Exact SMILES matches against the same reference set
python analysis/smiles_overlap_check.py --hits my_hits.csv

# Heavy-atom distribution: hits vs. known inhibitors vs. full library
python analysis/heavy_atom_comparison.py --hits my_hits.csv --out-dir results/
```

`novelty_check.py` reports how many hits fall below a Tanimoto of 0.3 to every
known inhibitor — the operational definition of "structurally novel" used
throughout. `heavy_atom_comparison.py` checks that the hits occupy the same
size range as known binders; a hit list that drifts to much larger molecules is
usually an artefact of the scoring function rather than a discovery.

Every script answers `--help`.

## Input formats

**Hit list** (`--hits` / `--your-hits-csv`) — CSV with at least a `smiles` or
`SMILES` column. A score column is picked up automatically if one of
`boltz_pred_value_mean`, `pred_score`, `prediction_score` is present, or can be
named with `--score-col`; it is only used for the novelty-vs-score plot.

**Reference set** (`--known`) — CSV with a `smiles` column. Defaults to
[`data/reference/known_ca2_inhibitors.csv`](data/reference/README.md).

## Requirements

Python 3.10, RDKit, pandas, numpy, matplotlib, seaborn, requests. No GPU for
the analyses; the screening pipeline needs one.

## Related repositories

| | |
|---|---|
| [`boltz2-dude-benchmark`](https://github.com/NikMibu/boltz2-dude-benchmark) | Retrospective validation of Boltz-2 affinity prediction on DUD-E |
| [`plip_constraints_pipeline`](https://github.com/NikMibu/plip_constraints_pipeline) | Interaction-constrained co-folding, benchmarked against crystal poses |

## Licence

MIT, see [LICENSE](LICENSE). ChEMBL reference data under CC BY-SA 3.0.
