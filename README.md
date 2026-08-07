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

The screening and re-ranking pipeline (scripts `01`–`14`) is being migrated in.
Currently available: the downstream analyses.

| | |
|---|---|
| [`analysis/`](analysis/) | Novelty and property comparison of a hit list against known hCA II inhibitors |
| `screening/` | *to follow* — YAML generation, Boltz-2 inference, ranking |
| `constraints/` | *to follow* — DiffDock poses, PLIP contacts, constrained re-scoring |

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
