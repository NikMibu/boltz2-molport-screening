# Demo set

**These are not screening results.** They are an arbitrary sample of input
compounds, shipped so that both pipelines can be run end to end without a
compound library of your own. Nothing here was selected, ranked or predicted.

| File | Compounds | Used by |
|---|---|---|
| `demo_100.csv` | 100 | `run_screening.sh` |
| `demo_10.csv` | 10, a subset of the above | `run_constraints.sh --smoke` |

## Provenance

Drawn from the MolPort natural-product catalogue (`Molport ID`, `SMILES`),
a commercial vendor catalogue of purchasable compounds. Selection:

1. every structure annotated with `screening/00_characterize_dataset.py`
   (RDKit 2026.03.5)
2. restricted to molecular weight ≤ 700 Da and ≤ 56 heavy atoms — Boltz-2's
   affinity module is not recommended above that size, and compounds outside
   it would only be dropped again during ranking
3. `pandas.DataFrame.sample(n=100, random_state=42)` over the 8,078 remaining
   compounds, sorted by MolPort ID
4. `demo_10.csv` is a second draw of 10 from those 100, same seed

Both draws are reproducible from the catalogue with the seed above.

## Columns

`molport_number`, `smiles`, then the properties from
`00_characterize_dataset.py`: `mw`, `heavy_atoms`, `num_atoms`, `logp`, `hbd`,
`hba`, `rotatable_bonds`, `pains`.

The property columns are what `04_analyze_ligand_predictions.py` filters on.
Two of them shift slightly between RDKit releases: recomputing the full
catalogue under RDKit 2026.03.5 instead of the 2025.9.2 used in the thesis
changes `hba` for 2.4 % of compounds and `rotatable_bonds` for one. `mw`,
`heavy_atoms`, `logp`, `hbd` and the PAINS flag are identical, and the size
filter admits the same 8,078 compounds either way.

## What the demo set is good for

Checking that the wiring works, seeing the input and output formats, and
producing a ranking whose *shape* is representative. It is not good for
anything else: 100 arbitrary catalogue entries screened against one target
produce a ranking of 100 arbitrary catalogue entries.

## Composition

Median molecular weight 377 Da, median 27 heavy atoms, 9 of the 100 carry a
PAINS flag — close to the 8.7 % of the full catalogue, as an untargeted sample
should be.
