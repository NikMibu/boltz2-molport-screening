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

## Four of the 100 are salts, on purpose

`demo_100.csv` contains four multi-fragment SMILES — 4 %, against 3.9 % in the
full catalogue, so the sample is representative there too.

**Boltz-2 2.2.1 does not predict them.** It reduces a multi-fragment input to
its largest fragment, and RDKit 2025.9.2 aborts partway through:

```
RuntimeError: Pre-condition Violation
        getNumImplicitHs() called without preceding call to calcImplicitValence()
```

Boltz reports `Failed to process … Skipping.` and continues, so the compound
does not appear in the output at all and quietly vanishes from the ranking.

This is not a demo-set artefact. In the screening run behind this pipeline,
**all 306 multi-fragment compounds that passed the size filter came back
without a prediction — every single one.** They account for 71 % of the 432
compounds that the run lost between the 8,078 that passed filtering and the
7,646 that produced predictions.

They are left in deliberately: a demo set that quietly avoided the one input
class the pipeline mishandles would be misleading. `01_generate_yamls.py` names
them in a warning before the GPU run rather than after it. If you need salts
predicted, desalt them first — the pipeline does not do it for you, because
choosing which fragment to keep is a chemistry decision.
