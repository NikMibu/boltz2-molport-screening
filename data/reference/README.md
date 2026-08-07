# Reference data

Public data only. Everything in this directory must be traceable to a public
source and reproducible from a documented query. No screening library and no
screening result belongs here.

## `known_ca2_inhibitors.csv`

1,472 compounds with reported activity against human carbonic anhydrase II,
used as the reference set for the novelty and heavy-atom analyses.

| | |
|---|---|
| Source | [ChEMBL](https://www.ebi.ac.uk/chembl/) target `CHEMBL205` (carbonic anhydrase II) |
| Retrieved via | ChEMBL REST API, `activity.json` and `molecule.json` |
| Activity types | IC50 and Ki, in nM |
| Threshold | activity <= 100,000 nM (100 µM) |
| Limit | 1,000 results per activity type, deduplicated by compound keeping the best value |
| Retrieved | 14 July 2026 |
| Licence | ChEMBL data is released under CC BY-SA 3.0 |

Columns: `molecule_chembl_id`, `activity_value_nm`, `activity_type`,
`relation`, `smiles`, `heavy_atoms`.

**On the threshold.** 100 µM is permissive on purpose. The set defines "has
been reported to bind hCA II", not "is a potent inhibitor". A stricter cutoff
would shrink the reference set and make any hit list look more novel than it
is — the comparison should not be flattered by its own baseline.

**Reproducing or updating the file.** ChEMBL grows, so a fresh query returns a
different set and the novelty numbers shift accordingly. The version checked in
here is the one behind the figures in the thesis. To rebuild it:

```bash
python analysis/heavy_atom_comparison.py \
    --hits <your_hits.csv> \
    --refetch \
    --out-dir <output_dir>
```

This writes a new `known_ca2_inhibitors.csv` into the output directory rather
than overwriting the file here.
