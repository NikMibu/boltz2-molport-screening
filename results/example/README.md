# Example outputs

Produced from [`data/demo/`](../../data/demo/README.md) — 100 arbitrary vendor
catalogue compounds — so that the input and output formats are visible without
running anything. **No screening result of the thesis is reproduced here**, and
nothing in this folder says anything about which compounds bind carbonic
anhydrase.

## Present

| | Produced by |
|---|---|
| `yaml_example/*.yaml` | `screening/01_generate_yamls.py` — two Boltz-2 inputs, one per compound |
| `demo_100_properties_distributions.png` | `screening/00_characterize_dataset.py --plot` |

### The Boltz-2 input format

One YAML per compound. The protein is identical across the file set; only the
ligand block changes.

```yaml
version: 1

sequences:
  - protein:
      id: [A]
      sequence: MSHHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKY...
      msa: /abs/path/to/repo/data/msa/P00918.a3m

  - ligand:
      id: ["L1"]
      smiles: "COc1cc(O)c2C(=O)CC(Oc2c1)c1ccc(O)cc1"

properties:
  - affinity:
      binder: "L1"
```

Three things matter here and are easy to get wrong:

- **`msa:` is an absolute path.** Boltz-2 resolves it as written, so a YAML set
  moved between machines has to be repointed — that is what
  `02_update_msa.py` is for. Without the field, Boltz-2 falls back to its MSA
  server and makes one round trip per compound.
- **The ligand `id` is a running counter (`L1`, `L2`, …), not the compound ID.**
  Boltz truncates chain names at five characters. The compound ID lives in the
  *filename*, which is what `03_make_boltz_predictions_summary.py` matches on.
- **The `affinity` property is what produces a prediction at all.** Drop it and
  Boltz-2 returns a structure and no affinity.

### The constraint variant

`constraints/07_generate_constraint_yamls.py` writes two YAMLs per compound: the
form above, and the same file with a pocket constraint appended.

```yaml
constraints:
  - pocket:
      binder: "LIG"
      contacts: [["A", 199], ["A", 121], ["A", 200], ["A", 198], ["A", 131]]
      max_distance: 6.0
      force: false
```

The contact list comes from PLIP, run on the top-ranked DiffDock pose at PLIP's
default binding-site radius of 7.5 Å. `max_distance` is a *Boltz-2* parameter —
how far the binder may sit from those residues — and is not a PLIP setting.

No constraint example is checked in: producing one needs DiffDock and PLIP, and
the poses would be specific to whatever hit list was fed in.

## Not present yet

The ranking outputs need a GPU pass over the demo set:

| | Produced by |
|---|---|
| `demo_ranked_full.csv`, `demo_ranked_top100.csv` | `04_analyze_ligand_predictions.py` |
| `boltz_predictions_distributions.png` | same |
| `constraint_comparison.png`, `old_vs_new_comparison.png` | `08_analyze_constraint_predictions.py` |

```bash
bash run_screening.sh --target ca2 --smoke   # 5 compounds, minutes
bash run_screening.sh --target ca2           # 100 compounds at 2/50/1/50
```

See [`SMOKE_TEST.md`](../../SMOKE_TEST.md) for the full sequence.
