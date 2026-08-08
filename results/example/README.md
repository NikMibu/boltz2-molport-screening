# Example outputs

Produced by running both pipelines over [`data/demo/`](../../data/demo/README.md)
— 100 arbitrary vendor catalogue compounds — so that the formats and the shape
of the output are visible without running anything.

**These are not the results of the thesis.** The compounds were not selected,
the target was not chosen for them, and nothing here says anything about which
molecules bind carbonic anhydrase. What the numbers demonstrate is that the
pipeline runs and that its filters do something.

## The run behind these files

| | |
|---|---|
| Input | `data/demo/demo_100.csv` — 100 compounds; `demo_10.csv` for the constraint stage |
| Screening | 2/50/1/50, QC threshold 0.6, per isoform |
| Constraints | 40 DiffDock poses per ligand, then 3/200/5/200 over both YAML variants |
| Hardware | single NVIDIA RTX 4060 Ti (16 GB), WSL2, `no_kernels: true` |
| Wall clock | ~5 h for the three isoforms, ~35 min for the constraint stage on 10 compounds |
| Versions | Boltz-2 2.2.1, PLIP 2.3.1, RDKit 2025.9.2 |

## Files

| | Produced by |
|---|---|
| `demo_ca2_ranked_full.csv`, `demo_ca4_…`, `demo_ca7_…` | `04_analyze_ligand_predictions.py` |
| `ca2_predictions_distributions.png`, `ca4_…`, `ca7_…` | same |
| `isoform_selectivity_analysis.png`, `isoform_overlap_bar.png` | `analyze_isoform_selectivity.py` |
| `constraint_comparison.png` | `08_analyze_constraint_predictions.py` |
| `demo_100_properties_distributions.png` | `00_characterize_dataset.py --plot` |
| `yaml_example/*.yaml` | `01_generate_yamls.py` |

There is no `*_ranked_top100.csv`: the demo list is 100 compounds long, so the
top 100 is the same file.

## What the numbers show

**95 of 100 compounds produced a prediction**, identically in all three
isoforms. The five that did not are listed in the CSVs with empty prediction
columns:

| | |
|---|---|
| 4 compounds | multi-fragment SMILES (salts) — Boltz-2 2.2.1 cannot reduce them, see [`data/demo/`](../../data/demo/README.md) |
| `MolPort-001-739-307` | a single-fragment polyketide that failed for another reason |

**Ten, nine and eleven compounds pass the confidence QC** for hCA II, hCA IV and
hCA VII. That is the point of the threshold: an untargeted sample of natural
products screened against one enzyme should mostly not clear it. The highest
binding probability in the hCA II set is 0.714, and the 90 compounds below 0.6
carry the sentinel score 999 regardless of their predicted affinity — including
one whose predicted affinity is the best in the whole set.

**The isoform overlap is 100 % and means nothing here.** Comparing the top 100
of three 100-compound lists returns every compound three times by construction.
`analyze_isoform_selectivity.py` says so before printing the number. The
selectivity statistics rest on the 8 compounds that clear QC in more than one
isoform, which is too few to read anything into. Demonstrating that analysis
needs a real library; these files only show that it runs and what it emits.

## The Boltz-2 input format

One YAML per compound, in `yaml_example/`. The protein block is identical across
the set; only the ligand changes.

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

Three things here are easy to get wrong:

- **`msa:` is an absolute path.** Boltz-2 resolves it as written, so a YAML set
  moved between machines has to be repointed — that is what `02_update_msa.py`
  is for. Without the field, Boltz-2 falls back to its MSA server and makes one
  round trip per compound.
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
default binding-site radius. `max_distance` is a *Boltz-2* parameter — how far
the binder may sit from those residues — and is not a PLIP setting.

No constraint example is checked in: the poses would be specific to whichever
hit list was fed in.
