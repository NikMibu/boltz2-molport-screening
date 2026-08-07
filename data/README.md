# Input data

Everything tracked here comes from a public source. The screening library and
all results derived from it are proprietary and are not part of this
repository — see the note at the top of [`../.gitignore`](../.gitignore).

## `fasta/` — target sequences (UniProt)

| File | UniProt | Residues | Used for |
|---|---|---|---|
| `P00918.fasta` | P00918 (CA2_HUMAN) | 260 | hCA II |
| `P22748.fasta` | P22748 (CAH4_HUMAN) | 312 | *provenance only, not used* |
| `P22748_mat.fasta` | P22748, trimmed | 266 | hCA IV |
| `P43166.fasta` | P43166 (CAH7_HUMAN) | 264 | hCA VII |

**hCA IV is screened against the mature sequence.** The UniProt entry carries an
N-terminal signal peptide and a C-terminal GPI-anchor signal; both are cleaved
in the mature protein, and both are trimmed here (312 → 266 residues) so that
prediction sees the membrane-anchored catalytic domain rather than two
disordered tails. `P22748_mat.fasta` keeps the original UniProt FASTA header, so
the header alone does not tell the two files apart — check the length.

## `msa/` — multiple sequence alignments

Generated from the sequences above with the Boltz-2 MSA server. One `.a3m` per
target, ~3–4 MB each.

| File | Target | Pairs with |
|---|---|---|
| `P00918.a3m` | hCA II | `fasta/P00918.fasta` |
| `P22748_mat.a3m` | hCA IV | `fasta/P22748_mat.fasta` |
| `P43166.a3m` | hCA VII | `fasta/P43166.fasta` |

An alignment must match the sequence it was built from. Pairing
`P22748_mat.a3m` with the untrimmed `P22748.fasta` produces YAMLs that Boltz-2
accepts and silently mispredicts.

There is no alignment for the untrimmed hCA IV sequence in this repository. One
was generated during the work (`P22748_mit_Signal.a3m`) and superseded; it is
deliberately untracked.

## `pdb/1CNC.pdb` — apo receptor (RCSB)

Human carbonic anhydrase II. Used as the DiffDock receptor in the constraint
re-ranking (`constraints/`), after stripping HETATM records. Not used by the
screening pipeline.

## `reference/` — known hCA II inhibitors (ChEMBL)

1,472 compounds with SMILES, the reference set for the novelty and heavy-atom
analyses in `analysis/`. Provenance and the activity cutoff are documented in
[`reference/README.md`](reference/README.md).

## `demo/` — public demo set

Not present yet. A ~100-compound set drawn from `reference/` with a fixed seed,
so that the pipeline can be run end to end without the proprietary library.
`screening/screening_config.yaml` points `input_csv` at it by default.
