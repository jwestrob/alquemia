# Stage B: pinned sources and prepared inputs

Approved by the agreement in `AGREEMENT.md`. Preparation began after all five
Stage A preparations were ready; the exact Stage A inventory is pinned in the
Stage B inventory. No electronic-structure jobs were submitted by this builder.

## Result

| Source | Deposited Ca selector | Qualifying O donors | Donor residues | Waters | v3 / extended atoms |
|---|---|---:|---:|---:|---:|
| 2FW0 | chain A, CA310/CA | 7 | 6 | 0 | 52 / 58 |
| 2FVY | chain A, CA311/CA | 7 | 6 | 0 | 52 / 58 |

Both declared 309-residue sequences exactly match 1GLG. Mapping uses deposited
polymer sequence positions and source connectivity, not residue-number
arithmetic. Both lack only terminal positions 1 and 307–309. The unmodified
3.1 Å qualification and 3.3 Å inclusion policies naturally recover positions
134, 136, 138, 140, 142 and 205. Asp138 contributes OD2 in both new structures,
whereas 1GLG uses OD1. This changes no residue membership or donor type.

Both local representations and their cap-neighbor atoms are complete and have
no alternate conformations. The established residue-level occupancy policy
selects remote alternate conformations before PDBFixer reads the monomer.
Original alternate identities, occupancies and coordinates remain in each
`source_selection.json`; the single selected conformer is then written with
blank alternate labels. No structure or conformation was selected from scores.

The nearest deposited water is 4.466492 Å from Ca in 2FW0 and 4.512841 Å in
2FVY. No nonstandard species intersects either frozen inclusion shell. The
2FW0 citrate/sodium/malonate and 2FVY glucose/acetate/CO2/glycerol remain in the
source and diagnostic inventories. The reported Glu149 radiation damage is
retained as a source caveat; its nearest selected atom is 26.884379 Å from Ca
and 22.937135 Å from the retained local heavy atoms. It is outside both cores.

One pH-7 protonation preparation per source is shared by v3 and extended
representations. PDBFixer 1.9.0/OpenMM 8.1.2/gemmi 0.7.5 ran in the existing
`/home/jwestrob/miniconda3/envs/fep` environment, with one CPU thread and GPUs
disabled for this preparation process. No missing residue or missing ordinary
heavy atom was rebuilt. The standard protonation wrapper added only terminal
OXT at residue306, outside both cores. Every original selected source heavy
atom and every retained QM source heavy atom has **0.0 Å displacement**.

All four La/Ca pairs have identical coordinates and ligand ordering, singlet
multiplicity, charges 0/−1, and the unchanged native r2SCAN-3c/CPCM(Water)/DefGrid3
recipe. No bands or calibrated zero are assigned. The eight inputs are ready;
their execution and results belong to the parent runner's records.

## Artifacts and accounting

- Source CIFs and HTTP/hash receipts:
  `workspaces/ggr_mechanism_20260915/stage_b_sources/`.
- Authoritative source audit and preserved implementations:
  `workspaces/ggr_mechanism_20260915/stage_b_audit_v2/source_audit.json`.
- Four preparation manifests/eight endpoint inputs, protonation provenance,
  source maps and timings:
  `workspaces/ggr_mechanism_20260915/stage_b_prepared_v1/inventory.json`.
- First audit v1 is preserved. Its diagnostic inventory mistakenly listed
  standard HIS as nonstandard because the topology uses HID/HIE names. v2 fixes
  that reporting alias only; the local gates, selected atoms and preparation
  settings did not change. No energies were inspected during either audit.

Preparation across both sources used 11.641096 seconds wall time and 11.247494
process CPU-seconds. Process-wide peak RSS was 213,520 KiB. Hydrogen placement
alone accounted for 9.197413 wall-seconds and 9.115599 process CPU-seconds.
These are local preparation measurements, not allocated cluster cost. Source
download CPU time was not instrumented. No GPU or ORCA calculation ran here.

Six tests on real pinned artifacts pass: natural composition recovery, retained
crystal/radiation caveats, sequence mapping after author-ID renumbering, explicit
failure of a corrupted copy missing the actual amide N, frozen shared protonation
and paired invariants, and source-heavy-coordinate identity. No synthetic
scientific result or executable substitute was used.

## Executed commands

From the repository root:

```bash
python scripts/ggr_structure_prepare.py audit \
  --source-2fw0 workspaces/ggr_mechanism_20260915/stage_b_sources/2FW0.cif \
  --source-2fvy workspaces/ggr_mechanism_20260915/stage_b_sources/2FVY.cif \
  --reference diagnostics/laca_benchmark_expansion_20260915/sources/pdb/1GLG.cif \
  --reference-repair workspaces/affordable_challenger_20260915/verified_repairs/ggr_1glg_GGR/repair_manifest.json \
  --agreement diagnostics/ggr_mechanism_plan_20260915/AGREEMENT.md \
  --output workspaces/ggr_mechanism_20260915/stage_b_audit_v2

OPENMM_CPU_THREADS=1 CUDA_VISIBLE_DEVICES='' \
  /home/jwestrob/miniconda3/envs/fep/bin/python scripts/ggr_structure_prepare.py prepare \
  --source-audit workspaces/ggr_mechanism_20260915/stage_b_audit_v2/source_audit.json \
  --stage-a-ready workspaces/ggr_mechanism_20260915/stage_a_prepared_v1/inventory.json \
  --output workspaces/ggr_mechanism_20260915/stage_b_prepared_v1

python -m unittest discover -s tests -p test_ggr_structure_prepare.py -v
```

The first two commands intentionally refuse to overwrite their completed
outputs. Continue by consuming the pinned preparation inventory; do not repeat
hydrogen placement to obtain another source state.

## Interpretation

These remain development conformations of **one GGR observation**. The pair
differs in crystal conditions, cleft occupants and reported radiation damage;
it tests structural robustness, not a causal sugar-only effect. The retained
sugar is outside both cores. Input integrity does not establish predictive
success, and this preparation changes neither the default scorer nor 1GLG's
archived records.
