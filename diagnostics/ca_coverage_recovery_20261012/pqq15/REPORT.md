# All fifteen failures are off-site Ca placements

Slurm job **2866** completed on `node-64-768g-5`. All fifteen exact saved
Ca-conditioned structures retain a complete, schema-valid PQQ molecule, but the
selected Ca is about 25 A away. **No metadata repair is justified.** The unchanged
preparation gate rejects these geometries correctly; all fifteen are recorded as
`requires_new_site_conditioned_structure` in the separate recovery ledger.

## Actual geometry evidence

All rows are the original seed-1 samples 0–4, with no source substitution.
Distances below are read from the completed Slurm audit.

| Protein | Enzyme | Ca to nearest PQQ atom across five samples, A | Residues around off-site Ca | Typed CN |
|---|---|---:|---|---:|
| C5AXV8 | ExaF | 25.016–25.209 | Asp36, Thr39, Asp42 | 6 in each |
| Q88JH5 | PedE | 25.087–25.162 | Asp53, Thr56, Asp59 | 6 in each |
| Q9Z4J7 | ExaA | 25.038–25.127 | Asp45, Thr48, Asp51 | 6 in each |

Each off-site site has backbone and side-chain O contacts from these three
residues. The six contacts satisfy the original chemically typed 3.1 A donor
definition, but none supplies the required PQQ/core site. The PQQ association
cutoff remains 4.0 A, with 3.3 A buffered QM donor inclusion. A count of six at
a different site does not admit the intended fixed core.

The exact PQQ selectors are `C:PQQ1`, using the existing
`protenix_generic_pqq_v1` schema, 24 source heavy atoms and explicit
`pqq_ox_3minus_v1` (C14H3N2O8, charge −3, singlet). The allowed PQQ donors remain
N2/O1/O5/O6. The raw Ca atom remains `B:LIG_B1:CA1`, residue-aliased to
`B:LA1:CA1`; its element, atom name and coordinates remain native Ca.

Declared PQQ-site roles are:

| Protein | Glu anchor | Asn anchor | Catalytic Asp | D+2 homolog | Cationic partner |
|---|---|---|---|---|---|
| C5AXV8 | Glu198 | Asn275 | Asp317 | Asp319 | Arg344 |
| Q88JH5 | Glu221 | Asn308 | Asp358 | Ser360 | Arg386 |
| Q9Z4J7 | Glu213 | Asn300 | Asp350 | Ser352 | Arg378 |

These intended site residues have their nearest O/N atoms about 27–33 A from the
selected Ca. The frozen 4.3 A context memberships and original state anchors are
retained in the manifest. Context/map compatibility was **not evaluated**, because
none of the fifteen sources passes core admission. No context was built around
the off-site ion or substituted for the intended PQQ site.

## Normalization and source preservation

For all fifteen sources, unchanged canonical normalization reproduces the saved
normalized PDB bytes exactly. The audit verifies zero displacement and identical
atom identities/occupancies for every raw heavy atom through normalization and
saved protonation; all copied protonated atoms are also unchanged. Raw heavy-atom
counts are 4542 C5AXV8, 4921 Q88JH5 and 4848 Q9Z4J7. Every repair-action list is
empty. The PQQ microstate helper instantiates its usual three ring hydrogens for
schema validation; no protein reprotonation or new prepared endpoint occurred.

The local enzyme-label source is the benchmark FASTA headers in
`/groups/banfield/projects/multienv/corkscrew/supplementary_structures/`:
`ln_binder_positive_controls.faa` (C5AXV8 ExaF),
`ca_binder_positive_controls.faa` (Q88JH5 PedE, Q9Z4J7 ExaA).

## Consequence and limits

These are genuine source-placement failures for the intended PQQ-site scoring
task. They do not establish an experimental second site, biological preference,
or a defect in PQQ normalization. Root owns the newly authorized, separately
declared structure-generation route. The consistent distant DxxTxxD-region
placement must be considered when defining its input state; simply translating
the ion or increasing the PQQ radius would change the scientific question.

Original 125-source outcomes and frozen canonical bands remain unchanged.
The metadata-repair and successful-core/context branches were not exercised,
because no input justified them. There were no energy, optimization, refold or
protein-protonation calls in this job.

## Evidence, execution and review

- Driver: `PREPARE.py`, SHA256
  `5a6027c294a4751d47658e706529cfe615f5f3bb26a27751d681d283efb0a4bf`.
- Frozen manifest: `workspaces/ca_coverage_recovery_20261012/pqq15/run_v1/manifest.json`,
  SHA256 `ff2b41b1db1294243870d0cd003c483578a9ae713836f2c75290aa8150b9a80d`.
- Actual output: that directory's `execution_v1/{AUDIT.json,PREPARATION.json,EXECUTION.json}`;
  per-source `cases/*/{AUDIT.json,ADMISSION.json}` retain all pins and exact contacts.
- PREPARATION SHA256:
  `95923f4da3e516e158b59cd93fa071d2233163348d8dc66218ed8642d23d492e`.
- Measured cost: **12.755131 s wall, 8.657201 CPU-s**, peak RSS 192424 KiB,
  one numerical thread. Root owns submission, scheduler accounting and wake.
- Login checks were AST syntax parsing and metadata/byte pinning only.
  All molecular parsing, normalization replay and geometry checks ran on Slurm.
- Independent bounded code review by `/root/expanded_launch_review` found no
  source/admission relaxation, hidden optimization or launch blocker. It was a
  read-only review, not an additional scientific execution.

The runnable entry is `RUN.sbatch ABSOLUTE_MANIFEST_JSON`; the completed v1 output
refuses overwrite. Do not resubmit this completed diagnostic.
