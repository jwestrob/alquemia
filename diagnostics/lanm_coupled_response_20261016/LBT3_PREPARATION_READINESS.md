# LBT3 whole-peptide preparation readiness

16 October 2026. Source/code inspection only; no geometry preparation or scoring
was executed by this agent. Root owns the frozen pilot and Slurm execution.

## Verdict

Both consumed model-1 sources support one narrowly specified preparation using
existing native helpers. **The deposited whole peptides are chemically
incomplete at both termini and cannot simply inherit the old carve's charge.**
The adapter is ready for a Slurm preparation check, not yet numerically admitted.
LBT3 is the opposing-preference peptide control for the monomer LanM objective.
Its success alone would not validate LanM specificity.

## Exact deposited material

| Source | Deposited target | Model-1 peptide | Other deposited material |
|---|---|---|---|
| 7CCO | La, A101, (1.219, −3.957, −2.119) Å | Chain A1–17, 252 ATOM records | No waters, spectators, buffer or counterions |
| 7CCN | Lu, A101, (1.503, 3.134, 1.764) Å | Chain A1–17, 252 ATOM records | No waters, spectators, buffer or counterions |

Both sequences are `FIDTNNDGWIEGDELLA`: Phe1 Ile2 Asp3 Thr4 Asn5 Asn6
Asp7 Gly8 Trp9 Ile10 Glu11 Gly12 Asp13 Glu14 Leu15 Leu16 Ala17.
Both are synthetic, author-assigned monomers, 20-model solution-NMR ensembles
with representative conformer1. There is no deposited Dy structure here.
NMR conditions are pH6 and283K. 7CCO used400µM peptide/800µM La; 7CCN
used500µM peptide/1000µM Lu. Dissolved solvent is not an explicit atom inventory.

Raw files are `lanthanide_series_capability/lbt3/raw/{7CCO,7CCN}.pdb`;
the corresponding `first_models/*_model1.pdb` retain their model1 coordinates.
The existing `benchmark_manifest.json` pins both. Raw SHA256 values:

- 7CCO: `0977b2592106845b73650ad142e0a781f0d5fd226d87801e7d616a2505785e72`.
- 7CCN: `4bea8d556854f8398d8749530a23ed86fd587893d9c1299e3109de5a4650c575`.

Each deposit has Phe1/H1 only, without H2/H3, and Ala17/C,O without OXT.
The incomplete atom list does not define a valid isolated neutral terminus.
No acetyl or amide capping residue is deposited. The locally archived primary
paper `raw/PMC7656248_fulltext.xml`, paragraph26, explicitly uses negatively
charged Asp/Glu and a net−2 metal complex. This is compatible with the proposed
standard free termini but does not establish the experimental terminal chemistry.
Root's primary-evidence review additionally recovered supplementary FigS5's
peptide charge−4 statement, an unresolved discrepancy retained in `LBT3_PLAN.md`.
The declared scout assumes NH3+ and COO− cancel, Asp3/7/13 plus Glu11/14
contribute−5, and target LnIII contributes+3. Whole peptide charge−5; whole
complex−2. Native ff19SB must reproduce−5. This is a conditional chemical
state, not a verified microscopic protonation assignment.

## Minimal recipe and exact scope

`LBT3_PREPARE.py` preserves every deposited H and heavy coordinate, removes
metal coordination LINK/CONECT edges from the peptide mechanics by deleting
only the target before native parameterization, and adds exactly three atoms:

1. Ala17/OXT from existing `mace_global_prepare.terminal_position`: reflect the
   C→O direction about C→CA, length1.25Å. This is a terminal completion, not a cap.
2. Phe1/H2,H3 through existing OpenMM `Modeller.addHydrogens` at pH6,
   ff19SB, Reference platform and the established fixed seed20260918; explicitly
   select deprotonated ASP/GLU. OpenMM canonically names deposited H1 as H.
3. Keep every existing atom fixed during this native missing-H placement.
   The API internally minimizes only the two added H, at most50iterations;
   its return value exposes neither the actual iteration nor internal force-call
   count. Count the attempted/completed API calls separately from MACE scoring.

The resulting expected inventory is255 peptide atoms plus target =256.
No radial normalization of deposited H, interior reconstruction, water,
counterion, arbitrary atomic charge or whole-peptide fragment cap is added.
H-bond lengths and native equilibrium lengths are recorded as diagnostics.
Failures stay unavailable; no parameter, geometry or protonation retry occurs.

Full native peptide bonds and `peptide_system.xml` are exported with ff19SB
charges, unconstrained bonded terms and nonperiodic NoCutoff nonbonded terms.
Their particle mapping remains explicit. This is complete peptide mechanics;
it does not invent metal force-field parameters or calculate a holo MM energy.
The target has no artificial covalent graph edges. No pre-existing full native
LBT3 Amber preparation was located; prior material is the48-atom carve or
separate consumed Protenix folds.

## Matched fragment control and eight future cells

For each source the adapter calls the original `lbt3/build_benchmark.py`
`build_cluster` on atoms mapped directly from the newly completed peptide.
The exact legacy recipe retains D3/N5/D7/E11/E14 sidechains and W9/C,O,
adds the five CB→removed-CA H caps and the original two W9 carbonyl caps,
all at1.09Å. It emits a separate48-atom matched carve with ligand−4/holo−1.
Every shared real atom must equal its whole-peptide coordinate exactly.
The source/cap map distinguishes all seven generated H caps from physical atoms.
Old cluster files are preserved and are not overwritten or used as new score inputs.

Both representations substitute only the target for paired LaIII singlet
(multiplicity1) and physical DyIII high-spin sextet (multiplicity6), preserving
identical coordinates within each pair. No spin-orbit qualification is implied.
The pending pilot is2sources ×2representations ×2metals =8 fixed MACE cells.
No quantum input, GFN calculation, optimization trajectory or series campaign
is prepared. Flat task metadata includes exact XYZ pins, charges,
multiplicities, atom counts and explicit metal indices (whole255, carve0).

The whole-vs-carve comparison changes represented scaffold, terminal/D13 charge
and capping together. It is a controlled omission comparison, not a unique
energetic decomposition into “mechanics.” Both dry fixed structures omit bulk
solvent; the paper's proposed water-mediated mechanism remains unrepresented.
Raw Asn6 Ramachandran warnings and 7CCN model1 Glu11/Glu14 carboxyl angle
warnings remain; neither heavy-atom repair nor equilibrium claims follow.

The failed old `LBT3_CAPABILITY_VERDICT.md` protocol remains closed. Its
unavailable series curves are not evidence for this newly declared model.

## Runtime contract

Use the existing `lanm_qmmm` Python with NumPy/OpenMM on Slurm. No new package,
network call or scheduler wrapper is needed. Entry point:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python ABS_LBT3_PREPARE.py \
  --repo-root ABS_REPO --config ABS_CONFIG_JSON --output FRESH_ABS_OUTPUT
```

Config requires pins for `plan` (root's `LBT3_PLAN.md`), `benchmark_manifest`,
`forcefield` (existing `amber19/protein.ff19SB.xml`), and
`sources.{7CCO,7CCN}.{raw,model1}`. The adapter refuses execution without
`SLURM_JOB_ID`. It records dependency pins, raw-to-native mappings, the native
XML, bonds, terminal additions, preserved-source checks, elemental inventories,
H diagnostics,8 endpoint tasks and availability in a fresh directory.
Root should pin the existing legacy builder and imported native helper files
alongside the adapter before execution. Native completion failure is retained;
there is no automatic fallback to rebuilding missing nonterminal heavy atoms.
