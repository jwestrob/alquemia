# Mex static whole-system additive model

Declared16October2026 before new molecular calls. Root owns submission.

## Finite question and Hamiltonian

Can the exact repaired Mex8FNS EF3 origin be assembled into one explicitly
defined energy and analytic physical gradient for both target metals?

`E = MACE_vac(Q) + [GFN2_emb(Q;MM) - GFN2_vac(Q)] + MM_ledger + G_ALPB(qFF,R)`.

Q is the existing209-atom capped local region, charge0; full R has2049 physical
atoms, including all retained waters, the target and three Nd spectators. Use
the exact12569 preparation and12570 four accepted GFN2 gradient cells. Preserve
all source coordinates, caps, charge shifts, field identities and occupancy.

Run exactly two new native MACE vacuum energy/force calls on the existing209
coordinates: LaIII physical multiplicity1 and DyIII physical multiplicity6.
Use the same pinned native checkpoint, float64, head, software and state batching
as12500; no external field is supplied to MACE. GFN2 retains its separately
declared original f-in-core effective singlet. No whole-protein MACE spin or SCC
calculation occurs. This combination is a declared approximate model, not a
qualified equivalence between the physical and effective electronic states.

Reuse `scaffold_origin_transfer_v1/Mex8FNS/{La_A,Dy_A}.json`: retained bonded,
shifted-charge MM Coulomb, MM LJ, and realQM–MM LJ, with their actual gradients.
These bind the same repaired INPUTS. Reorder strictly by source ID, including
the target alias, never by assuming metal append order. Preserve the existing
fixed-cap boundary policy: radial cross bonds retained, angular terms represented
by caps omitted, unrepresented exterior terms retained. No extra QM–MM Coulomb,
QM–QM classical terms, synthetic-cap FF particles or C4 induction. The embedded
minus isolated correction already supplies field interaction and response once.

## One fixed physical solvent charge model

Use the original unmodified ff19SB/TIP3P charges from all2045 records in
`lanm_ef3_preparation_v1/Mex8FNS/export_v1/atoms_charges.json`, joined by source
ID to the repaired2049-atom physical map. Append the already declared formal+3
charge at the EF3 target and each of three Nd spectators. Total charge is+2.
Do not use the ledger's zeroed QM charge entries or shifted MM field as a global
charge vector. Do not project synthetic-cap QM populations. No synthetic caps
enter this cavity. Hold these declared FF charges fixed; serialize with17digits
and require exact float round trip, without renormalization or fitted offsets.

This is explicitly a force-field-charge solvent approximation distinct from
the PQQ and LBT GFN-population solvent models. It supplies no charge-response
derivative. A later LBT exchange control must use the same stated charge rule;
no old LBT solvent scores, calibrated offsets, thresholds or affinity labels
transfer into this result.

This Mex result concerns EF3 at fixed three-Nd spectator occupancy, not a whole
protein affinity label. A later Mex/LBT preference comparison needs a declared
same-rule local-QM/physical-MM partition, cap/boundary policy, MACE plus local
embedded-minus-isolated GFN2 correction and FF-charge global solvent on LBT.
The previous whole-LBT MACE plus self-consistent-solvent energies cannot be
combined with this Mex partition by transferring offsets. No such additional
LBT preparation or campaign is included here.

Run exactly two native full2049-atom ALPB energy/gradient calls. Use the same
native fixed-charge functional, built static library and water settings whose
26-call interface test12566 passed. Compile a separate thin wrapper changing
only the declared atom-count restriction and charge-rounding bookkeeping;
no engine or parameters change. Native H/C/N/O/La/Nd/Dy tables are present.
Lanthanide generic nonpolar/HB entries and descreening defaults are not
independent validation of Nd or trivalent-ion transfer. The actual physical
metal element determines its cavity radii while its declared charge remains+3.

## Admission, assembly and limits

Before calls validate exact parent pins and successful states, source IDs,
coordinates, all four caps,2049 physical sites,1840 field sites, charge closure,
fixed-vector La/Dy identity and existing MM ledger/components. Local MACE uses
the existing native model/state checks; require finite energies and209×3 forces.
ALPB requires native successful return, finite components and2049×3 gradients,
and charge closure≤1e−5e. Preserve native energy units and convert exactly once.

Map `gMACE = -F * eV_to_kcal` through both existing cap-anchor Jacobians. Add
the saved mapped `gGFN_emb - gGFN_vac`, reordered classical gradient, and global
fixed-charge solvent gradient. Report each component, total, target Dy−La
contrast and full forces. Record translation residuals without inventing new
tolerances or claiming a finite-difference pass. Reconstruct the saved GFN
component difference from its actual paired rows before reuse. Require matching
saved classical total and component sums; do not silently reuse wrong geometry.

This static origin result is a model assembly diagnostic. No optimization,
motion, new geometry, Hessian, affinity inference, source ranking or classifier.
Existing coarse classical derivative failures remain preserved. A separately
declared total mixed/boundary directional finite-difference check is mandatory
before response optimization; it must include cap-anchor and external-site
motion, rather than testing only a local metal translation.

Existing source graph/bonds/XML and coordinate/cap/field-update primitives permit
later motions without new protonation, topology or force-field preparation.
The current Hans eight-mode adapter is source-specific; no Mex motion basis has
yet been prepared or qualified. Reuse primitives, not the Hans motion manifest.

## Execution

One GPU allocation for two warm MACE calls, with the two small native solvent
calls and wrapper compilation in the same allocation. Native solvent uses one
thread. No Slurm memory/time override. Record actual calls, timings, allocation,
runtime/model/build pins and all failures. New GFN2, MM, DFT calls are zero;
previous4GFN and2MM cells are counted as reuse. Root arms completion wake and
owns submission. No numerical work on login or production/default changes.
