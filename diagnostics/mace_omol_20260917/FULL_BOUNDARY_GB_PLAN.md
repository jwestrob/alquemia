# Full-boundary QM-charge solvent extension of the matched vacuum hybrid

Declared after NORMALIZED_CHARGE_REPORT.md, before assembly or solver energies.
Active-goal authorization applies. This is a new opt-in development candidate;
the failed vacuum hybrid and production baseline remain immutable.

Protocol: `normalized_QM_projected_ff19SB_full_OBC2_vacuum_hybrid_v1`.

## Scientific question and fixed systems

Does an inexpensive reaction field on the actual full-protein boundary improve
the matched vacuum hybrid without reintroducing severe partition dependence?
Use the four consumed normalized preparations and actual charges in
normalized_charge_v1/manifest.json (SHA31caf15249cbfab1ff381f1185c7243e77c69a003c735f52ca80e5ad63bc3fc6),
normalized_charge_report_v1/result.json (SHA2689bcdfa9b3f2d0358c5e21fa8b5b38533b36502c7e5c3f9826ee7487041d10)
and matched_H_report_v1/result.json. Require successful charge/projection gates
and verify actual receipts. Reuse every quantum/MACE endpoint; no new DFT,
MACE inference, geometry search, training or mechanical correction.

GGR extended/connected use identical full physical atoms. Alpha1F6S and6IP9
retain their distinct original inventories of two/three waters. Preserve source
coordinates, hydrogen rule, assembly, donor identity, protonation and disulfides.
No water, sugar, cofactor, residue or metal is added/omitted to rescue a result.

## Charge model: one local boundary rule

Use the actual projected QM charge distribution Pq_M from the completed charge
trial. Map it to the full physical atom array. P includes real source atoms,
selected metal, complete retained waters and the omitted real cap anchors.
No synthetic cap sphere remains. The entire projection support is excluded
from forcefield charges, including omitted anchors receiving cap charge.

Prepare ff19SB charges for exactly the already declared standard protein chain,
using the same pinned XML/source protonation/disulfides. Replay physical IDs,
connectivity and heavy coordinates against the existing preparation. Never
assign forcefield defaults to La or a cofactor. Retained waters are entirely
quantum; no forcefield charge is added to them. Unsupported topology fails.

Use the existing preparation's formal-charge ledger to assign each retained
charged sidechain to its source residue; neutral backbone/water entries add0.
Require the sum to equal the recorded ligand formal charge. A nonzero ledger
entry whose complete charged functional group is absent from the QM source
atoms is unsupported. Merge all contributions for the same source residue.

For each protein residue r touching projection support, let Q_r be its original
integer ff19SB residue charge, L_r its recorded selected-fragment formal charge,
and F_r the sum of untouched forcefield charges on atoms outside support.
The target exterior charge is Q_r-L_r. Redistribute
  delta_r = (Q_r-L_r)-F_r
uniformly only over atoms of that residue outside support that are directly
bonded to a support atom (including H). Retain exact recipients, bonds,
original charges and increments. If delta_r is nonzero and no such atom exists,
fail explicitly. Fully selected neutral/charged residues must close without
an external recipient. No residue-number adjacency, global neutralization,
endpoint-dependent redistribution or guessed charge for missing atoms.

The resulting q_env is identical for Ca/La, exactly0 throughout projection
support, and has the required exterior formal charge. The total is
q_M=Pq_M+q_env. Its small native printed-QM charge error is retained, not
renormalized; require full closure within5e-5e and exact paired environment.
This is a documented approximate boundary scheme, not a uniqueness claim.
Its energy sensitivity is assessed by the GGR partition test below.

## Energy expression and solver

H_M = E_DFT,vacuum(core,M)+T_mask(full,M)-T_mask(core,M)
A_M = H_M + G_full(q_M)
R_A = A_Ca-A_La.

G is ONLY the OBC-II electrostatic reaction-field charging energy, including
self and cross terms, from maintained OpenMM8.5.1. Reuse the existing isolated
CUDA12 solver environment and numerical implementation; do not install over
shared environments. Use all physical atoms with no cutoff or periodic box,
solute dielectric1, solvent78.5, salt0, nominal298.15K, surface-area term0.
Use the already recorded Bondi CHNOS radii and common Ca/La1.8A radius and
existing descreen scales(H.85,C.72,N.79,O.85,S.96,metal.8). Physical sphere
inventory is identical across endpoints and GGR partitions; zero-charge atoms
still contribute to the boundary. Metal cavity parameters remain unvalidated
and frozen; do not tune them after inspecting signs.

There is no core CPCM term and no bare core/environment Coulomb addition.
Learned vacuum interactions already enter T_mask(full-core). This extension
does not repair missing long-range vacuum response in that model. Forcefield
charges represent a fixed, already effectively polarized protein background;
no explicit induced-dipole energy is added. QM charges remain vacuum-derived,
not self-consistent with protein/solvent. These mismatches remain limitations.

For component audit, on the SAME full physical boundary retain G_total,
G_QM=G(Pq_M), G_env=G(q_env), and
G_cross=G_total-G_QM-G_env. G_env cancels in Ca/La because q_env/boundary are
identical; G_QM and G_cross generally do not. There is no isolated-core cavity
or subtraction. Raw contrasts have no compatible aquo reference, calibrated
zero or inherited baseline band. Aqueous binding affinity/class remains null.
Any forces returned by the solver hold charges fixed and are not gradients
of the combined model; relaxation remains disabled.

## Finite solver tasks and fixed gates

48 solver evaluations, no new high-level endpoints:
-8 primary full-charge states,8 quantum-only charge states,4 environment-only
 states (one per representation); every state retains the full physical cavity.
-8 identity states with solvent dielectric=solute dielectric=1 and full charges.
-8 jointly rotated and8 translated full-charge states using the existing
 rotation matrix and translation[10,-7,3]A, independent of labels.
-2 repeated GGR-connected full endpoints;2 corresponding CPU Reference solves
 compared to CUDA-double, on the same inputs. No alternative scientific model.

All values finite; state/projection/charge ledgers exact to declared roundoff.
Identity absG<=1e-6kcal; rigid/repeat/Reference-CUDA energies and contrasts
<=0.01kcal. Component and direct algebra closure<=1e-7kcal. No grid refinement
applies to this no-grid analytic approximation; cross-platform agreement tests
numerics, not the accuracy of OBC against PB or experiment.

Retain the pre-existing GGR connected-minus-extended total-score criterion
abs(deltaR)<=2kcal-scale and all four alpha-minus-GGR differences>0.02kcal-scale.
Both representations and both alpha preparations stay in all denominators.
No changing charge scheme, radius, dielectric, water inventory, threshold or
preferred structure after outputs. A failure stays a failure. These remain two
consumed biological groups with qualified cross-study directions; no broad or
blind accuracy claim even if all pass. Do not calibrate on these four results.

Use existing manifested runner, locks, cache/receipt machinery; add minimal
schema/worker dispatch if needed. One A5000/16CPU/64474MiB host allocation;
existing GB warm evaluation~0.09s/full system, process startup dominates.
Expect minutes, not hours; no project CPU/time stopping budget. Count all
attempts and resource overhead. Prepare and test source/boundary invariants
before submission; actual scientific tests remain explicitly unrun until
outputs exist. Baseline and unrelated active jobs remain untouched.

Sources: [OpenMM OBC expression](https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html#gbsaobcforce),
installed `openmm/app/internal/customgbforces.py`, and the pinned earlier
mace_gb implementation/report. These specify component equations; none
validates this new assembled La/Ca discriminator.
