# Proposed: whole-protein MACE with charge conservation by chemical group

**Status: declared pilot under Jacob's standing discretionary MACE-goal
approval; not yet executed.** The intervening generic AGENTS text does not revoke
his explicit instruction to pursue contained pilots without repeated approval.
The scope question sent during recovery is optional steering, not a new gate.
The previous coupled-response and collective-substitution studies are complete.

## Question and rationale

Can restricting charge redistribution improve the intact-protein discriminator?
The existing POLAR traces show a broad response to a one-unit metal-charge
change; they do not establish that the response is erroneous. The OMOL charge
mask avoids one global-conditioning problem but removes charge-state information.
Neither finding has yet been tested by the group-constrained POLAR model below.

This changes the learned density construction itself. It is not another additive
energy correction or a longer search along the failed donor-response path.
Keep all atoms and long-range field interactions in the intact protein.

Already tried: global POLAR with solvent, charge-masked OMOL, explicit QM fields,
AMOEBA/GK polarization, conductor corrections, and local mechanical response.
Do not describe those as new. Alpha 6IP9 is already a La-bound crystal structure;
using a La-shaped starting structure alone is not an untested remedy.

## One model, fixed before scoring

Use the existing medium POLAR checkpoint, analytic electrostatic kernel and
isolated environment. Medium is chosen for cost, not its benchmark scores.
Retain source geometry, protonation, assembly, PQQ state and water inventory.

At each of the three charge-restoration stages, replace the protein-wide sum
in the Fukui normalization by a sum over a fixed chemical group. For each spin
channel and atom i in group g:

    p_i = p_tilde_i + f_i / sum_g(f) * (Q_g / 2 - sum_g(p_tilde))

The paired endpoints are singlets. Preserve the learned Fukui weights, higher
multipoles, field updates and every energy readout. The endpoint's complete
energy must be recomputed with the changed densities; do not append a Coulomb
scalar to an old total. A single group containing the whole system must recover
the native model. Near-cancelling group weights fail explicitly (absolute sum /
sum absolute weights <= 1e-10); no denominator clipping or fitted damping.

Primary grouping: one joint metal/ligand block containing the selected metal,
PQQ where present, all retained site waters, and every complete protein residue
containing a source-mapped heavy atom of the existing extended-amide core.
Remaining residues are individual groups; join source disulfide partners.
No caps or new atoms are added. Map actual amide N atoms by the existing graph,
not residue-number arithmetic. The PQQ cases use their existing canonical
source mapping. Record every group before evaluating a model.

Use formal group charges from the existing protonation/template ledger,
including terminal charges; PQQ remains -3, waters neutral, metal +2/+3.
Do not use force-field atomic charges in the MACE energy. Missing or ambiguous
template/source mappings are unsupported. Verify that group charges sum to the
existing endpoint charge and groups partition every physical atom exactly once.

**Key limitation:** formal residue charges are a modeling constraint, not exact
integrals of the learned density. Peptide charge transfer across a group boundary
is restricted, and the pretrained model was not trained with these constraints.
The authors also caution against interpreting individual learned monopoles as
ordinary atomic partial charges. This is an empirical candidate to test, not a
verified physical repair or a constrained-DFT calculation.

Keep the existing frozen-monopole OBC-II solvent prescription unchanged:
dielectrics 1/78.5, salt zero, common Ca/La radius 1.8 A, original element radii
and descreen factors, surface-area coefficient zero. Use new endpoint monopoles.
Its unvalidated metal cavity and Gaussian-density/point-solvent mismatch remain.

    E_candidate(M) = E_group_POLAR(M) + G_OBC2(new_monopoles_M)
    R = E_candidate(Ca) - E_candidate(La)

No CPCM, aquo offset, absolute decision, entropy or relaxation correction.
Report vacuum and solvent components, but judge the predeclared total.

## Exact initial scope

Seven consumed structures: MxaF 1H4I, XoxF 4MAE, GGR 1GLG/2FW0/2FVY and bovine
alpha-lactalbumin 1F6S/6IP9. Existing physical preparations are inventoried in
SOURCE_INVENTORY.json. Keep both alpha structures and all three GGR structures.
This is four biological groups across separate evidence strata, not seven
independent experimental observations or prospective validation.

| Purpose | New MACE calls |
|---|---:|
| Seven structures, Ca and La, primary grouping | 14 |
| One-group identity replay, GGR 1GLG Ca and La | 2 |
| Existing rigid rotation, grouped GGR 1GLG Ca and La | 2 |
| Each GGR with the entire recorded connected-core residue set merged into the metal block, Ca and La | 6 |
| Total | 24 |

Run 24 matched OBC-II evaluations on the returned densities. No new DFT,
training, trajectories, structural optimization or label fitting. Preserve
failure receipts; do not expand the task inventory automatically.

Numerical checks: existing 0.01 model-kcal energy/replay/rotation tolerance,
1e-5 e group and total charge closure, exact paired atom/coordinate invariants.
The grouping sensitivity is an end-to-end representation test: absolute
change in each GGR contrast <= 2 model kcal, using the existing research target.
It is not expected to be mathematically invariant.

Predictive screen, in parallel with qualification: XoxF minus MxaF and all six
alpha-minus-GGR contrasts must exceed 0.02 model kcal. Report all seven raw
directions even if qualification fails. A gain only on 1GLG is insufficient.
No new bands or threshold fitting. Full canonical calibration, paper-domain
evaluation and model-size sensitivity would be later, separately stated stages.

Expected scale: tens of GPU minutes on one existing A5000 allocation, with
16 CPUs and 64,474 MiB host RAM; this is extrapolated from previous medium
whole-protein receipts, not a measured cost of the proposed model. Retain the
existing runner, caches and allocation policy; no installed environment edits.

## Sources and implementation feasibility

The installed `mace/modules/extensions.py` has explicit initial and subsequent
restoration blocks. A task-local adapter can replace those two blocks while
retaining the existing electrostatic kernels and readouts. It needs a one-group
identity test before interpreting changed-model results; no adapter exists yet.
Source-code and input hashes are in the inventory. Source inspection launched
zero model, solvent or quantum calculations.

- [MACE-POLAR equations 13-14 and 26-27](https://arxiv.org/html/2602.19411v1).
- [Model documentation and density interpretation](https://mace-docs.readthedocs.io/en/latest/guide/polar_mace.html).
- [Actual earlier charge tracing](../mace_response_trace_20260916/REPORT.md).
- [Existing whole-protein solvent/model experiment](../mace_global_benchmark_20260916/PLAN.md).
- [Latest coupled-response failure](../mace_site_response_20260918/COUPLED_PATH_REPORT.md).
- [Latest collective-substitution failure](../mace_collective_20260918/REPORT.md).

Baseline/default, all previous protocols and results remain unchanged.
