# ExaA source and packing tests — 14 October 2026

Jacob approved the two proposed tests with “yep. go”. These are contained
research comparisons; current production, calibration and historical results
remain intact. Root owns molecular submissions and scoring. Crystal preparation
and fixed-pool packing preparation have separate agent ownership.

## Questions and fixed population

1. Can the current frozen-charge candidate score both chains of the existing
   same-protein ExaA crystal 1FLG, after checking chemical-state compatibility?
   Compare with the saved ExaA AF3 sources. A crystal/prediction difference is
   source sensitivity, not proof that the prediction is wrong. Source redox,
   adduct, retained waters, occupancy, construct and alternate conformers must
   be explicit. Unsupported chemistry stays unavailable. The existing dry,
   oxidized-PQQ3− hypothesis is not an experimental redox determination.
2. Does explicitly represented surrounding packing change the differential
   work of the existing donor movements? Fix ExaA/Q9Z4J7 two-Ca samples **0 and
   3**, the first same-conditioning correct control and the ambiguous target.
   Use the exact saved origin, Ca proposal and La proposal; both metals see the
   same pool. No optimization or folding in this branch.

Packing membership retains the current polar context and adds complete mapped
fragments containing any heavy atom within **4.3 Å** of genuinely moved donor
heavy atoms in any saved candidate. Take the union of fragment identities across
both declared sources and all three candidates before evaluating energies. This
extends the old polar-atom envelope to swept all-heavy contacts. It is not a
Leu-only correction. Source heavy atoms, original donor motions, protons, PQQ,
water inventory and spectator treatment remain traceable and fixed; new caps
follow the source graph. Unsupported fragments fail explicitly.

## Energy and interpretation

Use the already qualified candidate expression

`E_M(q) = E_OMOL,vac,M(q) + G_nativeALPB(q_GFN2,vac,M(q); q)`.

Reuse the pinned OMOL weights, CPU float64 transport, published-parameter native
xTB6.7.1 and native ALPB replay. ORCA6.1.1 defaults are not an alternative:
any ORCA native-GFN2 invocation would require the full corrected parameter file
and existing guard. The planned cheap endpoint path uses the existing corrected
standalone implementation and its independent component/parity checks.

For each metal, report `W_M(q)=E_M(q)-E_M(origin)` separately in each context.
The paired response is `W_Ca-W_La`; report the large-minus-small change in that
response, MACE and solvent components, and row-minimum selections. Absolute
totals across differently sized contexts are not compared. Larger Ca−La
contrasts are La-like. Never apply the current bands to a new context and call
that a validated classification. New calibration, if warranted later, belongs
to the designated canonical population rather than these two cases.

All twelve packing matrix cells are required for the two complete pools. Reuse
the twelve archived small-context cells. Failed new cells stay unavailable.
Each successful energy must retain its actual coordinate/state/software pins,
receipt, unrounded components and existing numerical tolerances. Use the existing
0.01 kcal/mol CPU-MACE replay and 0.001 eV/Å force gates; retain the qualified
native component precision policy. Report continuous effects and limitations,
not a fitted success threshold. No new DFT is needed for these first tests.

If the packing effect largely cancels between metals, close this explanation
for these motions. An independently referenced, finite donor-motion check is
the contingent next question, not an automatic broad quantum campaign. Earlier
extra-Asp DFT tests do not establish the current Glu/coupled response.

## Execution and delivery

All numerical preparation, tests and molecular work use Slurm under the current
shared-node policy (at most12CPUs/job,60/user, no memory/time override). Reuse warm
MACE workers and independent one-thread solvent endpoints. Record actual work
and failures. No project-wide arbitrary compute cap or production rescore.
Arm explicit completion wakes before yielding with outstanding jobs.

Candidate products live under `workspaces/exaa_source_packing_20261014/`;
compact plans/code/results here. No existing results are overwritten. Deliver
both source/state audit and paired deformation tables, then recommend the next
scientific action from those results. No automatic production promotion.
