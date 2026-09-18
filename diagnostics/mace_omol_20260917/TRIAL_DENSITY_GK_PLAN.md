# Saved responsive quantum density in the same AMOEBA-GK/MACE functional

Declared after frozen-density hybrid1201074 failed partition and ordering,
before any new density queries/native scores for this protocol. Jacob's
autonomous MACE goal authorizes this contained experiment. Production baseline
and every earlier record remain unchanged. All cases are consumed development.

## Question and inputs

Does replacing the vacuum quantum density by an already available protein-field
responsive trial density improve the same full physical electrostatic hybrid?
This tests a missing quantum response component while retaining native protein
induction, permanent multipoles, cavity and learned context. It is not a new
fit, a dielectric/radius scan, or an assertion of self-consistency.

Use all eight matching real endpoints from
`workspaces/mace_omol_20260917/responsive_quantum_result_v2.json` and the actual
charge/potential outputs `responsive_charge_report_v2/result.json`. Those
wavefunctions responded to the previously declared permanent ff19SB field.
Its different source field is explicit provenance: it supplies a trial density,
not a purported self-consistent solution of the current AMOEBA/GK functional.
Use every preassigned responsive density; do not select vacuum versus responsive
per endpoint, choose a favorable SCF branch, or minimize over inspected scores.

All physical/QM coordinates, atoms/caps, formal states, mappings, assembly and
0/2/3waters must exactly match the frozen-density parent's corresponding case.
Same r2SCAN-3c method/basis/ECP/D4/gCP and electronic states. The saved embedded
quantum outputs already passed their original checks. Reverify actual inputs,
receipts and surviving wavefunctions; no recomputation is implied if missing.
Retain all prior successes/failures and previously inspected evidence status.

## Energy accounting: remove the old generating field exactly once

For the saved density rho*, reconstruct its intrinsic core energy from actual
native output and its actual old-field potential observations:

```
E_core[rho*] = E_DFT,old_embedded[rho*] - V_old_field[rho*]
V_old_field = sum(old_environment_charge * phi_rho*(old_charge_position))

A_M[rho*] = E_core,M[rho*] + V_current_AMOEBA[rho*]
           + G_current_GK[projected_charge(rho*)]
           + I_current[density_field(rho*), proxy_RF(rho*)] - I_environment
           + T_short(full,M) - T_short(core,M)
R = A_Ca - A_La
```

The original `%method DoEQ false` excludes external charge self-energy.
Native embedded DFT includes old direct electronic+nuclear interaction; leaving
it in would double count/retain the wrong environment. Audit the subtraction
against independently stored `intrinsic_core_kcal_mol` and raw native energies,
old point charges, exact potential arrays and coordinate order. Original
source atoms are only the core; the old external field is not part of the
bare-core `orca_vpot` output already qualified by that experiment.

Same-geometry intrinsic trial core energy should not be below the actual vacuum
SCF minimum by more than0.05kcal/mol. Freeze that diagnostic tolerance now,
matching the prior quantum variational allowance. No clamp or energy shift.
Because the current functional includes a fitted GK proxy and an approximate
MACE split, do not claim a rigorous variational bound for the complete hybrid.
The old-field wavefunction is not stationary in the new field. No combined
gradient, force, relaxation, entropy, absolute reference or decision band.

Use the exact same current AMOEBA residue charge ledger, source-frozen mask,
common Ca2018 metal cavity, native solver/library, dielectrics, exclusions and
MACE short receipts as1201074. Reuse actual environment-only static/response
outputs only after exact physical-state and dependency compatibility checks.
Replace all quantum-density terms together: intrinsic core energy, direct
potential/field/Hessian, GK source charges, reaction field and induced response.
Neither the old total GB score nor a frozen-density correction can satisfy a
responsive-density task. No missing correction becomes zero.

## Declared calculations and numerical checks

1. Read-only source audit:8embedded outputs,8saved-density/charge records,
   old-field subtraction, exact geometry/microstate/source support, compatible
   MACE receipts. No new DFT, charge fit, MACE, solve or geometry optimization.
2. Eight native `orca_vpot` queries from those saved densities. At every actual
   exterior site use the existing37point potential/Hessian stencil
   (center plus18offsets at0.01/0.005bohr) and the existing12point electric-field
   stencil (six offsets at0.001/0.0005bohr):49observations/site in one query.
   Check center potentials against the old responsive potential archive.
   Fine electric fields/Hessians supply the candidate; coarse values test
   numerics, never select favorable scores. No nuclear Hessians/DFT gradients.
3. Sixteen native initialization-only source/reference preparations, preserving
   the current boundary rules and using the existing responsive projected
   CHELPG charges. No new charge fitting or boundary redistribution scheme.
4. The same four variants and native checks as DENSITY_GK_HYBRID_PLAN:
   primary/rigid/common-metal-radius0.95/1.05, source states environment/Ca/La,
   primary response1e-7and1e-9Debye, plus three vacuum identity static states.
   Logical inventory51static energies,48field queries,60responses. Exact
   unchanged environment-only outputs may be reused with explicit receipts;
   report newly executed counts separately. No new high-level/model calls.

Keep existing center tolerance1e-8au, maximum field refinement1e-6au, paired
bare diagonal response refinement0.01kcal (diagnostic, not an environment
energy), endpoint/paired quadrupole refinement0.02/0.01kcal and rigid moment
contraction1e-8kcal. Keep all parent full-model convergence, identity, algebra,
radius-sensitivity, large-dipole, partition abs<=2 and four ordering>0.02 rules.
No relaxed criteria if the density change is inconvenient. Report how every
component changes, not only the final sign. A pass warrants broader evaluation,
not promotion or broad validation; a failure stays a failed challenger.

Native utility allocation8CPUs/16GB, one worker per endpoint; native response
allocation64CPUs/64GB with existing exclusions. Previous37point utility pilot
needed846wall seconds/6768allocatedcore-s and2951.569CPU-s.49point work should
be of the same order; measure rather than claiming this estimate as a result.
Current native full numerical pilot needed570wall/36480allocatedcore-s,
3659.722079CPU-s. No project budget or explicit wall-time cap. Retry technical
failures only with matching science and visible receipts. Do not launch new
DFT to replace a missing archive without declaring that separate computation.

New protocol:
`saved_responsive_trial_density_AMOEBA2018_GK_proxy_POLAR_short_hybrid_v1`.
Products under `workspaces/mace_omol_20260917/trial_density_gk_*`; reuse existing
utility/native execution machinery and pin all implementations. Original
vacuum protocol/fixtures must continue replaying unchanged. No new training,
long trajectory, global quantum solve or unsupported mechanical correction.
