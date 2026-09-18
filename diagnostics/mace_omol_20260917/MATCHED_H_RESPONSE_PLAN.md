# Matched vacuum hybrid plus bounded metal response

Declared2026-09-18 under Jacob's active autonomous MACE goal, before any new
response output. This is a separate candidate, not an addition to production
or a reinterpretation of the completed CPCM response experiment.

## Question and systems

Does structural accommodation improve the already implemented, representation-
consistent normalized-H vacuum hybrid? Its static alpha-minus-GGR margins are
+4.329632/+4.211135 for1F6S and−0.805099/−0.923595 for6IP9 (extended/connected).
The separate CPCM/POLAR response gave real, DFT-confirmed margin improvements,
but its gradients/geometry differ and its numerical corrections cannot be reused.
Reuse methods and archived compatible artifacts, not incompatible corrections.

Use all four existing normalized-H preparations from matched_H_prepared_v1:
GGRextended, GGRconnected, ALPHA1F6S and ALPHA6IP9; bothCa andLa. Keep each
source/core/assembly/protonation/water/cap inventory and all coordinates except
selected metal fixed. These are consumed development, two biological groups.
No new reference, threshold, label, fitted spring or dielectric. Preserve the
static hybrid, original-H and all CPCM experiments as separate protocols.

## Energy and approximation

Let T be the unchanged masked MACE-OMOL descriptor (same100M checkpoint,
float64, charge-feature mask and spin policy). For each representation,

    H_M(u) = E_DFT,vac,core,M(u) + T_full,M(u) - T_core,M(u)
    g_H = g_DFT,vac,core + g_Tfull - g_Tcore
    K_cheap = Hessian_u(T_full)

Thus K_cheap approximates the DFT core curvature by the matching learned core
curvature; the latter cancels its subtractive counterpart. This is a testable
curvature approximation, not an assertion that T has physical forces. No CPCM,
GB, GK, conductor, polarization correction or duplicate solvation is included.
The native DFT method stays r2SCAN-3c/DefGrid3/TightSCF vacuum/analyticEnGrad,
exactly matching this hybrid's saved centers. The energy is a mixed descriptor,
not an aqueous binding free energy. No fluctuation entropy is added.

Use the existing three Cartesian metal coordinates. K comes from all37points
(center,6axis+12mixed corners at0.02A, and18noncentral at0.01A), preserving
mixed terms. Match center analytic MACE gradients to the same physical XYZ;
source-to-metal mapping must be explicit, caps have no independent motion.

For every state solve min(g_H.u+0.5u.K.u),||u||<=0.20A using the already tested
positive3x3 sphere solver. Reject nonpositive coarse/fine K; no eigenvalue
clamping or radius changes. Fine curvature is primary. Require refinement
change at predicted point<=0.02modelkcal and coarse/fine displacement<=0.01A.
Require identical typed donor membership before admitting native validation.

## Finite execution and reuse

Reuse eight completed normalized vacuum native DFT energy/gradient centers,
eight learned-core scalar centers and six whole scalar centers after actual
receipt/geometry/method checks. Reuse the two completed normalized GGR whole
analytic gradients only if exact coordinates/state/checkpoint agree.

Initial new tasks: eight core analytic MACE gradients and four alpha whole
analytic gradients (12), plus216whole scalar grid points (36x3proteinsx2metals).
If either GGR gradient cannot be matched, report unsupported rather than silently
invent or recompute it under this inventory. GGR whole curvature is shared by
both representations because the physical full system is identical.

Only eligible predictions receive native validation: at most8new analytic DFT
endpoints and16new matching core/full analytic MACE endpoints. Every predicted
point is frozen before these outputs. Do not reuse a full displaced point merely
by proteinID: the two GGR representations may predict different displacements.
No new initial DFT, numerical DFT derivative, trajectory or quantum optimization.
No broader benchmark is launched by this plan.

Existing whole MACE energy calls cost~5s foralpha/~13s forGGR; GGR gradients
~41s. The initial216grid calls imply roughly28GPU-model-minutes total before
startup, distributable among existing allocations. These are development costs,
not routine throughput. Existing native8endpoint validation suggests minutes
on64CPU. Record actual complete costs, failed attempts, prep and reuse separately.
No project CPU/time cap. Use existing manifested runners and standard resources.

## Frozen checks and decisions

- Real source/paired/state invariants, exact accepted receipts and native model
  readout closure. Recomputed center energies agree with scalar archives within
  0.01modelkcal. Gradient adapter must retain its existing core qualification.
- Initial gradient odd-axis changes: error<=max(0.01modelkcal,1%ofprediction)
  at0.02A; half-step odd error must not exceed the same tolerance. Refinement
  conditions above apply before new DFT. All axes are included.
- Native check at each admitted point: actual H(u)-H(0)<0; prediction error
  <=max(0.05kcal,25%of|DFT(u)-DFT(0)-g_DFT.u|), exactly the preceding bounded
  response rule. Actual hybrid gradient combines native and learned derivatives.
  Interior residual norm<=max(2kcal/A,25%ofinitialnorm); boundary requires
  u.g_actual<=0 and tangent norm below the same tolerance. Keep full vectors.
- GGR final hybrid partition difference and response-correction partition
  difference must each have absolute value<=2kcal. Never select a favorable core.
- Report all four alpha-minus-GGR contrasts, actual exploratory and fully qualified
  fields separately. Primary usefulness requires every contrast>0.02kcal and
  all native/numerical gates; no threshold refitting. Both structures remain one
  alpha biological group; no prospective or broad validation claim.
- No newly calibrated decision or production promotion. If useful, the next
  independently declared step must test transfer/robustness, not declare the
  full goal achieved from this consumed pair.

## Status

Plan declared; no new model/DFT calculation yet. Next implement preparation,
frozen preflight, finite initial tasks, assessment and conditional native checks.
