# Adaptive candidates and common-pool expansion: independent checkpoint

2026-09-22. Read-only review of actual pilot1209857 and the new pool adapter.
No additional molecular calls, optimization or job submission by this reviewer.

## Adapter review

**No remaining blocker found.** Actual expansion manifest
`workspaces/adaptive_accommodation_20260922/common_pool_v1/manifest.json`, SHA
`db863273ac07019308fc53f5ab347d4c2d98a5965a362c5bf53bb36d79d25fb6`, passes validation:
four prepared cases, five physical geometries each (origin, previous Ca/La
proposals, adaptive Ca/La proposals). Sixteen composite cells need solvent work;
eight own-metal MACE values are reused exactly, leaving eight new MACE and32GFN2
singlepoints. Actual Ca/La origin components match the old pool exactly.

Reviewed physical source/charge/multiplicity and model reuse, mapped candidate
geometry, unchanged inactive coordinates, exact native receipt/XYZ checks,
deterministic candidate ordering and the existing Ca-minus-La algebra. The new
adaptive pool has its own protocol ID and no inherited new calibrated reference.
Unavailable candidates remain unavailable, with the valid old pool retained
separately. No missing solvent energy is replaced by zero.

Two edge cases were found and fixed by root before this manifest:

1. Numerical-copy deduplication permits1e-12-Angstrom coordinate differences;
   native receipt reuse requires exact representative coordinates. Nonexact
   aliases now receive a fresh representative evaluation, with explicit status,
   rather than incorrectly invalidating a case.
2. Partially staged cells from an unsupported case could otherwise reference
   discarded tasks during collection. Unsupported cases now short-circuit to
   unavailable before those task lookups, preserving the full denominator.

## Actual candidate physics

All eight endpoints report successful SLSQP and pass final source/cap/bond/fixed-
atom checks. Seven candidates are inside both bounds; PLM8344 La reaches the
0.8-Angstrom physical boundary at extra-Asp307 OD2. Its selected normalized
residual is2.460 kcal/mol/angstrom; all other selected residuals are at most
0.000691. Omitted angular loads remain1.546–13.083. These are vacuum-OMOL
derivatives, not composite gradients or complete stationarity qualifications.

The31 completed infeasible **trial requests** are retained, including a maximum
3.804-Angstrom intermediate displacement for1H4I La. Requests can reuse an existing
evaluation; this is not a claim of31 distinct extra model calls. The final-domain
rule explicitly permitted intermediate infeasibility. No final result was clipped
or projected into success. No angle-bound flag fired.

| Source | Final max displacement Ca / La, A | New native work Ca / La, kcal/mol | Additional native contrast change versus old proposal |
|---|---:|---:|---:|
|1H4I|0.627 / 0.542|-4.662 / -5.543|+1.102|
|4MAE|0.438 / 0.253|-6.061 / -2.516|-2.859|
|PLM8344|0.737 / 0.800|-40.852 / -106.339|+3.731|
|PLM07ab|0.660 / 0.728|-24.146 / -62.562|+0.901|

Work is measured from each original q0; the last column subtracts the already
completed terminal-proposal change. Most PLM relief was already present in that
earlier result. The new motions narrow the native vacuum4MAE-minus-1H4I contrast
gap by3.962 kcal/mol relative to the old proposal. Their actual largest heavy
movements involve1H4I Asp303,4MAE Asn256 and each PLM extra-Asp OD2.

**Next decision:** use the actual solvent-aware expanded common-pool result,
including old geometries, to judge these candidates. Lower native energies or
favorable shifts for unknown PLM labels do not establish improved discrimination.
Do not automatically enlarge the bound or add modes because one endpoint is
constrained. The current result qualifies candidate generation, not utility.

Sources: `proposals_v1/after_proposals_1209857.json` and the eight pinned
`proposals_v1/proposals/TASK_ID/result.json` records in the adaptive workspace;
previous proposal records remain under
`workspaces/accommodation_nonlinear_20260920/proposals_v1/proposals/`.
