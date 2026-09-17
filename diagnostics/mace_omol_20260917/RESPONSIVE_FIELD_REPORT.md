# Responsive quantum density: useful solvent effect, failed boundary gate

**Retain the production baseline.** The new candidate passes3/4 development
ordering tests, compared with1/4 for the frozen-density predecessor. It still
fails the unchanged GGR partition criterion. No calibrated decision, production
promotion or broad validation follows from this small consumed panel.

## Paired result

All inputs, waters, assembly, source mapping, permanent environment charges,
MACE short components and solvent parameters were frozen before execution.
Only the core density responded to the permanent protein field; the solvent
term was recomputed from that new density.

| Alpha minus GGR, kcal/mol-scale | Frozen density | Responsive density | New direction gate |
|---|---:|---:|---|
| 1F6S minus extended | -2.523559158 | +2.795344719 | pass |
| 1F6S minus connected | +1.615626846 | +7.233492066 | pass |
| 6IP9 minus extended | -7.739380959 | -2.107829947 | fail |
| 6IP9 minus connected | -3.600194955 | +2.330317400 | pass |

The predeclared requirement was all4differences>0.02. GGR connected-minus-
extended changes from-4.139186004 to **-4.438147347kcal**, failing abs<=2.
Selecting the connected core or the favorable alpha structure would hide that
failure. Alpha water inventories remain2versus3; do not call their difference
pure structural noise. These are two biological groups, not four independent
observations: GGR has direct same-assay evidence; the alpha comparison remains
qualified cross-study evidence. Every case was already consumed for development.

## What changed the result

| Representation | Electronic response contribution to R | Updated GB contribution to R |
|---|---:|---:|
| GGR extended | +0.026715730 | -3.176083552 |
| GGR connected | -0.176729340 | -3.271599825 |
| Alpha1F6S | +0.016360626 | +2.153175429 |
| Alpha6IP9 | +0.040415310 | +2.441767880 |

R=E_Ca-E_La; values are changes from the exact frozen-density predecessor.
The direct energetic response almost cancels between metals. But the updated
charge distribution changes its interaction with the solvent boundary enough
to shift alpha-minus-GGR comparisons by about5.3–5.9kcal. **A small differential
quantum response energy does not justify reusing the old reaction-field term.**
The saved component algebra independently reproduces each complete score update.
This identifies a useful contribution within this approximation, not a unique
physical explanation of experimental selectivity or a validated accuracy gain.

All8variational diagnostics are nonpositive within the fixed+.05kcal allowance.
Both endpoint and paired exterior CHELPG/projection checks pass. The actual
near-boundary projected-minus-exact Ca/La coupling errors are-0.500830,
-0.542104,-0.348571,+0.325649kcal foralpha1F6S,alpha6IP9,GGRconnected,GGRextended.
The GGR difference is-0.674220. All5fixed1kcal checks pass. No charge scheme,
reference, dielectric, cap, water or threshold was changed after these outputs.

## Model and numerical scope

Final protocol `embedded_QM_ff19SB_full_OBC2_POLAR_short_hybrid_v1`:

`B_M = E_DFT,embedded(core,M; q_env) + G_full(Pq_embedded,M + q_env) + T_short(full,M) - T_short(core,M)`.

The native quantum energy already includes direct core–environment interaction.
Do not add exact coupling a second time. `%method DoEQ false` excludes external
charge self-energy. Native r2SCAN-3c/DefGrid3/TightSCF, basis/ECP/D4/gCP and
analytic gradients are verified in the actual outputs. No CPCM was added.
Environment-only terms cancel in the endpoint contrast; all GB self/cross terms
remain auditable. The unchanged physical full boundary uses the existing OBC2
parameters and commonCa/La radius. The MACE short/full terms are exactly reused.

This is **one-way coupling**: the core responds to permanent protein charges,
but does not self-consistently respond to the added GB field. Classical protein
charges remain fixed. Projected monopoles approximate the embedded density;
MACE's learned short/field split is not a unique physical decomposition.
No claim of self-consistent QM/GB, complete protein polarization, or binding
free energy is made. Combined gradients, relaxation, absolute affinity score,
aquo reference and calibrated class remain null. Core analytic gradients are
exported; they do not differentiate the full candidate.

All50solver/component checks pass, including identity, rigid transformations,
repeat andReference comparisons. Maximum solver discrepancy0.000273732kcal
versus0.01tolerance. These are component checks, not newly executed quantum
rotation or combined-gradient qualifications. The failed partition gate is a
separate physical representation limitation.

## Real execution, recovery and tests

-1200983:8native quantum endpoints, noSCFretry,655s×64CPUs=41920core-s.
-1200984:utility preflight failed before any fitting/potential evaluation;
  25s×8CPUs=200core-s, retained.
-1200985:read-only replay on the same worker proved1.11e-16cap-weight roundoff;
  1s×1CPU. Source atom identities, recorded weights and scientific thresholds
  stayed fixed. Replay now permits1e-12arithmetic error in that scalar metadata.
-1200986:8CHELPG+8potential utilities completed;294s×8CPUs=2352core-s.
-1200989:44GBcalls completed,4identical environment-only terms reused;
  66s×16CPUs=1056core-s and **66GPU-allocation-s**. Zero newMACE/training calls.

**Total:45529allocatedcore-s;37362.387reported actualCPU-s;66GPU-s.**
Utility wall times sum1055.487128s over concurrent processes; peak individual
RSS354328KiB. GB solves sum3.301663s; peak process RSS241664KiB; GPU memory
was not measured. Quantum peak batch RSS16047704KiB. Local preparation,
failed preflights, testing and collection resource files are separately retained.
Prior quantum, source preparation and reused MACE costs are additional. This
is a development allocation accounting, not a matched end-to-end production
speedup. The quantum stage remains in the same cost range as the recent matched
vacuum8endpoint job (751s/48064core-s); no larger cost regime was introduced.

21distinct real-fixture/receipt/parser/component tests passed across the
focused runs:5responsive quantum,8responsive component,4existing vacuum,
4existing charge tests. Actual-result tests initially skipped until their
scientific outputs existed, then passed. No dummy scientific results.
A legacy error-message regression was fixed before launch; the first charge
preparation caught a local variable-name collision before execution. The initial
collector also mistook a conditional ORCA warning for numerical-gradient use.
Its exact-warning repair was checked against real analytic outputs and a
clearly corrupted copy carrying a numerical-gradient header. No scientific
rerun or acceptance-threshold relaxation was used for these technical repairs.

## Assessment and next action

**Numerics:** executed components and charge representation pass; partition
robustness fails. **Scientific information:** the consistent solvent update
moves the development ordering substantially, but robust predictive utility
is not established. **Affordability:** this contained experiment fits the
existing endpoint cost range; full production-pair cost is unmeasured.

Recommendation: retain baseline and do not expand this frozen candidate into
benchmark calibration. Further work should address the remaining model and
boundary limitations rather than tune this panel's charges or thresholds.

[Full expression and frozen scope](RESPONSIVE_FIELD_PLAN.md),
[native accounting and primary documentation](RESPONSIVE_FIELD_ACCOUNTING.md),
[executable operations](RESPONSIVE_FIELD_COMMANDS.md),
[compact unrounded result](RESPONSIVE_FIELD_RESULT.json).
All native outputs, wavefunctions, receipts and failed attempts are retained
underworkspaces/mace_omol_20260917. The active MACE goal remains incomplete.
