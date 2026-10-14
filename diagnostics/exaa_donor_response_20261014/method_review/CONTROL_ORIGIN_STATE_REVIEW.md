# Control La origin: energy-converged, gradient criterion not met

Read-only review of job 4119 and its frozen collector. **`complete_diagnostic` correctly applies the declared solver-specific admission rule. It does not establish gradient convergence equivalent to the two TRAH scouts.** The control has 16 passing recorded checks; it has no `final_TRAH_gradient` check. This is a qualification asymmetry in the frozen protocol, not a newly introduced parser error or a reason to silently rewrite the completed result.

## Actual ending

Raw file: `workspaces/exaa_donor_response_20261014/reference_v1/tasks/sample0_La_origin_vacuum/endpoint.out`.

SOSCF starts at line 1786. Its final iteration 47 prints MaxGrad **4.47e−5**, followed by **“Energy Check signals convergence”** and SCF convergence after 47 cycles (1801–1806). There is no subsequent TRAH macro or smaller final gradient. The final summary confirms **orbital gradient 4.4651e−5 versus tolerance 1e−5**, and orbital rotation **3.9244e−4 versus 1e−5** (1844–1845). Unlike a pre-TRAH DIIS summary, this gradient is direct evidence from the ending solver.

The configured native convergence mode is `Total+1el-Energy`, with `ConvForced 0` (1686–1687). The energy change passes its printed tolerance, while density, DIIS, gradient and rotation summary criteria do not (1840–1845). ORCA nevertheless reports normal termination (6232). The final scalar **−6073.032935684248 Eh** (6031), molecular charge −2/singlet/828 explicit electrons (1625–1627), identities, basis/ECP, populations and retained GBW pass the existing collection checks.

## Why the collector admits it

- Frozen `PLAN.md:41–47` requires native normal termination/convergence and identity/state/artifact checks. Its additional actual-gradient requirement is explicitly **“For TRAH”**; other density residuals and state limitations must remain visible.
- Pinned `implementation/OLD_REFERENCE.py:147–159` enforces native receipt success for every endpoint and adds `final_TRAH_gradient` only when a TRAH macro exists.
- Pinned `implementation/COLLECT_V2.py:51–63` corrects that TRAH threshold parser and skips non-TRAH endings. The ExaA wrapper adds provenance checks without adding a universal gradient gate (`implementation/REFERENCE.py:141–157`).

Thus this is an energy-converged conditional diagnostic under the recorded rule. Calling it “all 17 checks passed,” “all endpoints meet 1e−5 orbital-gradient convergence,” or numerically qualified to the same residual standard as the scouts would be incorrect. Its collection already retains the excess residuals and explicitly marks stability/repeat/gradient qualification unestablished.

The residual alone does not quantify an energy error or prove a different electronic state. The single control origin also supplies no completed control relaxation work. Keep any later response comparison conditional on this weaker convergence evidence; a uniform stricter qualification would be a separately declared follow-up, not a retrospective threshold change. No input, code, threshold, job or result was modified for this review.

## Pins

All paths below are relative to `workspaces/exaa_donor_response_20261014/reference_v1/`.

| Record | SHA256 |
|---|---|
| `attempts/4119_execute_sample0_La_origin_vacuum/COLLECTION.json` | `2385c29923a3d42cf056060676e38dbc66236aa6535cd9c66b4c84c7c86b14c0` |
| `tasks/sample0_La_origin_vacuum/endpoint.out` | `6c8b8d901bafc261fe89cd5de9bfac06e55a9871831195cdbe2e369a3e3cb850` |
| `tasks/sample0_La_origin_vacuum/endpoint.out.execution.json` | `8b5be431527a11cb172b0ce403adc17a73b2051714cb590f7da9d378888e3b13` |
| `implementation/OLD_REFERENCE.py` | `083e430cc80ad94cc71bde063085e26c3020b46174fcd7c23fb86b573d30c395` |
| `implementation/COLLECT_V2.py` | `24182c847f9fb76359a93b784eb1731a236b7306922fdfbcba54a6700509f412` |
| `implementation/REFERENCE.py` | `7ec03787224f05785741c611a3700e8ad470e0e7c7f693dc876dd6712238e894` |

Native reported runtime: **1 h 1 min 3.438 s** (6233). Root owns scheduler cost accounting and the unchanged running endpoints.
