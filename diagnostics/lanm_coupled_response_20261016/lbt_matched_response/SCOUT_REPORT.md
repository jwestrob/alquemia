# Matched-objective LBT scout: completed

Job **12762** completed successfully on 17 October 2026 UTC. All four exact
saved starts passed the zero-call preflight. The declared **7CCO/La** search
retained an admitted endpoint after **8 accepted steps**, including 7 fallback
steps, and **12 new MACE plus 12 native fixed-charge ALPB evaluations**.
No GFN2, MM, DFT or initialization molecular calls were made.

## Own-state response

The same whole 256-atom peptide, charge −2, La multiplicity 1, ten physical modes,
original-source metric, contact/covalent guards and bounds were retained.
The starting gradient was reconstructed from saved raw MACE forces and the
12720 fixed-FF ALPB gradient; the old SCC-solvent composite gradient was not used.

| Energy change, kcal/mol | From exact matched warm start | From matched original source |
|---|---:|---:|
| MACE vacuum | +0.1283229405 | −3.6855057921 |
| Fixed-FF ALPB | −1.0065728860 | −1.1282054282 |
| Composite | **−0.8782499455** | **−4.8137112204** |

Thus the additional own-state decrease is solvent-driven and opposed by a small
vacuum-energy increase. These are saved work values, not newly calculated scores.
The endpoint is `evaluation_010/POINT.json`; later evaluated trials are retained
as audit records and are not additional pool candidates.

## Stop and physical limits

Both the old and new searches report `no_admissible_decreasing_trial`.
At the new terminal iteration, all four normal trials violated the unchanged
0.20-radian torsion bound. The permitted target-x and A10 N–CA fallback trials
passed geometry but raised the composite energy by +0.2156595085 and
+0.0151163597 kcal/mol, respectively, and failed Armijo admission.

The retained endpoint passes all physical checks: maximum displacement from the
original source **0.2999694178 Å** within 0.40 Å, target displacement
**0.1602141683 Å** within 0.20 Å, and maximum torsion **0.1999154017 rad**.
No new/worsened severe contacts, chirality, covalent or protected-group failures
were recorded. The torsion is close to its bound despite the retained endpoint's
strict boundary flags remaining false.

The maximum active-span force falls from **12.75715608** to
**2.47629453 kcal/mol/Å**; final active-span RMS is **0.14575824 kcal/mol/Å**.
The low-force gate remains false. This is a finite constrained endpoint, not a
stationary minimum or proof that no smaller feasible step exists.

## Interpretation and cost

The result admits the unchanged remaining-three stage under the declared plan.
It does **not** establish a paired La/Dy preference gain, changed Mex/LBT ordering,
source consistency, affinity or biological accuracy. Those require the declared
four-state denominator and shared paired pool; no sign-based continuation or
expanded bounds follow from this scout.

Actual scheduler cost: **27 seconds, 12 CPUs and one H200 allocation**, hence
**324 allocated CPU-seconds and 27 allocated GPU-seconds**. Slurm reports
`TotalCPU=00:22.712` and batch `MaxRSS=1311276K`. Application elapsed time is
25.24395381 seconds, with 12 actual model forwards and one warm model load.
Allocated resources are not measured utilization.

## Evidence

- [Actual result](../../../workspaces/lanm_coupled_response_20261016/lbt_matched_response_scout_v1/execution/RESULT.json), SHA256 `a47736fe052c039dfde7843ab29adf83b01b72e872313e3b4e437f729730fcdb`.
- [Search and terminal trial records](../../../workspaces/lanm_coupled_response_20261016/lbt_matched_response_scout_v1/execution/7CCO_La/SEARCH.json), SHA256 `9f7df62f12d0516d8d9b5ca2f84d0ce7ed136bc21858474edac7da7a09044388`.
- [Accounting and exact receipt pins](SCOUT_ACCOUNTING.json), [raw scheduler receipt](SCOUT_SACCT.psv), [completion-event snapshot](SCOUT_COMPLETION_EVENT.json).
- Source `a6cc6de4`, actual configuration `c20b989b`, plan `5a0d805c`; full pins are in the accounting record.

This report reads immutable saved outputs and scheduler metadata only. No new
scientific computation, molecular evaluation, geometry change or submission.
