# Mex finite response: both own steps descend; substantial load remains

Job12717 passed all four declared cells and both own-descent checks. Each metal
selects its own proposal from the shared origin/La-proposal/Dy-proposal menu.
The conditional Dy−La accommodation change is **−0.102379 kcal/mol**. This is a
small finite-response result at fixed three-Nd occupancy, not a protein affinity
or stationary-minimum result.

## Actual paired work

All values are kcal/mol relative to the same metal's original state.

| Geometry | La work | Dy work | Common work | Dy−La work |
|---|---:|---:|---:|---:|
| La proposal | −0.252489 | −0.260656 | −0.256572 | −0.008167 |
| Dy proposal | −0.233258 | −0.354867 | −0.294063 | −0.121609 |

Own linear predictions were −0.293230 for La and −0.380224 for Dy; actual work
was less favorable by0.040741 and0.025357, respectively. The La-path differential
changes from a tiny positive tangent prediction to a tiny negative finite result;
the finite paired energies, not the linear differential, determine the result.

Own vacuum MACE gains (La−0.487047/Dy−0.690275) are opposed by fixed-charge ALPB
(+0.226540/+0.301963) and smaller field-response costs (+0.032253/+0.058328).
The existing MM ledger contributes small favorable work (−0.024234/−0.024884).
This supports coherent finite descent while retaining meaningful cancellation
between the model components.

Both physical proposals pass the original covalent, protected-group and contact
guards. Actual maximum displacement, at the target, is0.013166Å for La and
0.013929Å for Dy; largest torsions are0.039233° and0.104282°. Fixed atoms remain
exact, and no severe contact worsens under the declared source-relative rule.

## Remaining own-state loads

These are current-axis energy derivatives divided by the **original** physical
Jacobian column norms, in kcal/mol/Å. They are not full Cartesian force maxima.

| Own endpoint | Leu83 hinge | Asp84 | Asp86 | Glu95 | Target x | Target y | Target z |
|---|---:|---:|---:|---:|---:|---:|---:|
| La proposal / La | −5.704 | +2.005 | −3.525 | −0.275 | +9.645 | −9.493 | −1.972 |
| Dy proposal / Dy | −5.622 | −9.098 | −2.942 | −1.465 | +8.466 | −15.092 | +1.106 |

The seven-coordinate metric retains full rank. Residual target, hinge and donor
loads justify a bounded own-energy continuation in this same basis. They do not
justify adding modes or changing the physical model before testing that response.

## Finite next scope

Root has authorized a separate declared continuation from the exact accepted
12717 endpoints and saved gradients. Reuse the seven modes, original-source
reference/whitening and unchanged physical guards. Freeze the existing LBT
source envelope of0.20Å for all atoms/target and0.20rad per total torsion;
retain0.02Å/0.5° per-trial limits. Use original own-energy Armijo/backtracking,
at most12 accepted steps and24 new own-composite evaluations per metal, then
at most two opposite-metal terminal evaluations. Keep original common3 and the
full origin/two-pilot/two-terminal common5 denominator. No sign-triggered retry,
new mode, parameter adjustment or extra interface campaign follows.

Maximum new work is50 composite cells (50MACE,100localGFN2,50MM configurations,
50ALPB), usually reduced by stopping and exact reuse. A one-warm-GPU/12CPU run
lasting minutes is an extrapolation from the pilot, not a measured saving or
promised wall time. Root owns scout/remaining submissions and final resource
records; the continuation source is separate from this completed pilot.

## Accounting and immutable evidence

Actual new work:4MACE,8localGFN2,4MM configurations/16 group queries,4ALPB;
zero whole-protein SCC, DFT or optimizer iterations. Complete candidate gradients
were assembled from these same calls. Application runtime was30.240890s on one
H200 with12allocatedCPUs:362.890685 application-allocated CPU-s and30.240890
application-allocated GPU-s. These are application intervals. Scheduler accounting separately records32seconds,
384allocatedCPU-seconds and32allocatedGPU-seconds; measuredCPU57.830seconds,
MaxRSS1696468K. See `mex_response/ACCOUNTING.txt`. Observed GFN2 child CPU time was40.750293s.

All paths below are relative to the repository:

- `workspaces/lanm_coupled_response_20261016/mex_response_v1/execution/RESULT.json`, SHA256 `2c5c50010f502ee0a09262cefa8934b5237716dc7290fe6a8cbc187a103eef32`.
- Same directory: `MODES.json` (`b2f227d5`), `PROPOSALS.json` (`b7dc055a`), `PRECHECK.json` (`58ff88b9`), and `{La,Dy}_GEOMETRY.json`.
- Own endpoint gradients: `cells/La_proposal_La/ASSEMBLED_GRADIENT.json` (`a960fe6c`) and `cells/Dy_proposal_Dy/ASSEMBLED_GRADIENT.json` (`ed9b63a7`); both cross-state gradient receipts are also complete.
- Frozen source `MEX_RESPONSE.py` SHA256 `f11777c7c49706c6798f4534ed7bf70b908f1fc1dbf8d815de0d39e6283b35d0`; config `c84d7b29`; plan `ee17c5f4`.

This report transcribes and interprets saved results. No new numerical calculation
or molecular call was performed in preparing it. Earlier failures and all four
state/candidate slots remain preserved.
