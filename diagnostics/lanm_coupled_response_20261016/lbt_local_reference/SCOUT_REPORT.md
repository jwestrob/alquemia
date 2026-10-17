# PBE0 local-reference scout: SCF evidence

**Job 12826 completed and 7CCO/La is admitted under the declared energy-based
SCF policy.** This is one endpoint; no same-metal source strain or four-cell
method comparison is available from the scout alone. No gate was changed.

The admitted PBE0-D4 energy is **−1347.308628960856 Eh**. The exact prepared
48-atom, charge −1 source, La ECP46 and restricted 200-electron state
(100 alpha / 100 beta) passed the existing parser. The analytic gradient's
inventory, coordinates and energy matched. These checks do not establish
wavefunction stability or force consistency.

## Which residuals describe the final state?

The raw output has DIIS iterations 1–6, followed by SOSCF iterations 7–13.
Iteration 13 is followed by `Energy Check signals convergence` and native
convergence after 13 cycles. The configured policy is `Total+1el-Energy`,
`ConvForced 0`, with total-energy tolerance 1e-9 Eh and one-electron-energy
tolerance 1e-6 Eh. A separate final one-electron-energy **change** is not
printed in the terminal residual summary.

| Printed terminal quantity | Value | Printed tolerance | Interpretation |
|---|---:|---:|---|
| Energy change | 5.6480e-10 Eh | 1e-9 Eh | Matches final SOSCF iteration; passes. |
| Maximum density change | 8.0279e-5 | 1e-8 | Matches final SOSCF iteration; above tolerance. |
| RMS density change | 2.0598e-6 | 1e-9 | Matches final SOSCF iteration; above tolerance. |
| Orbital gradient | 2.7042e-5 | 2e-6 | Matches final SOSCF iteration; above tolerance. |
| Orbital rotation | 3.5520e-5 | 2e-6 | Terminal reported value; not separately shown in the iteration table. |
| DIIS error | 1.6338e-3 | 1e-8 | Historical DIIS-related bookkeeping, not the final SOSCF gradient. |

The last DIIS row actually prints `DIISErr 5.65e-2`. The terminal value
1.6338e-3 instead agrees at printed precision with the first SOSCF row's
`MaxGrad 1.63e-3`. Thus it should not be described as either the final active
SOSCF error or simply the last tabulated DIIS error. The raw output does not
establish the exact internal update responsible for that retained field.
Even setting this historical field aside, the current density and orbital
residuals exceed their printed tolerances. The saved
`all_printed_residuals_pass=false` is therefore a substantive limitation.

## Gap warning and final energy

The `Small HOMO/LUMO gap (0.039)` warning occurs after iteration 1, while
DIIS is starting. It triggers a full diagonalization instead of
pre-diagonalization. It does **not** describe the final frontier orbitals:
the final table prints occupied orbital 99 at −0.080834 Eh (−2.1996 eV) and
unoccupied orbital 100 at +0.036700 Eh (+0.9987 eV). No final gap was newly
calculated for this report. These entries do not themselves prove stability.

After SCF convergence, ORCA recomputes exchange on `gridx3`, reporting an
exchange-energy change of −0.000049102 Eh and total SCF energy
−1347.26256434794686 Eh. D4 contributes the printed −0.046064613 Eh.
The reference comparison must use the admitted final PBE0-D4 energy, not an
iteration energy. The PModel atomic La guess's assumed doublet is a separate
atomic initialization; the admitted molecular calculation remains the
restricted singlet stated above.

## Consequence for the source-strain comparison

The endpoint can enter the frozen same-metal energy comparison. Its tiny last
energy change is not a bound on remaining energy error, and no cancellation
of residual errors between sources is established. Keep the residual and gap
context beside each remaining endpoint. This scout cannot alone establish
better accuracy than MACE, a stable physical electronic state, or suitability
of the analytic gradients for optimization. The original finite comparison
continues under its declared gates; this report requests no rerun or expansion.

## Actual cost and evidence

Slurm: **COMPLETED, exit 0:0; 839 seconds, 12 CPUs, no GPU; TotalCPU 02:44:49**.
ORCA prints **835.758 seconds** total runtime. The runner's timestamp interval
is **835.90676 seconds**; application time is **838.359499 seconds**. These
measure different intervals and are retained separately. There was one new
DFT endpoint, with no new MACE/GFN calls or geometry optimization.

[SCOUT_RESULT.json](SCOUT_RESULT.json) copies the execution result exactly.
[SCOUT_ACCOUNTING.json](SCOUT_ACCOUNTING.json) pins the scheduler, runtime,
native output and analytic-gradient receipts. SCF sections are at lines
1049–1056, 1102–1175 and 1279–1283 of the
[raw output](../../../workspaces/lanm_coupled_response_20261016/lbt_local_reference_scout_v1/execution/7CCO_La/endpoint.out).
