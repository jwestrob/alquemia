# Mex combined-force interfaces pass the declared finite check

**All 16 displaced cells are complete. Both metals pass both physical directions
at both step sizes, including every energy component, MM group, total, and the
derived common and Dy−La total checks.** This admits the bounded own-energy
response pilot under the implemented Hamiltonian. It is a local derivative
consistency result, not a binding-affinity or chemical-accuracy validation.

## What was tested

The [frozen plan](../MEX_FORCE_CHECK_PLAN.md) tests centered displacements of
±0.0001 and ±0.00005 Å at the exact repaired Mex8FNS origin:

- **Leu97 O:** move one real MM/field atom, updating its point-charge position,
  whole-protein MM energy and full-cavity fixed-charge ALPB. The unchanged local
  MACE and isolated GFN2 components are reused exactly and have zero derivative.
- **Leu83 CA:** move the omitted cap anchor transversely, reconstructing cap205
  and its anchor Jacobians. The external field stays byte-identical; local MACE,
  both GFN2 terms, MM and ALPB follow the same physical displacement.

The combined energy is MACE vacuum + embedded-minus-isolated GFN2 + the retained
MM ledger + full-geometry native FF-charge ALPB. Caps are electronic boundary
atoms; the classical and solvent terms use the 2049 real atoms. The tests preserve
the original charge, occupancy, three Nd spectators and all source waters.

## Saved numerical result

Values below are signed **energy derivatives**, in kcal/mol/Å; forces have the
opposite sign. Residual means finite difference minus analytic derivative.

| Physical direction | Metal | Analytic total | Residual at 0.0001 Å | Residual at 0.00005 Å |
|---|---|---:|---:|---:|
| Leu97 O | La | 9.019881125 | −6.844624e−6 | +5.770399e−6 |
| Leu97 O | Dy | 9.006718251 | −4.092599e−6 | +8.461117e−6 |
| Leu83 CA | La | 4.166892401 | −2.232188e−5 | +1.403851e−5 |
| Leu83 CA | Dy | 4.240676704 | +7.201563e−7 | +1.570650e−5 |

The derived Dy−La analytic derivatives are −0.013162874 for Leu97 O and
+0.073784303 for Leu83 CA. Their saved residuals are respectively
(+2.752026e−6, +2.690718e−6) and (+2.304203e−5, +1.667993e−6) at the two steps.
Both common means also pass. All component and step-agreement gates pass the
predeclared additive tolerance 0.2 + 5% of the analytic magnitude; individual
MM groups pass their tighter 0.001 + 0.01% tolerance. No tolerance, step, fixture
or failed slot was replaced after inspection. Full component checks and exact
unrounded saved values are retained in [RESULT.json](RESULT.json).

## Decision and limits

Proceed to the declared small, covalent-preserving protein/target response:
one proposal from each metal's own total gradient, followed by both metals
scoring the same admitted geometry pool. These single-atom derivative fixtures
are not proposed protein motions. The existing
[archived PBE0 comparison](../mex_reference_response/REPORT.md) supplies separate,
limited electronic-response evidence; its quantitative and own-load disagreements
remain relevant.

Passing two interfaces does not qualify every finite displacement, a stationary
minimum, the fixed-charge solvent approximation, f-electron states, source
robustness or experimental affinity. Geometry and execution gates remain necessary
for the next finite proposal; no larger optimization or new DFT campaign follows.

## Cost and provenance

Job **12709 COMPLETED, exit 0:0**, in **69 scheduler seconds**, with **12 CPUs and
one NVIDIA H200**. Actual work: **8 MACE calls, 24 local GFN2 calls, 16 MM
configurations / 64 MM group queries, and 16 native ALPB calls**. No DFT,
whole-protein SCC or optimizer steps. The unchanged local components were reused.

Scheduler allocation was 828 CPU-seconds; measured total CPU was 02:39.977.
Application wall time was 67.855864 seconds; the saved GFN2 child CPU total was
128.009703 seconds. Allocation time is not measured GPU utilization. This report
copies saved results and reads terminal accounting; it performs no new derivative
calculations or molecular calls.

- [Actual result](../../../workspaces/lanm_coupled_response_20261016/mex_force_check_v1/execution/RESULT.json), SHA256 `8788a7738beaa141f7b83d5ea6212dab77511fc290f97ec61e78c9fe6a6a35f0`.
- [Compact result and complete component checks](RESULT.json), [terminal accounting](ACCOUNTING.txt), and [frozen plan](../MEX_FORCE_CHECK_PLAN.md).
- Executed source SHA256 `333c775e76f80b94b8c85069e4b520ae848f35dc37340d8dda866d80898b2bad`; config `ac00a61a`; admitted static parent 12685 `ff7b27d4`.

Production and shared status are unchanged by this report.
