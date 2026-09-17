# Intact MACE passes the three declared development comparisons

The complete five-structure experiment passes its numerical and relative-order
gates. Both alpha-lactalbumin structures now rank above GGR, where the primary
core descriptor failed, while XoxF remains above MxaF. This is an improvement on
consumed development cases, not independent validation of broad affinity.

| Whole prepared chain | Atoms | R_coord, kcal/mol |
|---|---:|---:|
| GGR, 1GLG | 4698 | 62.714845 |
| Alpha-lactalbumin, 1F6S | 1932 | 99.868407 |
| Alpha-lactalbumin, 6IP9 | 1898 | 90.803133 |
| MxaF, 1H4I | 9088 | 101.614126 |
| XoxF, 4MAE | 8854 | 110.647960 |

Frozen differences, each required above 0.02 kcal/mol: XoxF−MxaF **+9.033834**;
alpha1F6S−GGR **+37.153563**; alpha6IP9−GGR **+28.088288**. The two alpha
structures are one qualified biological affinity comparison. PQQ functional
association is a separate evidence stratum. All five scores are available.

## What this measures

For each metal, subtract the energy after moving only that metal beyond the
model's neighbor cutoff from its bound energy. Then subtract the La difference
from the Ca difference, converting eV once. Each subtraction has identical atoms,
charge, spin and protein coordinates. This cancels geometry-independent atomic
and global readout offsets. Component audits pass; full terms and unrounded
energies are retained in the [result index](INTACT_RESULT.json).

The method uses intact prepared chain A, fixed cofactor/waters and native
MACE-OMOL-0 100M float64 energies. Exact edge batching reproduces the native
reference; the large PQQ calls also require PyTorch's expandable-segments
allocator. Two actual OOM attempts remain archived. No DFT, solvent, optimization,
training, force or entropy calculation was added. Production remains unchanged.

This descriptor is not a solution binding free energy. Its finite-range model
does not provide full protein electrostatics or certified separated ionic states.
All proteins exceed reported training sizes. No absolute decision bands or
universal zero are justified. Whole context and total-charge conditioning changed
together, so the improved ordering does not isolate a physical mechanism.

## Measured cost and remaining work

For the four newly scored proteins, jobs 1200815/1200816 made 16 successful calls
and one failed call, reusing the four qualified alpha1F6S endpoints. They cost
**1387 GPU allocation seconds, 22,192 allocated core-seconds and 986.322 actual
CPU-seconds**. Individual PQQ calls take about 26 seconds of inference; four
calls require about 104 seconds, or 109–111 seconds including worker setup.
Peak allocated GPU memory is 21.72 GB on the A5000. These worker timings exclude
preparation and recursive verification. Final report verification alone took
803.92 wall seconds / 485.49 CPU seconds; reducing this overhead remains work
before a convenient routine scorer is ready.

Including the earlier one-time native/adapter engineering and the next product
core qualification, the checkpoint totals are 68 successful calls, two OOMs,
2007 GPU-seconds, 95,344 allocated core-seconds and 19,534.325 actual CPU-seconds.
Local preparation/report costs have separate receipts; other local work is not
fully profiled. See [cost/status record](INTACT_ENGINEERING_STATUS.json).

Next: finish [exact product batching](EXACT_PRODUCT_PLAN.md), then execute the
[declared canonical extension](INTACT_CANONICAL_PLAN.md). The readiness audit
matched 27 protein templates and identified missing terminal chemistry in 1KB0;
four calibration proteins have charges outside the reported training range.
Those limitations are retained explicitly, without charge changes or omitted
denominators. The broader goal is active; this candidate is not promoted.
