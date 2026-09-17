# Rotation sensitivity traced to the realspace representation

Completed 2026-09-16, job **1200396**, under the [approved investigation](AGREEMENT.md).
**The memory rewrite reproduces the original backend, including its rotation
error. The fixed-axis displaced-charge approximation itself breaks rotational
covariance on the actual archived densities.** No baseline/default was changed.

## Evidence

Eight new MACE-POLAR-1 medium core energy/force calls completed: the same four
consumed 1H4I qm33/qm36 La/Ca cores, rotated by the previous 37-degree transform,
with original and blocked execution. All weights, charge/spin, protonation,
water inventory, widths, offsets, precision and physical structure were retained.
The four unrotated results for each implementation were reused. No new DFT,
full-protein, large-checkpoint, relaxation or training calculation ran.

Original and blocked rotated results agree to maximum absolute differences:

- Energy: 7.276e-12 eV.
- Forces: 3.932e-12 eV/Angstrom.
- Density coefficients: 1.004e-13.

These pass the existing implementation-equivalence tolerances by a wide margin.
Original-backend rotation errors are:

| Core | Energy change (kcal/mol) | Maximum force change (eV/Angstrom) |
|---|---:|---:|
| qm33 La | +0.022984469 | 0.004146875 |
| qm33 Ca | +0.033738587 | 0.003593137 |
| qm36 La | +0.035177043 | 0.003066917 |
| qm36 Ca | +0.053098241 | 0.003952725 |

All four fail the previously frozen endpoint energy/force tolerances
(0.01 kcal/mol and 0.001 eV/Angstrom). The paired Ca-minus-La contrast changes
by **+0.010754118** for qm33 and **+0.017921198 kcal/mol** for qm36. All eight
new calculations pass charge closure. Learned densities also fail exact
rotational covariance, with maximum coefficient changes around 0.0005–0.0006.

## Isolating the numerical mechanism

The original realspace backend replaces each dipole with charges displaced
along the laboratory x/y/z axes, using 0.1 Angstrom offsets for field features
and 0.02 Angstrom for energy. That finite charge cloud has extra higher moments
whose orientation depends on those axes. Field projections use forward
finite differences too. [Source audit and derivation](SOURCE_AUDIT.md).

For each of the four archived, actually computed densities, we evaluated the
field and energy modules before and after the same rigid rotation. Dipole
coefficients and expected field vectors were rotated using the installed
backend's [y,z,x] convention. Then we rotated the auxiliary axes as well,
keeping components expressed in that co-rotating basis. This last calculation
rotates the *same* auxiliary cloud rather than rebuilding it on lab-fixed axes.

| Frozen density | Lab-fixed energy change (kcal/mol) | Maximum lab-fixed field error (raw feature units) |
|---|---:|---:|
| qm33 La | -0.006420914 | 0.070602119 |
| qm33 Ca | +0.001778290 | 0.073652243 |
| qm36 La | +0.004218112 | 0.070648652 |
| qm36 Ca | +0.016635701 | 0.073484048 |

Co-rotating the stencil reduces the maximum error across both implementations
to **1.055e-11 eV** for energy and **1.103e-11** in raw feature units. All such
identity checks pass. The 48 module evaluations retain original offsets and
self terms; no fitted parameter or offset sweep was used.

This establishes a representation-level source of orientation dependence.
The feature probes use frozen final total densities, not the successive
spin-resolved densities inside the neural update. They are not an additive
decomposition of the full model's energy error. We have not demonstrated that
this is the only possible numerical defect, nor rerun a corrected full protein.

## What this means for the discriminator

The memory backend is cleared of introducing this tested rotation error. Keep it.
Analytic evaluation of monopole/dipole fields is the next correction to test:
it can remove the finite-stencil approximation while retaining blocked memory
use. Merely rotating the stencil with a chosen reference orientation is a
diagnostic, not a general physical remedy or a new production score.

The rotation error is small compared with the existing approximately 5 kcal/mol
hybrid partition residual. This tested rotation changes the difference between
the two core contrasts by only 0.007167080 kcal/mol; it does not explain away
that partition result. Numerical consistency and predictive accuracy remain
separate questions. No compatible aquo reference, S, class, or accuracy gain
has been established. The baseline remains usable and unchanged.

The [next pilot](NEXT_ANALYTIC_PLAN.md) is proposed, not authorized or executed.
It replaces a numerical approximation and therefore receives a new research
method version; it does not inherit baseline bands or silently alter old MACE
results. Large-checkpoint evaluation remains a separate step.

## Cost, tests and artifacts

- One A5000, 16 CPUs, 64,474 MiB requested host RAM.
- Job elapsed **62 seconds**; **992 allocated core-seconds**; reported actual
  CPU time **87.961 seconds**. No retries or failed GPU attempts.
- Eight core model evaluations total **7.267436 seconds**; subprocesses including
  imports, model loading and representation probes total **59.380493 seconds**.
- Peak GPU tensor allocation **1,283,105,280 bytes**; maximum host RSS
  **1,614,984 KiB**. CPU preparation was not fully profiled.
- **21 tests passed** in 42.880 seconds before submission: the 17 existing
  real-artifact/kernel tests plus four rotation-manifest integrity tests.
  These are distinct from the eight actual GPU calls and 48 component evaluations.
  After adding the rotation-specific report renderer, all five targeted tests
  passed, including a new real-result denominator check (22 unique tests total).

Protocol: `mace_polar_1m_realspace_rotation_attribution_v1`.
Workspace: `workspaces/mace_rotation_20260916/audit_v1/`.
[Compact result with exact values and pins](result.json).
`comparison.json` SHA256:
`185608b291a68287cb8ffdf4c08bdc311419a31009418a3891209c610da636a6`.
Both backend directories retain manifests, implementation snapshots, terminal
collections, forces, densities, raw probe arrays and execution receipts.
`launch.json`, `sacct_final.tsv` and `job_1200396_accounting.json` retain resources.
No job or continuation remains active.

## Runnable operations

From the repository root, using the preserved environment:

```bash
MACE_PY="$PWD/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python"
ROTATION_ROOT="$PWD/workspaces/mace_rotation_20260916/audit_v1"
"$MACE_PY" scripts/mace_hybrid.py dry-run --manifest "$ROTATION_ROOT/original/manifest.json"
"$MACE_PY" scripts/mace_hybrid.py dry-run --manifest "$ROTATION_ROOT/blocked/manifest.json"
"$MACE_PY" scripts/mace_hybrid.py collect --manifest "$ROTATION_ROOT/blocked/manifest.json"
```

The full original and blocked human-readable reports are `original/REPORT.md`
and `blocked/REPORT.md`. `scripts/mace_rotation.py prepare` builds versioned
manifests from explicit original/blocked collections and agreement; `compare`
checks both manifests and writes to a new output path. Execution uses the
existing MACE runner; the exact submitted script is `run_rotation.sbatch`.
No rerun is needed to read or verify these completed results.
