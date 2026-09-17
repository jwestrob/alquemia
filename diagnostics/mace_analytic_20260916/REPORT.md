# Analytic dipoles fix the tested MACE rotation error

Completed 2026-09-16, job **1200470**. All 12 real MACE calls completed on one
RTX A5000. Analytic monopole/dipole evaluation passes the frozen core and
full-protein rotation checks, with essentially unchanged runtime and memory.
The approximately 5 kcal/mol hybrid partition discrepancy remains.

The baseline/default, its references and old experiments are unchanged. This
is an opt-in research protocol:
`mace_polar_1m_analytic_multipole_vacuum_r2scan3c_pilot_v1`.
No calibrated S, class or predictive-accuracy improvement is claimed.

## What changed and what ran

Replaced fixed-axis displaced-charge sums with analytic derivatives of the
same regularized Gaussian radial kernel. Preserved widths, normalization,
self terms, exclusions, charge constraints, float64 weights and all physical
states. The field inputs and resulting learned densities do change; this is a
new numerical model, not a silent migration of old results. [Formulae](KERNEL.md).

Eight primary/rotated core calls used the archived 1H4I qm33/qm36 La/Ca states.
After their numerical gate passed, four primary/rotated full-protein calls used
the same 9,141 atoms. The rotation stayed at the previously declared 37 degrees.
No new DFT, relaxation, solvent variant, biological case or threshold fitting.
Four converged vacuum DFT endpoints were reused for the subtractive comparison.
Repeat/translation were not rerun in this analytic pilot; the previous finite-
displacement campaign's corresponding results remain separate.

## Numerical and scientific results

| Check | Finite displacement | Analytic dipoles | Frozen gate |
|---|---:|---:|---:|
| Maximum core rotation energy error, kcal/mol | 0.053098241 | 0.000002929 | 0.01 |
| Maximum core rotation force error, eV/Angstrom | 0.004146875 | 4.306e-8 | 0.001 |
| Full La rotation energy change, kcal/mol | -1.867693692 | -4.081e-7 | absolute 0.01 |
| Full Ca rotation energy change, kcal/mol | -1.811157196 | -2.148e-7 | absolute 0.01 |
| Full paired Ca-minus-La rotation change, kcal/mol | +0.056536496 | +1.933e-7 | absolute 0.01 |
| Maximum full rotation force error, eV/Angstrom | 0.059996344 | 8.011e-7 | 0.001 |
| Hybrid qm36-minus-qm33 partition shift, kcal/mol | -4.998674456 | -5.023220584 | absolute 2 |

Core/full rotation checks **PASS**. Maximum charge closure error across the 12
calls is 1.2861e-12 e, passing 1e-5 e. Partition check **FAILS**, essentially
unchanged. No tolerance or state was adjusted after inspecting results.

Removing the stencil changes the full La/Ca endpoint energies by -4.679784 and
-4.411285 kcal/mol, respectively. Their common part largely cancels: the full
raw Ca-minus-La contrast changes by +0.268499 kcal/mol. Core endpoint changes
are +0.126 to +0.211 kcal/mol. Maximum full force changes relative to the old
model are about 0.0944 eV/Angstrom; force-based research must use the matching
new energy/gradient pair rather than mixing implementations.

For the new model:

| Core partition | DFT core R | MACE core R | MACE full-minus-core R | Hybrid R |
|---|---:|---:|---:|---:|
| qm33 | -405349.0793037981 | -405420.1762964037 | 86.7191697813 | -405262.3601340168 |
| qm36 | -405287.3407004409 | -405353.4144724627 | 19.9573458403 | -405267.3833546005 |

All entries are kcal/mol; R = E(Ca)-E(La). These large raw elemental offsets
are not affinity labels. The expression is E_full,MACE + E_core,DFT - E_core,MACE
in vacuum, without a compatible aquo reference or inherited baseline bands.
The common full term cancels in the partition difference.

## Actual costs and tests

- One A5000, 16 CPUs, 64,474 MiB requested host RAM; no host offload.
- Full primary evaluations: **58.187702 s La**, **58.630392 s Ca**; pair
  **116.818094 s**, compared with 116.461473 s for the old blocked pair.
- Peak full GPU tensor allocation **10,251,166,720 bytes** (~9.55 GiB);
  peak worker host RSS **1,717,780 KiB** (~1.64 GiB).
- Whole 12-call allocation: **329 GPU-allocation seconds**, **5,264 allocated
  core-seconds**, **363.761 reported actual CPU seconds**. No failed model calls
  or retries. A read-only GPU monitor printed an I/O warning but returned its
  device measurement; it did not affect the model tasks.
- Subprocess total **324.737393 s**; model evaluation total **242.718734 s**.
  Queue wait is excluded. Allocated GPU time is not a utilization integral.
- **29 tests pass**. Four independent analytic/autograd tests first ran in
  12.726 s, then the complete regression suite ran in 26.674 s. They cover
  real-core fields/derivatives, energy/force/feature rotations, mixed-batch
  exclusions, linear saved storage, unchanged weights, manifest admission,
  stale-kernel rejection and missing-result gates. Preparation was not fully
  profiled. These tests are separate from the 12 actual GPU model evaluations.

## Interpretation and next work

**Retain the analytic memory backend for research.** The tested numerical
rotation defect is fixed at the same practical cost. This improves the basis
for global scoring and structural-response experiments; it does not establish
better biological classification or resolve the hybrid boundary discrepancy.
Continue baseline predictions unchanged.

Under Jacob's standing approval for contained pilots, the next comparison is
MACE-POLAR-large on the same consumed states: test architecture and core
rotation/partition behavior before full-system execution. This changes model
weights and receives its own protocol/results. Do not infer a successful large
run from this medium-model report. Solvent/reference compatibility and broader
experimental-label comparisons still need separate explicit experiment records.

## Artifacts

[Exact summary and pins](result.json); [runbook](RUNBOOK.md);
[authorization](AGREEMENT.md). Workspace:
`workspaces/mace_analytic_20260916/pilot_v1/`.

`collection_job_1200470.json` contains unrounded endpoint energies, components,
forces/density references, comparisons and all receipts; `REPORT.md` is the full
machine-rendered report. `sacct_final.tsv` and `job_1200470_accounting.json`
retain terminal accounting. The four-core kernel test receipt pins source,
checkpoint, tests and log. No analytic-medium job remains active.
