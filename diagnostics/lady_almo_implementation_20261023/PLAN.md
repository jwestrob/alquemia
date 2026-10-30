# Implement and verify molecular ALMO

## Authorization and scientific question

The active user goal is “implement the ALMO constraint as you described.” This
authorizes implementing the published closed-shell AO-block variational method
and testing it, superseding the former diagnostic's ban on new implementation.
The prior CP2K attempt remains terminal and preserved. No force-field fitting,
new geometries for the La/Dy experiment, or production changes are included.

At fixed geometry, minimize the same electronic functional while each fragment's
occupied orbitals retain support only on that fragment's atomic basis functions.
Allow intrafragment polarization and retain full electrostatics, exchange and
interfragment occupied overlap. Do not mask the Fock or density matrix as a
substitute for this constraint.

For real closed-shell occupied coefficients C and AO overlap S, use
D = 2 C (C^T S C)^-1 C^T. The backend supplies the full-system PBE0 energy and
Fock matrix at this density. The optimizer varies only allowed coefficients.
Document analytic derivatives, rank/conditioning checks and convergence.

## Deliverables and evidence required

1. Reusable solver, readable theory, input loader, execution interface and pinned
   isolated runtime. No dependence on a commercial license.
2. Actual Slurm tests of orbital finite differences, electron count/idempotency,
   fragment support, occupied-gauge invariance, one-fragment recovery of ordinary
   SCF and variational energy ordering. Tests must fail on incorrect results.
3. Independent implementation comparison on a fixed two-water fixture: native
   CP2K molecular GAPW/PBE versus PySCF with identical literal 6-31G* coefficients.
   Two predetermined CP2K grids, 400 and 600 Ry at 16 Angstrom, assess numerical
   sensitivity. This is implementation validation, not a change to quartet PBE0.
4. Load and hash-check exactly the archived La/Dy aquo8 and hydrated monodentate
   acetate quartet. Preserve coordinates, PBE0, orbital bases, La ECP46/Dy ECP55,
   effective singlets, metal/whole-shell grouping and archive D4 bookkeeping.
5. Reproduce full molecular electronic energies in PySCF; original absolute
   compatibility tolerance is 1e-4 Eh per endpoint. Record electronic and D4
   components separately. Failed comparison is not silently admitted.
6. Run the validated constrained solver on the real systems, preserving fragment,
   frozen, polarized and full states, raw logs, convergence and actual costs.
   Reuse completed ordinary/fragment results when proceeding. Admit scientific
   double differences only after implementation and numerical gates pass.
7. Final report states what is implemented, tested, unresolved and reproducible;
   distinguish software correctness from physical accuracy of the archived model.

## Numerical and execution choices

PySCF 2.14.0 in a new isolated environment. PBE0 direct Coulomb/exchange, no
density fitting or seminumerical exchange. Molecular quadrature initially
150 radial by 974 angular points per atom, unpruned. Ordinary SCF energy
tolerance 1e-10 Eh, orbital-gradient tolerance 1e-7, maximum 300 cycles.
This is an explicit backend change; it must pass cross-code energy checks.

All numerical work uses Slurm, account banfield, partition standard,
QoS standard_expand. Root owns submissions. Initial ordinary quartet runs
concurrently with eight CPU threads and 32 GiB per system. Independent solver
tests and native comparison run in parallel as ready. Native comparison may
use up to 64 GiB per allocation. Finite scheduler limits are recorded with each
job; no job is cancelled merely because its output is quiet. Successful prior
work is reused, failed attempts retained, and code changes create new snapshots.

## Responsibilities

- Root: integration, runtime, submission, evidence review and final reporting.
- Optimizer agent: almo.py and theory.
- Archive agent: exact basis/ECP/geometry loader and molecular runner.
- Independent agent: tests and native comparison fixture.

No scientific result or implementation success is claimed by this plan.

## Explicit Dy initialization recovery after actual startup failures

Jobs 22328, 22330 and 22335 failed in PySCF MINAO's nonsquare ECP/ANO overlap
determinant before an SCF iteration. Retry exactly those two Dy baselines and
one shared Dy metal fragment using PySCF's native `1e` core-Hamiltonian starting
density. Each retry retains 8 CPUs, 32 GiB, 2 hours and every electronic/input,
quadrature, convergence and archive-parity requirement. Bind the original failed
RESULT, INPUT and START plus unchanged PySCF source in the retry declaration.
Other jobs continue and their completed work is reused. This is an initialization
recovery, not a different functional, ECP, basis, spin state or relaxed gate.

## Bounded diagnosis of the actual La archive mismatch

La-aquo 22327 converged with orbital gradient 1.1585e-8 and all native density
checks passing, but differs from the archived electronic energy by
+0.000825580291 Eh (0.518059454 kcal/mol), failing the unchanged 1e-4 Eh gate.
Inspect the saved density rather than repeat the SCF. The question is whether
the explicitly documented PBE constants account for this cross-code difference.

One standard_expand job, 8 CPUs / 32 GiB / 1 hour, will evaluate exactly three
energy/Fock pairs on that same density and grid: native PBE0, its explicit
LibXC component expansion, and that expansion with ORCA's documented
mu=0.219520, kappa=0.804, beta=0.066725. The first two must agree and reproduce
the saved result. There are zero SCF optimizations, no fitted constants and no
automatic corrected-protocol admission. Record the remaining archive residual,
XC and orbital-gradient response. This only measures the documented parameter
substitution; it cannot establish complete ORCA functional equivalence. Keep
unit-conversion, quadrature and RIJCOSX distinctions separate.

The first diagnostic job 22416 completed and saved its native PBE0 evaluation
(energy replay difference 9.09e-13 Eh), then failed serializing NumPy integer
metadata keys. A separate recovery converts those metadata keys and reuses the
exact saved energy, Fock and gradient. It verifies the failed submission/source,
terminal receipt and saved arrays before executing only the two remaining
declared variants. Same 8 CPUs / 32 GiB / 1 hour; zero SCF calls. Preserve the
failed log, first evaluation and original admission thresholds.

22419 completed the explicit-component control, then the custom functional name
failed PySCF's case-sensitive registry lookup after its dispersion parser
lowercased the name. Registering the same definition under a lowercase name
repairs this dispatch issue. The final recovery reuses both completed energy/
Fock evaluations and executes only the last declared parameter substitution.
Installed-backend metadata dispatch checks pass without molecular evaluations.
Both earlier failures and their actual attempted/completed call counts remain.

## Resource choice for any subsequently admitted constrained quartet

The archived mixed systems have 490 AOs. PySCF's `SCF._is_mem_enough` tests
nbf^4/1e6 MB plus current RSS against 95% of max_memory; its cached `s8`
two-electron integral tensor is about 57.6 GB at this size. The existing 32 GiB
allocations with max_memory=28000 MB cannot select that cached route. Aquo391
does fit, consistent with the observed runtime difference. This is a source-based
performance explanation, not a measured speedup claim.

If the scientific and implementation prerequisites pass, use 24 CPUs / 96 GiB
per independent constrained endpoint, max_memory=84000 MB, with a finite
4-hour limit. Preserve direct J/K, functional/grid and all accuracy gates. Record
whether the backend actually caches integrals and actual time/memory use; do not
change or cancel already running baselines. This resource declaration does not
admit the unresolved cross-code comparisons or submit the quartet.

## Declared native saved-orbital J/K diagnosis

Execute the finite two-call protocol in ORCA_FROZEN_SPEC.md: archived La-aquo
GBW and geometry, independent RIJCOSX replay and direct J/K one-Fock evaluations
in parallel, each 8 CPUs / 32 GiB / 1 hour, native cap 900s, standard_expand.
Prepared manifest: workspaces/lady_almo_implementation_20261023/orca_frozen_inputs_v1/MANIFEST.json.
No SCF optimization, changed physics or automatic baseline admission.

## Observed MPI startup failure and scheduling-only recovery

22432/22433 failed in orca_startup_mpi before a Fock evaluation: one Slurm
task with eight threads provides only one MPI slot. Preserve both failed
attempts. Retry the identical two native inputs with eight tasks, one CPU each,
32 GiB, one hour and the same native 900-second cap; standard_expand. A fresh
manifest and recovery record bind the originals. No oversubscription, changed
method, new geometry, extra successful energy evaluation or gate change.
