# Four-call isolated SCC damping diagnostic

## Actual failure and scope

Hybrid pool job12992 stopped on new composite59: 7CCN/baseline/La,
`cartesian_03`, isolated172-atom GFN2. Its embedded companion was accepted.
The original failed POINT, native outputs, full104-cell denominator and all
completed components remain immutable. This new solver diagnostic has its own
protocol and receipts; it cannot relabel the old failure or resume the pool.
No optimizer, MACE, MM, ALPB or DFT calculation is included.

The failed isolated calculation ran26.983669seconds. Its printed SCC trajectory
reaches RMSdq0.00557 at iteration30, then oscillates. Iterations81–87 alternate
large electronic-energy jumps approaching10Eh; the final500th iteration has
dE−3.10611Eh and RMSdq0.128. These are unconverged SCC electronic energies, not
accepted total energies. Printed thresholds remain1e−10Eh and1e−8e. The same
geometry's embedded calculation converges in60 iterations, final RMSdq1.36e−9
and dE6.82121e−13Eh. Immediately preceding point058 isolated converges in48,
final RMSdq2.53e−9 and dE7.95808e−13Eh. This supports testing convergence-path
controls rather than increasing iteration count; it does not identify a unique
physical electronic state or establish that the geometry is wrong.

## Supported native control

The executed binary is installed xTB6.7.1, conda buildh8876d29_4,
SHA481a1e7d8151fc08e0341b8e06ab60ac09cf1c2628738eeb83346271274bf8f4.
The package manifest explicitly hashes both this binary and its installed
`share/man/man7/xcontrol.7`; that manual documents `broydamp` in the `$scc`
group. Archived6.7.1 source `set_module.f90:set_scc` reads it into
`set%broydamp`; `scc_core.f90` passes it to Broyden mixing, separate from the
Hamiltonian and convergence tests. The packaged patches touch DIPRO/tblite,
optimizer/build/version files, not this legacy SCC parser/mixer path. All source,
manual, package and recipe evidence is pinned in CONFIG.json.

The source's default is0.4; its comment names0.05 for difficult cases. This is
support for the solver control, not evidence that it will converge this case.
The source also contains native gap-dependent damping behavior. We do not change
it or claim the requested coefficient is held throughout every iteration. Every
new run must print the requested **initial** Broyden damping, establishing actual
parser uptake. Unknown/ignored controls fail the diagnostic.

The only added native input is:

```text
$scc
  broydamp=0.05
$end
```

The second target call changes only that value to0.10. Existing command flags are
unchanged, with `--input control.inp` appended. No charge-guess override is added.
No restart file is seeded; `--norestart` stays present.

## Frozen four-call matrix

| Case | Exact archived isolated geometry | Broyden request | Comparison |
| --- | --- | ---: | --- |
| failed_005 | failedpoint059 | 0.05 | target0.10 |
| failed_010 | failedpoint059 | 0.10 | target0.05 |
| preceding_005 | acceptedpoint058 | 0.05 | exact accepted default output |
| origin_005 | accepted7CCN/baseline/La original, origins point004 | 0.05 | exact accepted default output |

All calls use unchanged172 atoms, original coordinate bytes, charge−2,
538 valence electrons, UHF0/physicalLa multiplicity1, 300K, accuracy0.0001,
maximum500iterations, no external field, no solvent and the full qualified
published standalone parameter filef6f05c6c. Corrected La/Dy profile lineage635e111b
remains explicit. Raw native analytic gradients and populations are required.

## Admission and agreement

Reuse exact existing GFN_ADMISSION.py, gradient reader and backend.call. Require
normal native termination, original SCC/temperature/accuracy/state/charge checks,
finite172-atom analytic gradient, original energy/coordinate/unit consistency,
no external point-charge gradient, and actual printed requested damping.

All four native calls and all three comparisons must pass:

- Absolute energy difference <=1e−7Eh.
- Maximum Cartesian gradient-component difference <=0.005kcal/mol/A across all
  172 native atoms and across all256 physical atoms after the exact archived
  cap scatter. No coordinate or map is regenerated for this comparison.
- Report complete charge vectors and their maximum difference descriptively;
  original charge-closure/population admission remains unchanged.

Two dampings agreeing is a finite solver-branch consistency check, not a proof
of the stable ground state, force accuracy, unique SCF solution or general
recovery policy. An energy-lowering sign is not an acceptance criterion.

## Execution and accounting

`SCOUT.py --stage preflight|run --config CONFIG.json --output FRESH`.
Both stages require Slurm. Preflight replays the two accepted raw controls,
verifies failed/native/source/package lineage, and writes all controls and exact
case bindings with zero new calls. `run` performs the same preflight first and
then the four declared calls sequentially with one native thread each. It needs
no GPU and never constructs HybridEvaluator or OpenMM. Root owns submission.
Each admitted call also emits flat `NATIVE_RESULT.json` using the exact admission
payload, native task/execution and explicit changed-solver provenance. The
diagnostic `RESULT.json` remains separate. This allows later declared component
reuse without recomputation or silently rewriting the original protocol.

An individual SCC failure is retained while the other already-declared controls
may still run. There is no automatic retry or fifth call. A surviving native
process group stops further scheduling. All four slots exist before preflight;
failures preserve unavailable siblings. Each new call has explicit TASK/control,
START/EXECUTION and raw output pins, and separate solver-diagnostic provenance.
Actual attempts, native wall/CPU time and allocation are reported. Original
failed and accepted artifacts are never overwritten or removed.
