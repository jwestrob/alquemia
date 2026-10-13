# Fixed-origin PedE DFT reference: prepared; first scout running

The eight-cell matrix is frozen at
`workspaces/pede_solvent_response_20261013/reference_v1/manifest.json`, SHA-256
`1d5bff48eec698c83f2e5a643d13061b33d1b21752e78e7193592e81a77c669f`.
Root submitted **only sample1 La/CPCM-water**, job3354, with its completion and
health wakes. Its actual Slurm preflight passed; the molecular single point is
still running. Vacuum and the other six cells remain unexecuted. This is a new
diagnostic reference, with no classifier fit or production change.

## Frozen question and state

Compare the solvent response at the exact saved **origin** geometries of PedE
two-Ca samples1 and3. The matrix is2samples × Ca/La × vacuum/CPCM(Water).
Every task receives a byte copy of its original192-atom XYZ. No relaxation,
coordinate alignment, new cap, proton, water, spectator, gradient, or electronic
restart is introduced. The full-source second Ca remains outside the already
prepared local context; this calculation does not test its missing interaction.

The two sources share the same state signature
`ae9b123d0a1c907bd91b6fef222cbf96339fd58639084a4b26b599a9dad2dd57`.
Actual Slurm preflight verified192 atoms,14 caps, zero explicit waters and exact
Ca/La nonmetal coordinates. Composition is M1 C66 H83 S2 O24 N16. Ca charge−3,
singlet,838 explicit electrons; La charge−2, singlet,874 total/828 explicit
electrons with the46-electron ECP. All mapped nonmetals1–191 are retained.
The manifest pins original preparation/admission/state and all four XYZ sources.

The intended comparison is
`S_z(sample)=E_CPCM(z,sample)−E_vacuum(z,sample)`, then
`D(sample)=S_Ca(sample)−S_La(sample)` and `D(1)−D(3)`.
CPCM's printed dielectric energy is a separate diagnostic already included in
the total. It is not an extra additive correction or a replacement vacuum energy.
No favorable subset, cutoff, calibration band, or class label is introduced.

## Method, executable and solver semantics

Templates use `! r2SCAN-3c NoAutostart DefGrid3 TightSCF`, with `CPCM(Water)`
only in the water cells. The TightSCF policy follows the existing native
`scripts/mace_omol_vacuum.py` reference, omitting its gradient request. Historical
NormalSCF DFT files remain unchanged and are not silently pooled with this policy.
No native-GFN2 JSON or MACE weight enters these native DFT inputs.

ORCA6.1.1 is pinned at
`/groups/banfield/users/jwestrob/bin/ORCA/orca_6_1_1_linux_x86-64_shared_openmpi418_nodmrg/orca`,
SHA `38b5f057452fef275c0a1b98d270ad03d0411dab70453820d1c678bd732f6c83`.
The manifest also pins the serial/MPI startup, guess, leanscf, prop and propint
helpers plus bundled libraries. Actual compute-node MPI, loaded libraries,
environment and affinity are recorded by preflight, rather than assuming that
login-host system libraries equal compute-node libraries.

Actual historical ORCA outputs establish native `def2-mTZVPP`,
`def2-mTZVPP/J` and La Def2-ECP replacing46 core electrons. Their CHNOS/Ca/La
coverage is directly observed. Current3354 output confirms192-atom La NEL828,
native basis, singlet−2, MaxIter125, TolE1e−8, TolG1e−5 and DIIS threshold5e−7.

The [ORCA6.1 SCF manual](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html)
documents TightSCF's energy-based mode2 and TRAH's gradient-norm stopping
criterion. TightSCF alone therefore does not establish small final density
changes. The collector retains all printed density/DIIS/orbital summaries and
the actual last TRAH macro residual; it checks the latter against printed TolG.
It does not reject final TRAH from pre-transition density/DIIS summaries.

The prior `metal_environment_response_20260926/SCF_CHECK_REPORT.md`,
`scf_check/FAILURE_DIAGNOSIS.md` and `GRID_CHECK_REPORT.md` showed that strict
stopping and denser grids did not repair a distinct rigid-rotation discrepancy.
Those failed qualification records remain failed. This matrix does not repeat
that sweep or inherit force/density-response qualification from it.

## Execution, evidence and first-scout cost

`REFERENCE.py` calls a byte-pinned copy of `run_orca_task_manifest.py` unchanged.
`metal_environment_runtime.py` is copied as its runtime renderer, with the
unchanged base renderer beside it. One task uses12 MPI ranks, one thread/rank;
the scientific template is untouched. Actual3354 granted48000MiB, yielding
**MaxCore3100 decimalMB/rank** after the existing25% reserve. There are no memory
or time submission overrides. The runner's native MPI binding fix is preserved.
A workspace lock serializes execution; do not launch parallel subsets against
this manifest. Completed matching tasks are reused; failed/partial tasks are
preserved and never silently restarted.

Actual3354 preflight also passed three archived parser checks:139-atom Mulliken,
139-atom Loewdin, and the known final TRAH macro residual5.528534e−7. Source review
by the independent agent found no launch blocker. Local checks were AST and
shell syntax only; no molecular import, geometry analysis or scientific test ran
on login. Root owns every submission. Unrelated2669/2832/2833 remained untouched.

The nearest explicitly inspected native runtime precedents are139-atom PedE
CPCM Ca/La endpoints at16 MPI:25m22.378s/27m8.843s,33/37 SCF cycles. The largest
reported module memories were503.5/470.4MB; these are **not aggregate process
RSS**. Two related155/161-atom systems took about29–34min per endpoint at16 MPI.
These used historical NormalSCF and different coordinates. No exact192-atom
TightSCF runtime or speedup is claimed;3354 measures it directly.

The bounded metadata cache search inspected four known native-DFT manifests,
listed in the new manifest. None mentioned the four exact XYZ hashes. This is
not an exhaustive archive search, and old139-atom PedE energies are incompatible
with these new192-atom origins.

## Completion review

Read `attempts/3354_execute/EXECUTION.json`, its collection, the native
`endpoint.out.execution.json`, runtime memory sidecar and raw output. Check
normal termination/SCF convergence, ORCA6.1.1, NEL828, native basis/La ECP,
complete192-atom Mulliken/Loewdin populations and retained `.gbw`. Inspect the
actual ending solver, final residuals and printed numerical settings before
expanding. Normal termination alone is not physical density-response validation.

Collections retain atom/group populations (metal, PQQ, caps, protein), dipole,
CPCM dielectric term, final solver evidence and hashes of wavefunction/property
artifacts. Grouping uses the full192-atom preparation map, not the189-atom
kinematics subset. Population charges are method-dependent diagnostics. A
subsequent matched restart, if specifically justified, needs a separate declared
input; `NoAutostart` prevents one from happening silently here.

## Exact continuation interface

After root reviews the water scout, the next declared cell is sample1 La vacuum:

```bash
REFERENCE_WORK=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pede_solvent_response_20261013/reference_v1
sbatch "$REFERENCE_WORK/implementation/RUN.sbatch" \
  "$REFERENCE_WORK/manifest.json" \
  1d5bff48eec698c83f2e5a643d13061b33d1b21752e78e7193592e81a77c669f \
  execute sample1_La_vacuum
```

Use a healthy node/output destination and completion wake as root's scheduler
wrapper requires. This command has not been submitted. Subsequent explicit IDs
are `sample1_Ca_water sample1_Ca_vacuum sample3_Ca_water sample3_Ca_vacuum
sample3_La_water sample3_La_vacuum`; none launches automatically. `preflight`
performs checks without molecular calls; `collect` parses existing artifacts.
Both also require an actual Slurm allocation. Each job writes a new attempt
receipt; no final collection is presented as complete while cells are missing.

## Minimal optional parallel layout after scout review

The active snapshot and its global lock stay unchanged. If measured cost warrants
parallel execution, run the existing metadata-only `PREPARE.py` into three new
versioned workspace directories. Each contains the same eight scientific inputs;
only its manifest/snapshot/output paths differ. Use explicit disjoint task lists:

| Workspace | Newly executed cells | Reuse |
|---|---|---|
| original `reference_v1` | sample1 La vacuum | its actual sample1 La water scout |
| new `reference_sample1_Ca_v1` | sample1 Ca water, vacuum | none |
| new `reference_sample3_La_v1` | sample3 La water, vacuum | none |
| new `reference_sample3_Ca_v1` | sample3 Ca water, vacuum | none |

At most four12-MPI jobs use48 CPUs; root must check the live aggregate60-user
limit and other owned jobs before submission. Each manifest retains its own
lock, native outputs, wavefunctions and per-job receipts. Pin a dispatch record
joining original matrix task ID to exactly one workspace/attempt, and require
template/XYZ/protocol/source-state hashes to agree with the original matrix
before scoring. Do not copy the scout output into another manifest or rerun it.
Root can join these eight actual records for comparison; an individual shard's
partial collection does not imply missing matrix science. This is a preparation
plan only: no shard directory, new runner, new submission or snapshot edit has
been made for it.
