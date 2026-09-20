# Composite prepared-context scanner: working, 29 seconds/site

The opt-in operation completed all four fixed sites in **116.654seconds** from
prepared input, including process startup, model loads, native solvent runs,
validation and final collection. Allocation elapsed117seconds. The actual
scientific results reproduce the archived composite scores: maximum difference
6.403e−10 model-kcal/mol; all eight native MACE scalar energies are identical.
No new accuracy claim, preparation, threshold, DFT or default change.

| Case | MACE stage seconds | GFN2 stage seconds | Endpoint path seconds | Composite R, model-kcal/mol | Published decision |
|---|---:|---:|---:|---:|---|
| 1H4I | 10.327 | 20.668 | 30.998 | −405506.6479819513 | Ca-supported |
| Q9Z4J7 | 9.240 | 18.858 | 28.102 | −405464.18774828064 | Ca-supported |
| 1F6S | 9.907 | 19.044 | 28.954 | −405464.8337448340 | relative contrast only |
| 1GLG | 9.961 | 17.014 | 26.980 | −405489.1383748444 | relative contrast only |

Median endpoint path28.528seconds; full score-path mean29.164seconds/site.
Individual paths exclude final shared collection; the total includes it.
Alpha-lactalbumin−GGR remains **+24.3046300104**, one known biological direction.
The two PQQ timing controls are both Ca-family. This is execution repeatability
on four consumed cases, not another full25-reference accuracy test or prospective
validation. The earlier full-panel outcome remains in
[the solvent report](../compact_solvation_20260920/REPORT.md).

## Same energy, working interface

E(M) = E_OMOL,vac(M) + E_GFN2,ALPB(M) − E_GFN2,vac(M).
R = E(Ca)−E(La). eV and Hartree converted independently exactly once. Aquo S is
unavailable. Nothing adds GFN2 solvation to CPCM DFT. These are compact cluster
solvent corrections; no whole-protein dielectric boundary is claimed.

The native100M OMOL checkpoint and native GFN2 recipe are unchanged. The existing
MACE energy-only worker requests no forces; comparison to the archived
force-enabled endpoint verifies the scalar output, rather than assuming that
changing execution mode is harmless. All16 successful GFN2 endpoints terminate
normally with actual charge/state checks. Fixed bands remain exactly
Ca<=−405464.18774828, La>=−405459.1113819997. Q9Z4J7 remains on the Ca side of its
literal edge; no numerical tolerance pads that edge.

`prepare-pairs` accepts hash-pinned prepared Ca/La XYZ, formal charge,
multiplicity, model/software and context-preparation provenance. Archived MACE
energies are optional, demonstrated by dry-running the real1H4I inputs with
those fields removed. The operation uses existing worker/executor/parser code;
it does not launch unmanifested tasks or calibrate against the incoming cases.
PQQ decisions require **the canonical PQQ complete-context protocol plus PQQ
scope**. A generic fragment cannot acquire PQQ bands merely by changing its
label. Other chemistry gets raw R. Supported preparation IDs are explicit;
this is neither arbitrary-XYZ acceptance nor a raw-fold scanner.

## Real cost and the retained failure

Successful job**1203319**: one H200,32CPUs,200000MiB. 117allocated wall/GPU-seconds,
**3744allocated core-seconds**. CPU solvent work uses four concurrent8-rank
endpoints. MACE is one32CPU/GPU step; GPU idle during CPU stages is charged.

Initial job**1203299** used the same resources but a1×32CPU Slurm task layout.
All8MACE scalars succeeded; all16ORCA launches failed in startup because
OpenMPI saw only one task slot. No SCF energy was accepted. Recovery changed
only task layout to32×1CPU, with an explicit1×32CPU MACE step. No
oversubscription flag, changed scientific parameters or shared permissions.
Failures/receipts remain immutable; see[recovery](RECOVERY.md).

Total development here: **5120allocated core-seconds/160allocated GPU-seconds**,
16MACE calls and32GFN2 launches, including16pre-SCF failures. No DFT. The complete
recovery repeated8MACE calls to measure genuine integrated latency; it did not
claim cached work as freshly timed. Counts and allocations include both jobs.

Peak observed MACE VRAM2,785,492,480bytes (2.594GiB), worker RSS1,923,424KiB
(1.834GiB). These are process observations, not an aggregate whole-job host
memory peak; that value remains unavailable. Reserved host memory was195.313GiB.
The eight native inference operations themselves took3.120seconds total;
process/model startup and solvent execution dominate this initial thin path.
No lower-memory allocation or alternate GPU was tested here. No matched DFT
speedup is inferred from different hardware or archived timings.

## Qualification and tests

Execution repeatability: **pass,4/4 scores and8/8 native scalars**. Independent
native solver precision: parent reports6/6contexts pass tightening, maximum
correction shift0.025545kcal/mol, below the frozen0.2 limit. Independent ordinary
SCF: no usable complete pair; explicit-GBW restart converged1/8endpoints (that
one agrees within0.000505kcal/mol). These are distinct findings; see the
[numerical qualification report](../compact_qualification_20260920/REPORT.md).
The original frozen collection retains its conservative contemporaneous status.

**Eight real-fixture tests pass.** They cover source/state corruption, exact band
edges, preserving the GPU venv invocation, explicit input without old energies,
actual failed MPI components staying null, actual complete unit/sign replay, and
rejecting generic preparation for PQQ bands. No fabricated successful scientific
output or skipped executable substitute. A post-run collector review tightened
PQQ protocol gating and consumed-evidence metadata, then replayed existing
outputs with unchanged scores; no chemistry reran for this review.

## Recommendation

Keep this as a **working opt-in prepared-context research scorer** alongside the
promoted baseline. About half a minute/site is a practical measured score-path
cost for these compact cases. It does not include folding/context generation,
prove throughput at metagenomic scale, or establish broader predictive validity.
No further science was needed for this timing phase. Use the explicit pair
operation for compatible future preparations, with separate evidence roles and
frozen calibration, before considering any default promotion.

Sources: [RESULT.json](RESULT.json), [allocation receipts](SACCT.txt),
[tests](TESTS.txt), [commands](COMMANDS.md). The two original job collections,
source/model pins and every endpoint receipt remain under
`workspaces/compact_scanner_20260920/`.
