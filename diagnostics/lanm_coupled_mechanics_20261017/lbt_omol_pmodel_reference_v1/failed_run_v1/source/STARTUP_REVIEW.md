# Independent startup observation:16611

Observed21October2026 at03:36:00UTC. **Startup allocation and input identity match the released PModel variant. No parsing error is visible. Actual PModel initialization and beta rotation have not yet been reached in the inspected output.** This is a bounded text/hash/metadata review, not a terminal result or numerical replay.

## Exact bindings

| Artifact | SHA256 |
| --- | --- |
| CONFIG_RUN.json | `b654d846a9e000cb7adb0c66718b38df7ab64451f5b242283ff843c176df5a0d` |
| ROOT_EXECUTION_DECISION.json | `6df38e02f9c41d2e088e065d7916dee848e42e3d47937bd531e80c9f181ca982` |
| execution/START.json | `0cae9f1199296a18e9774848464eecfefc72762df77f88d73a3c841d67bd66b5` |
| execution/MANIFEST.json | `d1416ce14b472bd5c92668ab47892cd297b2471f9a5de328083de4d7074edafd` |
| execution/SCONTROL_JOB.txt | `ddc6e5a8bac9cc3f5fe0ea2ea0b61cfc89ad1fff8d8eb486e7ad11dcbe76c7b2` |
| central_La/endpoint.runtime.inp | `b2893d69ad00d4d4b8998628209c6c3bd99fbd2474a3d2d70a0ddc0b59c7a4b9` |
| central_La/endpoint.runtime.json | `78c29fed0cefd5dfdaaf0927bb3b74e0dbafacd9324f6d0c20fdabb9a3e42652` |

Selected source/configuration/decision/manifest/evidence/input/XYZ/basis hashes match. The task signature exactly matches actual preparation16609. Input evidence remains `48ee5022…`, template `4359537c…`, XYZ `55e07d81…`, external basis `7389c56f…`; source remains `27ac4e80…`. The actual runtime input is byte-identical to16602's runtime input after the sole PAtom→PModel substitution, including unchanged MaxCore and PAL settings.

## Observed execution

The saved allocation and runtime sidecar report one node/worker,24 CPUs and local MPI slots,96000MiB,25% reserve and MaxCore3100MB per rank on `node-64-768g-13`. No prior orbital file is supplied; the input retains NoAutostart, UKS, beta rotation `{514,515,20,1,1}` and all prepared scientific controls.

The live output was103743bytes/1929lines, modified03:35:59.618492UTC. This changing file is not given a final immutable hash here. Native line54 identifies ORCA6.1.1;205/208 identify the external orbital basis and def2/J. Lines223–227 again switch EnGrad to SCNL and AutoTRAH off. Lines1835–1837 identify external-file La ECP46,lmax3 at atom255. The latest memory check says four symmetric matrices of dimension6239 fit. No explicit parsing/native failure was observed.

Completion watcher66459 and health watcher66460 have matching RUNNING/watching receipts; health was armed03:35:52UTC. These are monitoring registrations, not convergence evidence. RESULT.json and a native terminal receipt were absent.

The initial-guess, molecular state and SCF stages remain pending. In particular, the PModel keyword is verified in the runtime input, but successful La atomic-density generation, MODEL POTENTIAL execution and beta rotation cannot yet be claimed. The early memory check does not establish peak resource use. No intervention is indicated by this snapshot; retain unchanged execution and the existing bounded review process. Native method/state/gradient admission, actual internal atomic work and cost must be reviewed later. No additional endpoint or displaced call follows from startup.
