# Independent startup review

**Observed startup inputs and allocations match the released protocol. Actual molecular SCF/PModel headers and the plus health-arming confirmation were pending at the bounded inspection.** No startup mismatch or evidence requiring intervention was found.

Observation: 20 October 2026, approximately 20:47–20:48 UTC; final immutable-input metadata check at 20:48:11 UTC. This was one output/receipt inspection, without polling, numerical replay or scheduler changes. Mutable output and watcher receipts are described as observations, not final immutable artifacts.

## Allocation and selected tasks

| Task | Job | Host | MPI ranks | Allocated memory | MaxCore per rank |
| --- | --- | --- | ---: | ---: | ---: |
| minus | 15847 | node-64-768g-16 | 12 | 48000 MiB | 3100 MB |
| plus | 15848 | node-64-768g-16 | 12 | 48000 MiB | 3100 MB |

Saved scheduler receipts show both jobs RUNNING, partition `standard`, QoS `standard_expand`, 12 tasks with one CPU/task, and no GPU. Launchers contain no memory or time override; resolved scheduler limits are 4000M/CPU and five days. Runtime sidecars retain one worker, 12 MPI slots, `SLURM_MEM_PER_NODE` memory authority and the original 25% reserve. These are allocation/settings records, not measurements of rank utilization.

Both START records select only their requested singleton. CONFIG_RUN `dd4875fabd6d17a12d8c73ad7254068ba83a7070b33d079c62b3b7c2dcd74ff7` differs from the prepared configuration only in the preparation-result and execution-decision pins. Release `62cd2ead92154df6bfcd853be08d81ccdd5da929039c4872af459c7658aef35c` binds the reviewed preparation `779412e9`, the unchanged source `cd4396c0` and the two-call ceiling.

## Input and native observations

Each singleton's input, XYZ, physical metadata and position hashes exactly match its prepared task signature. Both tasks retain all 256 atoms, charge −2, physical/effective La singlet, ECP46, 1030 explicit electrons and target index 255. The runtime decks preserve PBE0-D4/def2-TZVP/def2-J/RIJCOSX/DefGrid3/VeryTightSCF, `NoAutostart`, `Guess PModel`, `%pal nprocs 12` and `%maxcore 3100`. EnGrad is absent and `engrad_path` is null; no seeded-wavefunction input is introduced.

At the output inspection, minus/plus files were 103687/103685 bytes and ended at the printed La ECP block replacing 46 core electrons. The input declarations and printed ECP agree. Molecular charge/NEL/multiplicity, actual PModel execution and SCF iteration evidence were not yet printed in the inspected tail. Their eventual native admission remains pending; declared settings are not substituted for observed convergence.

## Watcher observations

Both completion receipts were `watching`, with RUNNING scheduler records checked at 20:47:31 UTC and no query error. Minus HEALTH_ARMED and health receipt existed; its watcher was `watching`, bound to the exact minus manifest, with no alert and automatic cancellation disabled. Plus HEALTH_ARM_SUBMISSION existed with alert-only armer PID 2466404, but HEALTH_ARMED and health_receipt were not yet present at this inspection. The plus armer was therefore submitted, with successful arming confirmation pending. No additional wait or poll was performed.

## Verified startup pins

| Artifact | Minus SHA256 | Plus SHA256 |
| --- | --- | --- |
| START.json | `1fd2e19c85b0004a8500dc22b5f510836d032c8187a2375925aa899a4615f5c0` | `01cd4eb3419acccd8b33cffa38ea08b204457496df57272f5c1ede8d30cce2e3` |
| MANIFEST.json | `a10fe4f3b7d29443407fb46f6a9567c2a2ba6f63c466c5759a87e9c6c3c868ce` | `0dcd90f6bacb0689d67bb4eef3759f52a985aa53dc22215559e28ee5704c7029` |
| Runtime sidecar | `58aab5d10eb8bcf55434464ac62f9d3e6e55e03da35517947217ffc7e1c3b4cd` | `d3e2f7e759aa0d51fddf305234011af4684c46b4b1338c536178d439b0825840` |
| SUBMISSION.json | `58d50155fe1f694123d054ca658c7a39d59cf065c2ed3fe24402a0e705cc7b2a` | `f219dc4c5a5ca43bc70e79bbafeaed215c734333b09bfd73ba4fcf1fc8853609` |
| HEALTH_ARM_SUBMISSION.json | `cc1b41e9824dff8a5292ea2a7ca65ec6e38109ab8be82088bca9b5783ff23b76` | `95c551e29e9e0d3097344406910c34386106b6d403c20bb5abdfee9c308bc369` |

The shared runtime input is `17749f48c64fdd266a142c925aca71ec4f7cc3337594c3aeb6ac12047f559f23`; its template is the reviewed `e99306dc` singlepoint deck. Native outcomes, final costs and the directional-consistency verdict remain pending. This review admits no extra call, retry, state selection or scientific qualification.
