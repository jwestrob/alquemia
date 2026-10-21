# Independent actual preparation review:16609

21 October 2026 UTC. **PASS: actual Slurm16609 is `prepared`, with zero engine calls. No preparation blocker to the separately released single24-MPI PModel endpoint.** This review inspected saved fields, hashes and literal text/byte equality only; it performed no numerical or scientific replay on login.

## Exact pins

| Artifact | SHA256 |
| --- | --- |
| preparation/execution/RESULT.json | `95e1852b789ed69976ee85ab5a616a59b487c058107cfb279adcabe4ede586e3` |
| INPUT_EVIDENCE.json | `48ee50224d32618cfece6b6ad8534461a4d370c86f89da8324ae7da4c3d81078` |
| MANIFEST.json | `d0ed827a3bf924b5f3b1d66d1df96a251eb5cc6abe7e2fc6b71c465af91fe1e8` |
| TARGET.py | `27ac4e80b5fb53c049adc6f67822f1629cd3bd7cf372edf1d4708cedcc829856` |
| CONFIG_PREPARE.json | `a1007083e3d3139848b18bd3157d2cf0746f52e429c42f44dc34872fcfc159ad` |
| central_La/endpoint.inp | `4359537c4a569f2a68793f60a509bda509c854a15a7aaf1060790a9407086896` |
| central_La/physical.xyz | `55e07d8146adc13c09203243849b65e80a6afef89b181cfe4c22b265a994b72d` |
| central_La/def2-tzvpd.bas | `7389c56f0cff9e76098b0757a34d1abdd79c48e9249597cea074a0106a3462a6` |

Selected pins and all four frozen failed-parent bindings match. The new deck is exactly16602's template after its sole `PAtom` token becomes `PModel`. XYZ and external basis are byte-identical to16602. Configuration/result/manifest scientific bindings agree; the only manifest task equals its TASK.json. Physical-parent, positions, all physical IDs, atom count, target index, charge and multiplicities exactly match the failed parent.

The saved on-Slurm preparation admitted the unchanged step725 whole256 source, q−2/UKS singlet,1030 electrons,515/515 occupations, La ECP46 at index255,6239 orbital functions and beta rotation `{514,515,20,1,1}`. Input evidence records PModel and actual expected SCNL treatment with external-file ECP metadata. The manifest declares only `central_La`, one worker and24MPI; no additional task is present.

All runner/DFT/MACE/ALPB/GFN/MM, displaced, optimizer and other-metal counts are zero. The task directory contains only TASK.json, input, XYZ and basis. Native `unavailable/not_attempted`, comparison unavailable and null atomic-guess artifacts are correct for preparation. Application elapsed time is5.918838836019859s; scheduler cost remains a separate accounting record.

This is input preparation, not native success. The run still requires root's actual-result release, fresh files, one requested molecular endpoint and retained internal atomic-guess accounting. Actual PModel/rotation/SCNL/method evidence, native convergence, full gradient and measured cost must be reviewed. The changed initial guess can change the electronic solution; historical ORCA6.0.0 reproduction, physical-state and force qualification remain unproven. Failed16602 remains immutable; no automatic retry, displaced pair or optimization follows.
