# Repaired Mex pair is complete; Hans/Mex preference remains unavailable

Worker1220336 and collector1220337 completed both repaired Mex8FNS endpoints normally. Native energies, atomic gradients and point-charge gradients passed the existing collector. La energy−5636.761993646614Hartree (33SCF cycles); Dy−5642.249821148190Hartree (43cycles). State/source/water/spectator conventions remain those frozen in the original comparison plan.

The electronic Dy−La contrast is−3443.663748916275kcal/mol; adding the declared matched classical difference gives−3443.597047554222. These cross-element totals contain a large element-dependent offset and are NOT an affinity or a Dy-preference classification. The balanced Hans-minus-Mex contrast cancels the common offset but still requires the missing accepted Hans pairs. All Hans contrasts remain null; no unconverged trial energy is used.

Both physical gradients assembled successfully with their exact repaired classical partners, using the existing source/cap/point-charge force mapping. Assemblies and source-table hashes are in MEX_RESULT.json. Mapping consistency is not complete hybrid derivative or biological validation. No new model, displaced structure or quantum call was introduced during collection.

The worker cost6002s×224allocatedCPUs; collector2s×1CPU. Total1,344,450allocatedCPU-s,zeroGPU. La finished earlier; Dy completed in the same original allocation. These are actual allocated resources, not measured utilization or a production speed claim. The successful partial alert was reviewed; retaining the live Dy calculation avoided a disruptive restart.

The recovery comparison now contains3/6 accepted cells (Hans8DQ2La plus MexLa/Dy). Primary failedFNR and separate recovery identities remain intact. Full electronic/finite tables: workspaces/metal_environment_response_20260926/repaired_exchange_Mex_complete_v1.json. Original full-panel automatic postprocessing1220344 still waits for the existing Hans collectors; do not duplicate it. No PQQ/default change or reserved-label access.
