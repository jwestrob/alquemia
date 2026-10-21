# One resource-only recovery of the central PModel reference

21 October 2026 UTC. Root authorized one fresh resource-only recovery after job16611 exhausted its 96,000-MiB allocation during main LEANSCF. Its PModel La atomic initialization and beta rotation had completed. Preserve failed16602 and16611 source, inputs, outputs and costs.

## Exact scientific scope

One requested whole256 7CCN La step725 energy/gradient endpoint. The scientific template must have SHA256 `4359537c4a569f2a68793f60a509bda509c854a15a7aaf1060790a9407086896`, physical.xyz `55e07d8146adc13c09203243849b65e80a6afef89b181cfe4c22b265a994b72d`, and external basis `7389c56f0cff9e76098b0757a34d1abdd79c48e9249597cea074a0106a3462a6`, exactly as16611.

Keep PModel/NoAutostart, charge−2 UKS singlet, 1030 electrons, external-file La ECP46, beta rotation `{514,515,20,1,1}`, wB97M-V EnGrad and native SCNL/libXC treatment, def2-TZVPD, def2/J, RIJCOSX, DEFGRID3, Tight/300, DIIS/NOSOSCF, THRESH1e−12, TCUT1e−13 and all other scientific input bytes. No seed/GBW or integral restart. This remains a PModel derived-recipe variant on6.1.1, not exact historical OMol25 reproduction. No molecular geometry change, optimizer, other metal, MACE, ALPB, GFN, MM or displaced DFT call.

## Sole resource change and actual policy

Failed16611 used24 MPI ranks,24 allocated CPUs and native MaxCore3100MB/rank. It printed an estimated4650.8MB requirement and Slurm recorded OUT_OF_MEMORY. Use **12 MPI ranks, two allocated CPUs per task,24 allocated CPUs total**, on one standard node under standard_expand. All numerical thread variables remain1. No memory/time override; current default4000MiB/allocatedCPU supplies96,000MiB. The unchanged allocation renderer reserves25%, converts MiB to decimal MB, divides by12 ranks and rounds down to100MB, producing **MaxCore6200MB/rank**.

The current standard_expand QOS caps allocated CPUs/job at24; the rejected24×2/48CPU suggestion is not used. Extra allocated CPU capacity supports the memory entitlement; it is not a claim of24 active MPI processes or measured utilization. Record actual scheduler allocation, wall time, allocated CPU cost and native memory evidence. Fewer MPI ranks may alter floating-point reductions or SCF trajectory; no unchanged-result or memory-sufficiency guarantee is made.

Source basis: unchanged `render_orca_runtime_input.py::allocation_policy` distinguishes allocated CPUs from local MPI slots and does not replace slots with CPUs-per-task. The unchanged runner accepts one worker×12 ranks within24 CPUs and forces one numerical thread. The adapter requires observed NumNodes1/NumTasks12/CPUs-per-task2, allocated24 CPUs, local12 MPI slots,96,000MiB and the derived6200MB setting. No policy or numerical kernel is modified.

## Finite execution and admission

`TARGET.py --stage prepare|run --config CONFIG --output FRESH_DIRECTORY`

One CPU zero-call preparation admits pinned16611 result/config/source/input/XYZ and the original16602 lineage, exact saved physical/electronic state, all method assets and byte equality. Root must bind that actual result and a decision with status `approved_single_memory_recovery_central`, permitted task `central_La`, mpi_ranks12, allocated_cpus24, cpus_per_task2 and maximum_new_DFT_calls1 before the single fresh native endpoint.

At most one whole-molecule runner invocation; native PModel atomic-density subcalculations remain allowed and included in its cost. Preserve their surviving `*_atom*` files through the inherited artifact manifest. One requested endpoint is not one total internal solver invocation. No automatic retry, alternate settings, intermediate scout, SP release or additional preflight beyond the declared preparation.

All original scientific admission and comparison arithmetic remain: native termination/convergence/state/ECP and finite256 gradient,1e−8Eh energy agreement,1e−6Å gradient-coordinate tolerance, frozen full-vector direction and Cartesian-component RMS. Success remains `native_core_admitted_pending_method_review`. Native method/state/rotation/grid/density/S²/gap evidence and failures remain explicit. No force-accuracy, stable-state, equilibrium, absolute cross-method deformation work or metal-preference claim.

Only source/text/hash/AST work is performed on login. Numerical imports, saved replay and native work are Slurm-only. Root owns launchers, actual execution decision, submission and monitoring. Future preparation/result pins are null until actual evidence exists.
