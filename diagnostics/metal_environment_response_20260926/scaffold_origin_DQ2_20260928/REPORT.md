# Repaired Hans8DQ2 classical force bridge complete

Both La/Dy repaired source-A origins now have full component and total Cartesian gradients. The classical Dy-minus-La contrast is **+0.058069241204521 kcal/mol**, entirely the target cross-LJ term (+0.058069241204770 before total summation roundoff). All common bonded/MM-LJ/MM-Coulomb terms and gradients cancel exactly between target metals.

Saved origin projections were insufficient to reconstruct full gradients. Two authorized new classical configurations/eight component queries filled that specific gap. Actual source atoms, all four caps' exclusion policy, shifted field charges, La/La/Na spectators, native XML, bonded terms and exceptions match the earlier repaired coupled-scaffold ledger exactly. No original-H arrays were reused. Compared with archived repaired origin energies, maximum discrepancy is4.09e−12 kcal/mol, below the predeclared1e−9 tolerance; cross-LJ energies reproduce exactly. No quantum calculations or optimization occurred.

Full bridge paths:
`workspaces/metal_environment_response_20260926/scaffold_origin_DQ2_v1/Hans8DQ2/{La,Dy}_A.json`.
Each records1891 physical source IDs, coordinates, gradient_kcal_mol_A, component_gradients_kcal_mol_A and component energies, input/ledger pins. Gradients are negative OpenMM forces, in kcal/mol/angstrom. Actual target identity is `A/203//LA` for both endpoint substitutions; no synthetic caps are present in this physical array. Same interface as the Hans8FNR/Mex8FNS transfer.

Slurm1220340 COMPLETED:3elapsed seconds on1allocatedCPU = **3CPU-s**, mem=0, zeroGPU. Script wall2.130062s, process peakRSS167088KiB. All real-artifact identity/energy/gradient tests pass. Scheduler MaxRSS0 is not used as a memory measurement. This brings the newly saved six-origin classical bridges to8CPU-s total (otherfourcost5CPU-s).

This completes mechanical records for later matched electronic integration. It is not new force numerical qualification, bulk solvation, affinity prediction or parameter validation. Old coarse/halfstep failures and quarterstep numerical pass are preserved. No old output, production setting or molecular state changed.

Recheck without molecular evaluation:
```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/metal_environment_response_20260926/scaffold_origin_DQ2_20260928/TESTS.py
```
The batch file records the completed execution, not a request to rerun it. Root owns integration and further work.
