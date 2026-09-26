# Whole-protein LanM native feasibility

Job 1216564: FAILED; 2/4 native origin cells passed collection.
Node node-112-1500g-1; 112 allocated CPUs; 21141 elapsed seconds.

These are the same Hans8DQ2 EF12 La/Dy vacuum/ALPB endpoints. No new MACE or DFT ran.
Both earlier relaxed geometries remain rejected for covalent distortion. This result tests native whole-protein execution; it is not a validated within-series preference or accommodated score.

Resources: {"concurrent_tasks": 4, "maxcore_mb_per_rank": 10800, "mpi_ranks": 28, "scheduler_memory_mib": 1546754}
- Hans_8DQ2__EF12__origin__La__vacuum: unavailable; energy Eh=None; reason=no successful native execution receipt
  [file orca_leanscf/orca_leanscf.cpp, line 306, Process 1]: Error (ORCA_LEANSCF): unfortunately, the SCF has not converged. There may be a way out but we have to stop here
  [file orca_leanscf/orca_leanscf.cpp, line 306, Process 1]: Error (ORCA_LEANSCF): unfortunately, the SCF has not converged. There may be a way out but we have to stop here
  [node-112-1500g-1:4177135] [ 4] [file orca_leanscf/orca_leanscf.cpp, line 306, Process 6]: Error (ORCA_LEANSCF): unfortunately, the SCF has not converged. There may be a way out but we have to stop here
  [file orca_leanscf/orca_leanscf.cpp, line 306, Process 6]: Error (ORCA_LEANSCF): unfortunately, the SCF has not converged. There may be a way out but we have to stop here
  [file orca_leanscf/orca_leanscf.cpp, line 306, Process 13]: Error (ORCA_LEANSCF): unfortunately, the SCF has not converged. There may be a way out but we have to stop here
  [file orca_leanscf/orca_leanscf.cpp, line 306, Process 7]: Error (ORCA_LEANSCF): unfortunately, the SCF has not converged. There may be a way out but we have to stop here
- Hans_8DQ2__EF12__origin__La__alpb: complete; energy Eh=-3040.435409827431; reason=none
- Hans_8DQ2__EF12__origin__Dy__vacuum: unavailable; energy Eh=None; reason=no successful native execution receipt
  [file orca_leanscf/orca_leanscf.cpp, line 306, Process 13]: Error (ORCA_LEANSCF): unfortunately, the SCF has not converged. There may be a way out but we have to stop here
  [file orca_leanscf/orca_leanscf.cpp, line 306, Process 13]: Error (ORCA_LEANSCF): unfortunately, the SCF has not converged. There may be a way out but we have to stop here
  [file orca_leanscf/orca_leanscf.cpp, line 306, Process 1]: Error (ORCA_LEANSCF): unfortunately, the SCF has not converged. There may be a way out but we have to stop here
  [file orca_leanscf/orca_leanscf.cpp, line 306, Process 1]: Error (ORCA_LEANSCF): unfortunately, the SCF has not converged. There may be a way out but we have to stop here
  [file orca_leanscf/orca_leanscf.cpp, line 306, Process 18]: Error (ORCA_LEANSCF): unfortunately, the SCF has not converged. There may be a way out but we have to stop here
  [file orca_leanscf/orca_leanscf.cpp, line 306, Process 18]: Error (ORCA_LEANSCF): unfortunately, the SCF has not converged. There may be a way out but we have to stop here
- Hans_8DQ2__EF12__origin__Dy__alpb: complete; energy Eh=-3039.193803581625; reason=none

Collection: /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_global_occupancy_20260923/native_feasibility_retry_v2/COLLECTION.json
No automatic resubmission or expansion is attached. PQQ production is unchanged.
