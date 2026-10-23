# One-time tighter Dy-acetate progress review

At 06:25–06:28 UTC on 22 October, two roles pass exact saved raw admission (32 native checks each), while the ligand-with-ghost-ion role has no terminal result. The three-role CP cycle remains unavailable; the original 0.25 kcal/mol gate and cluster hold are unchanged.

| Role / job | Saved admission | Native wall | Allocated CPU |
| --- | --- | ---: | ---: |
| Complex / 18516 | Complete | 4241.074 s | 33928.593 CPU-s |
| Ion with ghost ligand / 18517 | Complete | 59.203 s | 473.627 CPU-s |
| Ligand with ghost ion / 18518 | Incomplete | Not terminal | Not terminal |

The completed outputs retain Energy Check stopping and failed printed MAX/RMS density and DIIS criteria. CC residuals are 8e-9 and 3e-9; terminal triples energy changes are 5.32e-9 and zero against 1e-8. The complex T1 diagnostic is 0.020818158. All declared cutoff, basis, core, electron and producer checks replay exactly; this does not qualify physical accuracy or electronic-state uniqueness.

Job 18518's captured output ends at `Calculating RI/C V-Matrix + Cholesky decomp....`, before any SCF/CC/triples residual. It had not changed for about 113 min 45 s and stderr was empty. No explicit native failure appeared. The native input reports 189 orbital and 718 auxiliary-C functions and MaxCore 3000 MiB.

The single scheduler/process observation confirms eight live `orca_startup_mpi` ranks, four on each of node-64-768g-5 and node-64-768g-6. Every rank is runnable with about 99.98–99.99% lifetime CPU and 72–82 MiB RSS. The ORCA/mpirun/srun/orted parent chain and exact endpoint working directories match. The allocated two-node layout is real; rank count is correct. Launcher-step `NTasks=1` does not mean only one MPI rank ran.

This excludes idle or dead workers at the observation, but cannot distinguish useful RI/C work from a busy MPI wait. Low observed RSS gives no indication of memory pressure. No cancellation, retry, new call, threshold change or repeated polling was performed. Root retains the execution decision; no incomplete CP contrast is reported.

Completed native allocated cost is 34,402.221 CPU-s (rounded scheduler cost supplied by root: 34,448 CPU-s). The eight running ranks had accumulated 55,933.46 CPU-s at their process snapshots; this is observed process CPU, not a terminal allocation cost. Raw snapshots, source pins, exact process parentage, all printed SCF flags and the transcribed one-time scheduler response are retained in `REVIEW.json` and its linked evidence.
