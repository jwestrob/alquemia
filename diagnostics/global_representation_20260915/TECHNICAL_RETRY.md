# Native xTB MPI startup recovery

The 172-rank-per-endpoint setting was too aggressive for native ORCA startup. RSS rose through 7 TB before SCF. All ranks were active, but that did not establish efficient parallel scaling; startup appears to replicate large per-rank arrays. This is an execution configuration problem, not evidence against whole-protein chemistry.

Only this task-owned job, 1198999, is stopped for technical recovery. Its partial output, inputs, parameter export and cost remain preserved. The retry uses 64 MPI ranks per endpoint, with every scientific input/XYZ byte unchanged. No method, source, assembly, protonation, water, geometry, charge or biological comparison changes. This is within the approved pair and routine technical-retry scope. No compute/time budget applies.

The physical-memory evidence and per-attempt accounting are in workspaces/global_representation_20260915/technical_restart_1198999.json. The new manifest is retry_64rank_v1/global_manifest.json; submission script is run_global_retry.sbatch. Runtime and memory behavior must still be measured.

## Corrected retry planning

Final observed batch peak RSS was 8,053,902,056 KiB. One native startup process reached 32,055,664 KiB RSS. Slurm cancellation could not clear all processes and drained node-344-8t-1 with “Kill task failed”. No node restart, undrain, permission change or other user job was attempted. Job 1198999 used 297,904 allocated core-seconds on the top-level receipt (866 s x 344); later batch cleanup accounting extends beyond the top-level cancellation time and is retained separately, not double-counted as another allocation.

The initial fixed-64-rank retry, job1199003, was cancelled while still pending before any scientific execution: that configuration would also be unsafe on a smaller node. It is superseded by retry_memory_v1 and run_global_memory.sbatch. The new launcher preserves the existing 16-rank endpoint maximum, allows 64 GiB per rank (about twice the observed startup maximum), and plans against 75% of actual node RAM. On a 768GB/64CPU node this gives four ranks per endpoint; on a 1.5TB/112CPU node, eight. This is explicit hardware-fit planning, not a compute-time or accumulated-cost budget. Native full-system setup may still be expensive; no affordability claim is made.
