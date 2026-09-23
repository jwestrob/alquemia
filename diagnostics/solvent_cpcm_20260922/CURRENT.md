# Checkpoint — pilot and matched vacuum qualification running

The CPCM backend is real, but ORCA6.1.1 automatically replaces the special native
xTB mixer with ordinary ORCA SCF. Thus the original archived native vacuum cannot
silently supply the reference. The first actual Ca output confirms the expected
Gaussian-vdW cavity, water epsilon80.151, Ca radius2.772Å and nonzero printed
reaction-field energy. This does not by itself qualify the complete transfer.

- Initial eight CPCM calls: job **1209970**, `pilot_v1/manifest.json`,
  SHA `8d6d6c69faa32fbc99f4e8629676b71a70d45d8b727f02997abc9fc25468d5dc`.
- Eight matched ordinary-vacuum controls: **1209980**,
  `matched_vacuum_v1/manifest.json`,
  SHA `e7cab4f09ee010bf44a2e73eb45a5ea87dd32865a02e658d257693a05af758e0`.
  Genuine explicit MOREAD uses each matching native vacuum GBW. Zero CPCM reruns.
- Both use64CPU/128GiB,8concurrent×8MPI, explicitnode224, noGPU. Scheduler default
  time limit remains; no added project CPU/time limit.
- Initial observations: all four Ca CPCM complete; La still oscillating. One of
  eight ordinary vacuum controls complete, with the others continuing. These are
  operational observations, not a final denominator or discrimination claim.

Finite local observer PID3554671 waits only for these two jobs. Its exact script,
command and receipt are `matched_vacuum_v1/finish_frozen.py` and `MONITOR.json`.
It performs no molecular calls: final matched collection, actual costs/scheduler
accounting, `REPORT.md` here and a dedicated vault note. The wrapper also collects
its phase when execution ends nonzero. No failed endpoint becomes a baseline
substitute. No additional seeded CPCM attempt is automatic.

All paths above are under `workspaces/solvent_cpcm_20260922/`. Final authoritative
output will be `matched_vacuum_v1/final_collection.json`. Do not duplicate the
observer or change running snapshots. If it fails, rerun the recorded command
after checking whether its final output already exists.

No full30 calibration or225fold transfer is submitted. Proceed to those only if
this matched pilot is coherent; numerical qualification currently remains open.
Released MACE+GFN2ALPB and preserved DFT remain unchanged.
