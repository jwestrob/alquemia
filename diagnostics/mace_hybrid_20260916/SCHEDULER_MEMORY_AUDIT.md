# Scheduler memory-request investigation, 2026-09-16

Later execution amendment: Jacob accepted200000 MiB and explicitly instructed
the pilot to proceed. See [RESOURCE_ACCEPTANCE.md](RESOURCE_ACCEPTANCE.md).
The investigation below records the earlier observations; its pending user
question is now resolved, and the initial exact-share guard has been removed.

Observed, before any MACE inference:

1. Job1200196 requested one H200/28 CPUs and2063701/8 MiB of host memory.
   Its record showed200000 MiB. Updating our pending job to257962 MiB initially
   changed the record, which subsequently returned to200000 MiB.
2. Cancelled1200196 while pending; accounting confirms zero elapsed execution.
   Replacement1200197 explicitly passed `--mem=257962M` on the command line.
   The immediate scheduler response showed257962 MiB; later responses and
   accounting again show200000 MiB. This is not a unit-rounding difference.
3. No visible partition maximum or standard-QOS MaxTRES memory restriction
   explains the rewrite. JobSubmitPlugins and CliFilterPlugins are null.
   Jacob has no user crontab; no relevant user timer or visible policy process
   was identified. These observations do not establish a particular cause.
4. Controller log `/var/log/slurmctld.log` is mode0600 owned byslurm. System
   journal access is also denied. No attempt was made to elevate permissions.
   Identifying the caller or internal operation that changes the request may
   require administrator access to those records. No messages sent to anyone.

Separate confirmed configuration finding:

```
SelectType = select/cons_tres
SelectTypeParameters = CR_CPU
gpu_h200 SelectTypeParameters = NONE  # inherits global selection setting
```

[SchedMD's memory-management documentation](https://slurm.schedmd.com/cons_tres_share.html)
states that `CR_CPU` without `_Memory` does not track host memory as a consumable
scheduling resource. Accordingly, requesting memory does not guarantee exclusive
reservation of those bytes against other jobs. Cgroup/job memory enforcement is
a separate question; local cgroup configuration was not readable. Reserving extra
GPUs may reduce competing GPU work, but does not by itself guarantee a proportional
host-RAM reservation with this configuration.

This setting does **not** establish why the request becomes200000 MiB. Keep the
two findings distinct. The worker's allocation guard refuses the smaller share
under Jacob's current exact-one-eighth instruction. A question about accepting
the smaller share is pending. No model throughput, GPU-memory result, numerical
check or La/Ca prediction is available from this pilot yet.
