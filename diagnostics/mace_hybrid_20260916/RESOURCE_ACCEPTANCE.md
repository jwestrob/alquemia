# Resource amendment approved, 2026-09-16

Jacob: “then that's an enforced cap. that's fine. let it run with the 195GB.”

This supersedes the initial exact-one-eighth host-memory requirement for the
first allocation. Request and accept **200000 MiB (195.3125 GiB), 28 CPUs and
one H200**. The original twelve scientific tasks, pinned model and implementation,
coordinates, charge/spin, precision, vacuum boundaries and numerical checks do
not change. Original AGREEMENT.md and the scientific manifest remain immutable.

The observed scheduler rewrite is treated operationally as the accepted cap;
the exact server-side mechanism has not been identified. Further reservations
remain subject to the earlier memory-recovery authorization. Additional GPU
shares cannot be assumed to yield additional host RAM: a host-memory retry must
actually receive more host memory before inference. No production/default change.

The older queued script contains the exact-share guard in Slurm's saved copy,
so updating the local script alone cannot unblock it. Replace only that owned,
still-pending job and its task-owned watcher, preserving all prior receipts.
