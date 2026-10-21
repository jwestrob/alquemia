# Published-recipe central scout — running

21 October UTC: **16602 is RUNNING**, 24 MPI ranks, 96,000 MiB, MaxCore3100MB/rank,
standard_expand, scheduler default memory/time, no GPU. One central EnGrad only.
Completion42618 and alert-only health42619 are live. Exact launch and runtime
receipts are in the matching workspace/run directory. Run submission archive
c1520d88 preserves the observed SCNL correction alongside original inputs.

Preparation16601 passed unchanged on1CPU,7CPU-s,zero engine calls; result6c4200a5,
independent actual review1ceefb4d and archive5fb4489b. Sourceaa62db4f/config1c9f22df/
SPEC0f11bb6a/root release7c7e0dcc are frozen. Root startup review296af6a0 and
independent startup reviewdd8dba98 record the actual allocation and parsed input.

## Actual method discovery

The native executable automatically **switches to self-consistent VV10 (SCNL)**
for EnGrad and disables AutoTRAH. It activates the libXC functional and labels
La's ECP as def2-tzvpd.bas,46core,lmax3 at atom255. Our earlier post-SCF and
built-in-ECP descriptions were incomplete/wrong for this actual gradient run.
Preserve the frozen original SPEC/source/evidence; do not silently rewrite them.
Read ROOT_NATIVE_METHOD_OBSERVATION.json (16e1c029), STARTUP_REVIEW.md and
NATIVE_NL_METHOD_DISCOVERY.md (8fcb2d2e). Historical6.0 automatic switching is
unverified. Compare the same RMS convention: TARGET explicitly uses Cartesian
component RMS, which differs from atom-vector RMS used in some old reports.

Root continues the unchanged exact published input because this still tests
what the published recipe executes under ORCA6.1.1 at the fixed geometry.
It must not be reported as a post-SCF-VV10 result or exact ORCA6.0.0 label replay.
No native scientific result is available yet; initial parsing is not admission.
The implemented success status is deliberately pending actual method review.

After terminal evidence, review full method/basis/ECP/rotation/grid/spin,
SCF residuals, finite forces and actual cost. Any future energy-only pair must
be separately declared to match the actual SCNL Hamiltonian; simply removing
EnGrad is not released. No duplicate call, optimizer or automatic retry.
The old PBE0 directional gate remains failed. Root owns all launches.
