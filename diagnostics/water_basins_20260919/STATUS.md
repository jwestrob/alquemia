# Finite water basin study status — 2026-09-19

Approved supporting research; production/default and the main PQQ/PLM scanner are
unchanged. All 5 new jobs are complete:1202081and1202085–1202088. Four paired development endpoints,
8192 candidate draws;6250 fall inside the largest physical domain. Rejected draws
remain in the integration denominator and do not need scientific evaluations.
Thirteen real-fixture tests pass, including exact rigid coordinate mapping, H-exchange
symmetry, paired physical domain identity, Haar Jacobian and rejection denominator.

Four native representatives per endpoint ran as each center completed,
overlapping the remaining GPU sampling. This preserves the agreed16 checks and
selection rules. No extra optimization, new water inventory or entropy offset.
All 16 native calculations completed;0/4 endpoints pass wider-domain nativeenergy
checks. No occupancy/entropy correction is available. See REPORT.md and export_v1.
The previous harmonic and minimum qualification failures remain unchanged.

The two1F6Swaters had hydrogen-label assignments differing betweenCa/La.
Uncorrected orientation differences were3.118/3.134rad; physical orientation
differences after equivalent-H mapping are 0.0748/0.0607rad. The mapping preserves
all actual endpoint atom order and coordinates. Common physical domains retain
endpoint internal shapes; largest inherited paired-shape serialization difference
is4.26e-11Angstrom.

## Older comparator collected, no rerun

Native orientation job1201825 completed in18490s on64allocatedCPUs. It belongs to
the earlier hydration-network experiment, not this sampling cohort. Collection:
`workspaces/hydration_network_20260918/orientation_6ip9_v1/collected_opt_v2.json`.

- 6IP9 sourceLa: qualified.
- 6IP9 sourceCa: frozen-coordinate drift failure.
- 6IP9 radial-awayLa: optimization did not converge.
- 6IP9 radial-awayCa: qualified.

Only2/4 qualify; no failed result is silently accepted or rerun. Its long native
optimization remains a development comparator, not a proposed scanner operation.
