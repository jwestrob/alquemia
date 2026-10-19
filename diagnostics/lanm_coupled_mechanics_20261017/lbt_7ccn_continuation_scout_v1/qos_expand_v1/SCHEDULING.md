# User-authorized alternative Slurm QoS

19 October 2026, 23:22 UTC. Jacob relayed Rohan's authorization to try
`--qos standard_expand`. Live association inspection permits that QoS for
jwestrob/banfield. Root changed only the QoS of the already-pending owned
preflight15192 using `scontrol update JobId=15192 QOS=standard_expand`.
The update succeeded and was verified in the scheduler. No duplicate job was
submitted; the existing watcher, job ID, command and one-CPU request remain.
Slurm automatically changed its default time limit from three to five days.

| QoS limit | standard | standard_expand |
| --- | ---: | ---: |
| CPUs per job | 12 | 24 |
| CPUs per user | 120 | 500 |
| Jobs per user | 10 | 20 |
| Wall time | 3 days | 5 days |

These are configured caps, not reservations or guaranteed availability. Both
QoS entries have priority0. `standard_expand` carries the administrator-set
`OverPartQOS` flag; partition policy still applies where not overridden.
Exact policy/association and before/after scheduler output are in UPDATE.json.

For the already-declared continuation scout, pass `--qos=standard_expand`
explicitly to sbatch before the preserved launcher path. Keep12CPUs/oneH200,
the original pinned scientific configuration and all admission gates. Existing
launcher snapshots and SPEC/READY remain immutable; this scheduling addendum
records the explicit command-line override. Actual preflight still gates that
single molecular run. No account-wide default, node priority, other user's
allocation or scientific protocol is changed.
