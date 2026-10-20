# Bounded 7CCN-La continuation: lower energy, still nonconverged

19 October 2026. Scout15254 completed the separately authorized12 additional
accepted steps and12 distinct composite evaluations. Cumulative counts are
92 accepted steps and92 evaluations. It stopped at `accepted_step_limit`;
**full-force convergence was not reached**.

The run resumed the exact80-step checkpoint with its eight history pairs and
original source geometry. All recorded geometry/state guards passed, without
rebasing the source or changing the energy, states, topology or force thresholds.
This is one La scout, not a completed paired-metal response campaign.

## Increment from the saved checkpoint

| Contribution | After − checkpoint, kcal/mol |
| --- | ---: |
| Whole MACE vacuum | −1.705988 |
| Same fixed-charge ALPB | −2.520990 |
| Total | **−4.226978** |

| Full-system residual force, kcal/mol/Å | At80 steps | At92 steps |
| --- | ---: | ---: |
| Maximum atom-vector force | 9.739698 | 5.674648 |
| Atom-vector RMS | 1.987798 | 1.086248 |

These endpoint values are lower, but the force trace oscillates. No monotonic
force decrease is claimed. The final values remain above the unchanged
maximum0.10/RMS0.03 convergence thresholds. Energy lowering under this model
does not demonstrate improved accuracy, a minimum or metal preference.

The prior saved PBE0 slope supported a conditional local direction at the old
endpoint. It does not independently validate this finite12-step displacement.
The new endpoint needs separate reference evidence; this report does not launch
or authorize further optimization, another metal or a reference calculation.
Conditional terminal chemistry, vacuum electronic treatment with post-added
solvent, Lu-conditioned7CCN and electronic-state/force-accuracy limitations remain.

## Execution and evidence

Actual new calls were12 whole-MACE and12 fixed-charge ALPB evaluations. There
were no DFT, GFN, MM, opposite-metal or other-source calls. One model was loaded;
the saved application records6.373613s loading,4.073442s across MACE calls and
0.725242s across ALPB calls. These timings are distinct from scheduler allocation.

| Job | Allocation | Elapsed | Allocated CPU-seconds | Allocated GPU-seconds |
| --- | --- | ---: | ---: | ---: |
| Preflight15192 | 1CPU, no GPU | 10s | 10 | 0 |
| Scout15254 | 12CPUs, 1GPU | 26s | 312 | 26 |
| Total | | | **322** | **26** |

Actual `standard_expand` QoS, H200 partition, allocated resources and scheduler
records are preserved in [COST.json](actual_result_v1/COST.json) and its linked
preflight evidence. Allocation cost is not measured utilization.

The [lossless result](actual_result_v1/execution/RESULT.json.gz) has original
SHA256 `c445cf4b622fa2dbc7d7c47456e4b20b7c8db1a0b7718365fab2310c3c90a8c6`.
The archive preserves all12 new POINTs, full gradients, raw MACE and ALPB outputs,
requests, coordinates, all14 checkpoint files and the complete trace:210
execution files, each roundtrip-verified. Historical parent data remain linked
through their exact pins. [EVIDENCE_PINS.json](actual_result_v1/EVIDENCE_PINS.json)
also records source, optimizer, configurations, decision, reviews, launcher,
submission and receipt snapshots; [RUNTIME.json](actual_result_v1/RUNTIME.json)
retains the executed environment metadata.

All reported quantities are copied from completed artifacts. Archival/reporting
performed no scientific recomputation. Root owns completion acknowledgement,
integration and any subsequent scientific decision.
