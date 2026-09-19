# Coupled water response complete — 2026-09-19

All9new allocations complete:27native DFT endpoints and2828recorded MACE objective
calls. All11trial steps lower native energy and pass energy checks;3/4final
endpoints satisfy water-coordinate minimum criteria.1F6SCa remains nonstationary
and retains its initial soft-mode failure;7/8direction checks pass overall.
Thermal extent exceeds the local domain; entropy/occupancy remain unavailable.

Read REPORT.md, export_final_v1/TABLES.md, COSTS.json and COMMANDS.md. Exact new
allocation178496core-s/840GPU-s.71real-fixture tests pass, zero skips. Baseline/
default/PQQ unchanged. Two recentering rounds complete; no further identical
round is scheduled.6IP9La is explicitly reused from round1. Preserve all earlier
collections/snapshots; strict coordinate equality recovery added no DFT calculation.
Final analysis implementation_v4 is frozen. Older comparator1201825 remains live;
preserve and collect when finished. No production rescore, promotion or push.
