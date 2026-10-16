# Selected PQQ decision policy — 16 October 2026

Jacob selected the **canonical midpoint**, then authorized documentation, preservation,
commit and push. This closes the present PQQ development round and returns research
priority to La/Dy profiling. Inclusion in the PLM East River manuscript is **undecided**;
the geometry-method alternative remains available. This decision adds no new comparison
of the geometry method.

## What is preserved

- Exact policy: [SELECTED_POLICY.json](SELECTED_POLICY.json).
- Calibration and leave-group-out fits: [FROZEN_THRESHOLDS.json](FROZEN_THRESHOLDS.json).
- Every scalar context result: [CONTEXTS.tsv](CONTEXTS.tsv).
- Existing overlapping protein summaries: [TRIPLES.tsv](TRIPLES.tsv).
- Completed comparison, implementation and replay command: [REPORT.md](REPORT.md).

The raw-score threshold is **−405457.73794054205 kcal/mol**. Lower is Ca; higher is
La; exact equality is inconclusive. Missing or inadmissible inputs remain unavailable.
This threshold belongs only to the pinned corrected-parameter, frozen-charge-solvent
common6 PQQ score. Keep the two calibration-sensitivity flags beside the midpoint calls.
Do not transplant the threshold to another scorer, LanM or La/Dy.

The completed test gives **257 correct / 0 wrong / 0 inconclusive** on the consumed
source contexts. Two calls remain calibration-sensitive. All 94 available existing
three-source summaries are correct; six preparation exclusions remain unavailable.
These contexts and overlapping summaries are not independent biological replicates.
No broad or prospective validation follows from this benchmark.

The selected rule changes the interpretation of existing scores; this preservation
commit changes no energies, geometries or historical production entrypoints. The
archived band-based results and all earlier failed experiments remain intact. Original
molecular outputs stay in their existing workspaces; the portable scalar evidence and
its hashes are committed here. No molecular calculation was rerun for preservation.
