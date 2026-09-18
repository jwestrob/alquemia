# Technical re-execution after the parser-cache fix

This continues the approved goal's implementation/performance work. It does not
add a scientific model or biological analysis: same25cases, physical preparations,
checkpoint, scientific adapters, two-call algebra and frozen calibration. No new
threshold or case selection. Original end-to-end campaign is retained in full.

Matched report-only evidence found62.06s without the cache versus44.67s with it,
with identical scientific output. All54supported endpoint coordinate arrays match
the original strict parser bit for bit. Four real-artifact/cache-failure tests pass.
These checks justify a fresh end-to-end timing of the repaired implementation;
stage timing alone is insufficient evidence of full-workflow speed.

Execute50fresh MACE bound endpoints (two per canonical case, same original order)
using isolated interface_source_v5. No new DFT calls: reuse the unchanged matched
DFT campaign1201562, verifying identical case/source/reference records. One A5000,
16CPU,64474MiB on the same node and pinned software. The manifest schedules MACE
only. No existing energy result satisfies a fresh timing task; parsed-coordinate
reuse within an operation is explicitly distinct from reusing a computed energy.

All extra development cost is counted separately and in the development total.
Report original and repaired execution results separately. Retain the original
median-time criterion, score reproduction tolerance and literal class decisions.
Do not replace a failed repaired run with the original and call it successful.

This is a technical performance regression run under the agreed workflow: no
claim of fresh biological validation, no scientific protocol change, no altered
production/default path. Runtime improvement is a hypothesis until measured.
