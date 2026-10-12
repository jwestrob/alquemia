# Frozen Ca-conditioned preparation extension

12 October 2026. Prepare the 107 source/context requests already inventoried.
All125 Ca-conditioned sources remain in the denominator: reuse three exact
expanded probes, preserve15 original PQQ exclusions, attempt107 new preparations.
No model energies, optimizations, force evaluations, folds or new protonation.

For each protein use exactly the archived canonical La-three-source4.3Å union,
including its original first-member state anchor. Do not discover a new union,
change memberships, add waters or select a context from a favorable result.
Reuse current `motion_envelope.prepare_task` and its underlying graph builder,
physical maps, original-context replay and paired-coordinate checks. Their
archived implementation pins must match. Pin the current imported source
modules and preparation runtime as well.

Before107 new preparations, reconstruct the three already admitted Ca probes in
new validation directories. Compare exact state/protein keys, atom/cap counts,
charge/multiplicity, source graph/cap maps and endpoint coordinates against the
saved artifacts. Preserve the existing1e-12Å coordinate-copy tolerance; every
discrete mapping field remains exact. A failed replay blocks the107 batch and
is reported without loosening tolerance or editing the archived data.

Use at most two workers, each restricted to one numerical thread. New rows are
admitted only when the original preparation, physical mapping and the archived
group-anchor state/protein identity pass. Retain failures explicitly. Stage
outputs and process timing are immutable under the new workspace.

The future feasibility scout is the lexicographically first newly prepared,
admitted source. This source-selection rule is frozen before coordinates are
prepared and without inspecting new model outputs. No scout molecular call is
authorized by this preparation operation. Retain original class labels only as
known-reference metadata; no calibration or decision is made here.

Version2 fixes a coordinate-comparison omission revealed by the first real-probe
replay, before any new-source preparations: original-cap vectors in core_links
are coordinate-valued copies, so receive the same1e-12Å copy tolerance as other
coordinate arrays. The observed discrepancy was below9e-16Å, with identical
state keys and endpoint XYZs. All noncoordinate link metadata stays exact.
No scientific admission tolerance, geometry or chemistry is changed. Preserve
failedv1 and its exact implementation snapshot; do not reuse its partial output.

Version3 completes that same geometric comparison: inferred cap lengths are
also coordinate-valued (observed3e-16Å discrepancy), and successful numerical
map-check residuals are derived diagnostics rather than graph identity. Their
observed cross-runtime residual difference is2.3e-16; compare these at1e-12
while retaining both original pass flags and unmodified builder check gates.
This correction changes only replay comparison; source builders and coordinates
remain untouched. Failedv2 and its immutable implementation snapshot are kept.
