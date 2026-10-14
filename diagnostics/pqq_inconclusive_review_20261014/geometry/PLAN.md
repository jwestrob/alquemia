# Inconclusive-source geometry comparison

Use all48 source/context cases in root's pinned REVIEW_INPUT, including all nine
target inconclusives and every declared same-protein peer. This is a descriptive
analysis of existing structures and proposals, not an optimization or a new
geometry-admission/score rule. Source conditioning, occupancy and context identity
remain separate comparison dimensions; related folds are not independent proteins.

Reuse the existing PedE geometry functions for actual donor identities/distances,
typed coordination/denticity under the unchanged3.1Å definition, PQQ bond lengths
and ring planarity, donor torsions, mapped origin-to-proposal motion, covalent-bond
changes and AF3 confidence. Use the passed LOCAL_ONLY physical-source mapping
reconstruction before these functions: the capped QM XYZ and full physical graph
do not have the same indexing. Report saved-vs-kinematic residuals explicitly.
Saved coordinates determine mapped atoms; no new coordinate is optimized.

Inspect origin, Ca proposal and La proposal for every case, identify which each
metal actually selected under the new scorer, and retain recorded optimizer
boundary/residual/convergence fields. Missing flags stay missing. No new bound
hit is inferred from an unavailable optimizer record. Rank the twelve closest
mapped nonbonded heavy-atom pairs excluding1–2/1–3 graph neighbors and the metal;
these distances are not a new clash criterion or a source-exclusion rule.

Reuse exact previously completed PedE raw/confidence metrics when raw source,
selectors, roles and helper implementation agree. Parse remaining original AF3
sources; use the actual audited active/spectator selectors for two-Ca sources,
including Q9 sample3's active chainD. Spectators remain present in full sources.
Confidence supports structural description; it is neither a metal-affinity label
nor a calibrated uncertainty. Predicted PQQ bond lengths alone do not establish
redox/protonation state. Do not transfer a mechanism reported for another sample.

All numerical geometry work and scientific imports run in a one-CPU Slurm
allocation. Metadata-only preparation merely pins real saved artifacts. No
molecular calls, new folds, score changes or production edits. Root owns score
margins/energy decomposition and integrates this report with those results.
