# V3: reuse4088 and validate exactly one remote terminal OXT

Both4088 attempts reached completed standard protonation, then failed the frozen
AF3 wrapper's zero-heavy-completion assertion: internal missing count0, terminal
count1. Actual outputs identify terminal A:TRP616/OXT, mapped crystal residue582.
The original failed admissions, protonation outputs and truthful manifests stay
unchanged. This is a separately named crystal-only preparation policy.

Reuse those exact outputs; no repeated hydrogen placement. Before carving, require
the raw manifest counts0/1, no added residues, pH7, the exact pinned normalized
source/output, and a heavy-atom set difference consisting solely of A:TRP616/OXT.
Every original heavy atom must retain its exact stored coordinate. Retain both
native Ca and PQQ, and zero H on those three nonstandard residues. No other
terminal/internal heavy completion, PQQ repair or source change is allowed.

Call unchanged pqq_fast_prepare.build_core, complete_expansion and
motion_envelope.prepare_task. Preserve the original protonation manifest with
terminal_count1; do not make it pretend to satisfy the AF3 zero-completion rule.
The separate completion receipt is referenced by source and final admission.

Before admission, require terminal residue616 absent from all context source,
bond and cap mappings, and >4.3Å from all retained physical heavy atoms. Record
actual terminal and OXT distances. Require the existing exact Q9 local state,
paired maps, exact mature construct and spectator checks. The existing structural
graph detects any unsupported covalent geometry. Failed checks remain explicit.

Same two chains, exact original Q9 fragment union, conditional oxidized PQQ3−,
singlets, dry water policy and no scoring/optimization. Both real chains exercise
all recovery checks on one Slurm CPU. No scientific tests/imports on login.
