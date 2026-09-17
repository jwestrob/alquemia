# Inspect the pinned checkpoint's charge representation

Declared2026-09-17 after the disconnected-sodium test failed. The intact MxaF
score was effectively invariant whereas the other four shifted4–15kcal/mol.
Investigate that difference by inspecting the existing checkpoint, without
changing it or evaluating any new molecular energy.

Read the pinned native100M checkpoint and installed GenericJointEmbedding
implementation. For every representable charge category(-100..100), record the
raw charge embedding's norm, zero/nonzero count and hash. At the already fixed
singlet input, compute the native joint charge/spin embedding for one graph
index and record its hash, norm and differences between adjacent categories.
Record exact equivalence groups, without inferring training frequencies from
weight values. Inspect all categories, not only those producing useful scores.

Pin the checkpoint, embedding source, inspection script, specs and output arrays.
These201 small embedding-block evaluations are parameter/representation
inspection, not201 MACE molecular forward/energy evaluations. No coordinates,
forces, scores, fitting, parameter changes, new checkpoint or class decisions.
Use the existing environment on CPU and record time/memory. Explain any exact
degeneracy as an observed property of this checkpoint; do not label it an
implementation bug or invent its training history without evidence.

The installed source and [maintainer explanation](https://github.com/ACEsuit/mace/discussions/1206)
establish that total charge/spin are globally embedded into node features.
Whether individual categories carry distinguishable information is the question
here. The ongoing canonical calculation and original five-case results remain
unchanged and retain their frozen criteria.
