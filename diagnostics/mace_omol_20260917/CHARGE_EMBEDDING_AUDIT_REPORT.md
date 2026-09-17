# This checkpoint collapses several distinct input charges

The pinned100M checkpoint was inspected directly, with no molecular energy
forward and no weight changes. Its categorical table accepts201 charge inputs,
but contains only18 distinct joint charge/singlet representations:

- All charges from−100 through−6 have exactly identical raw charge embeddings
  and joint charge/spin features.
- Charges−5 through+10 have distinct joint features.
- Charges+11 through+100 form a second exact equivalence group.

These are byte-identical checkpoint values, not a threshold applied by our
preparation code. The native input receipts retain the actual requested charges.
The audit does not establish why these weights were produced, how many training
examples supported them, or a different biological label. The dataset's
reported charge range and this checkpoint's actual distinguishing capacity are
different properties.

This explains the apparently successful MxaF spectator result: its original
La/Ca charges−9/−10 and spectator charges−8/−9 all fall in the same embedding
equivalence group. The other four cases cross distinguishable categories.
Their4–15kcal/mol shifts, including the XoxF/MxaF ordering reversal, therefore
expose a consequential global-conditioning ambiguity. A near-zero shift for
MxaF is not independent evidence that the representation handles remote ions
physically.

Combined with the [saved-output locality audit](INTACT_LOCALITY_REPORT.md),
this identifies a concrete limitation: remote coordinates cancel, while the
chosen system's net charge can alter all local features, and some very
different net charges are indistinguishable. This is not fixed by a new aquo
offset or recalibrating the old production bands. No arbitrary charge clipping,
neutralization, per-protein checkpoint choice or fitted correction was applied.

The [native architecture explanation](https://github.com/ACEsuit/mace/discussions/1206)
describes global charge/spin conditioning. Exact numerical findings above come
from our pinned checkpoint and installed implementation, not an inferred
training history. Raw arrays, all201 categories, source hashes and execution
timing are in `workspaces/mace_omol_20260917/charge_embedding_audit_v1/`.
The201 embedding-block rows are parameter inspection; zero new molecular
endpoint evaluations, gradients, solver calls or training jobs occurred.
