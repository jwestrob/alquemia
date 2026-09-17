# Trace learned charge updates on the existing real states

Status: approved under the active MACE discriminator goal and standing pilot
authorization. No new labels or physical input selection.

Question: where do the full-protein medium/large predictions diverge? Inspect
initial learned charge sources, the learned weights used to restore total
charge, normalized field inputs, each field-dependent update, and the final
density. In particular quantify cancellation in the sum of charge-restoration
weights, the amount of globally redistributed charge and update amplification.
These are internal diagnostics, not new affinity scores.

Run the same primary qm33/qm36 La/Ca cores and full La/Ca pair for each analytic
checkpoint: **six calls per checkpoint, twelve MACE calls, zero new DFT**.
All states/coordinates/weights/settings match completed jobs 1200470/1200525.
Read-only hooks save detached arrays; they must not change model tensors or
gradients. Verify each call against its uninstrumented archived counterpart:
energy <=1e-6 eV, maximum force <=1e-6 eV/Angstrom, density <=1e-8. These are
the existing memory-implementation equivalence tolerances. Four core calls
must pass before the corresponding full pair executes. Preserve failures.

Capture all atoms and both spin channels, not selected favorable residues.
For each update report raw and restored total charges, absolute/signed sums
of restoration weights, normalized weight extrema/L1 norms, charge correction
norms and field ranges. Retain arrays for source-mapped follow-up. Compare
endpoint changes at the same stages. No clamping, damping, charge neutralization,
new solvent, atom movement or force interpretation as a thermal covariance.

Use existing runner/allocation policy, one A5000/16 CPUs/64474 MiB. Expected
model time from matched receipts: about 2 minutes per medium full pair and
4 minutes per large pair, plus core/startup/trace overhead. No artificial
stopping budget. Separate manifests and cache identities for each checkpoint.
