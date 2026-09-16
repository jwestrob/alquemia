# Local neural-activation memory recovery

Continuation of Jacob's approved memory work in MEMORY_AGREEMENT.md; no new
scientific analysis/model. The pair-blocked version passed all kernel and full
core checks. Full-system jobs1200368(native) and1200369(host offload) failed
before electrostatics in `RealAgnosticResidualNonLinearInteractionBlock`:
`tp_weights * cutoff` requested another10.95GiB tensor on the24GB A5000.
This identifies a separate live neural-activation memory requirement.

Preserve the exact learned interaction, process4096 neighbor edges per block,
and checkpoint intermediate activations for backward recomputation. All edge
messages and densities accumulate before node normalization/nonlinearities.
No neighbor cutoff, weights, precision, physical input or model is changed.
The same four core endpoint retries must pass the existing1e-6eV energy,
1e-6eV/A force and1e-8 density thresholds before full inference. Small edge
tiles128 are recorded for core verification so those actual inputs cross block
boundaries; full-system edge blocks are4096. Both values are in the manifest.
All failed attempts and earlier implementation snapshots remain immutable.

The eight original full-system tasks remain the execution target, initially
on oneA5000/16CPU/64474MiB with the same technical recovery policy. No new DFT,
large-model inference, biological case, relaxation or training is introduced.
