# Exact MACE short-component engine: validated

All 34 real evaluations completed in job1200736. All 50 declared checks pass:
34 energy-equality checks against saved full-forward components, four rigid
transformations, and12 physical-direction gradient checks on the archived GGR
geometries. No new DFT calculation was needed. The production baseline is unchanged.

The engine stops the pinned PolarMACE forward immediately after `scale_shift`,
sums its local readouts, and differentiates that exact scalar. It uses the same
weights, float64 arithmetic and memory adapter. Maximum component-energy error
is4.547473508864641e-12eV; maximum rigid-force error is6.497156901374979e-8eV/A.
The temporary hook is removed on exit. This component has no explicit net-charge
or spin response. Its forces are derivatives of the short scalar, not full MACE
forces. Missing density and charge-response fields are null, not fabricated.
A short-component receipt cannot satisfy a total-energy task.

This is an engineering and derivative result, not a successful discriminator.
The earlier saved-component ordering screen still fails both alpha/GGR tests.
The component is now usable as a separately tested mechanical contribution in
an explicitly defined DFT-anchored model. Relaxation and classification remain
unavailable until that model passes its own physical and predictive tests.

## Actual cost

One A5000,16 allocated CPUs and64474MiB host RAM (one eighth of host RAM).
Job wall481s; allocated GPU481s and CPU7696core-s; actual reported CPU583.343s.
Summed measured inference257.83916721120477s; process/model startup is additional.
Peak CUDA allocated9,080,467,968bytes; reserved11,211,374,592bytes. Slurm batch
MaxRSS1,128,916KiB is its reported host-memory statistic, not CUDA memory.
Full-protein inference per endpoint: alpha5.52–5.81s, GGR13.50–13.57s,
canonical PQQ26.44–27.01s. These times include the short-component analytic forces.

Source collection, receipts, checkpoints and raw accounting are pinned in
[result.json](result.json). Three actual-artifact tests passed in2.159s before
this report. Scientific integration checks above actually ran on the GPU;
parser/cache tests do not stand in for scientific execution.
