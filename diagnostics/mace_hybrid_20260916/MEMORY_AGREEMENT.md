# Approved blocked-electrostatics implementation

2026-09-16. Jacob: "Go for the memory work. That will transfer to
MACE-POLAR-large and help if we have trouble fitting that into H200 later."

Question: can the already approved full 9,141-atom 1H4I system execute on
the A5000 when realspace electrostatic pairs and their derivatives are processed
in bounded blocks? Same medium model, float64, coordinates, charges/spin,
Gaussian widths, dipole offsets, epsilon, exclusions and isolated boundary.
No physical cutoff, truncation, changed precision or new scientific model.
Large-checkpoint inference remains outside this execution step.

Implement pair sums with a custom analytic first-order backward which recomputes
each block. Preserve both coordinate and learned-charge derivatives through all
field updates. No full pair graph, pair distance matrix or saved quadratic
autograd tensors. Start with512 by512 displaced-charge blocks; tile size is a
recorded implementation parameter, not a physical parameter. Higher derivatives
and training are explicitly unsupported until implemented/validated.

Verification before full inference: original vs blocked field, energy and
derivatives on the four actual completed core densities/coordinates, with tile
sizes17 and64 to force block boundaries. Numeric tolerances before outputs:
absolute1e-8 plus relative1e-10 for these kernel comparisons. Run all four
complete core energy/force calls and compare with their archived GPU results:
absolute energy1e-6eV, force1e-6eV/A, density coefficient1e-8. These are much
tighter than the approved physical invariance criteria; no tolerance tuning.
Check parameter/buffer identity. Rewritten kernels are tested against original
autograd, not fabricated successful outputs or numerical DFT derivatives.

Then try the same full protein on one A5000/16 CPUs/64474MiB host RAM. The eight
remaining primary/repeat/translate/rotate tasks retain their original manifest
and frozen physical checks. Native autograd first; host-offload retry if needed,
under the existing approved technical recovery. If CPU memory must grow, use
the already authorized appropriate allocation, with measured actual RAM.
No project compute/time stopping budget; cluster QOS imposes its7-day limit.

Expected cost: core validation minutes; full-system runtime/memory unmeasured.
Blocking reduces memory, not quadratic pair work. Preserve all failed attempts,
reference receipts and isolated software. New versioned implementation/cache;
baseline/default and old experiments remain unchanged. H200 job1200309 remains
queued initially; once the replacement is validated, coordinate only our own
job/watcher to avoid duplicate full-system execution. Report actual costs and
whether full-system affordability and existing numerical gates pass separately.
