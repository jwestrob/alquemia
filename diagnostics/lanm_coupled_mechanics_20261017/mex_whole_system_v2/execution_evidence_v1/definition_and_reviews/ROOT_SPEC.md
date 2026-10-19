# Coherent whole-system mechanics: finite intervention

Declared 19 October 2026 before new evaluations. Root owns submissions under
the active coupled-mechanics goal. Preserve the failed hybrid responses and the
closed Mex DFT derivative candidate. This is a research candidate, not a
production change.

## Question and model

Does replacing the partitioned internal-energy model with one whole-system
MACE energy remove the unsupported scaffold work, without introducing a large
dependence on the declared aggregate spin?

The candidate is exactly:

`E(R) = MACE_OMOL_vacuum(all physical atoms; charge, multiplicity) + ALPB(fixed FF charges; R)`

No MM, GFN embedding difference, extra dispersion, fitted correction, synthetic
cap or additional interaction term is added. This removes the partition and its
mechanical boundary; it also changes the treatment of long-range interactions
and spectator ions. Those effects are not separable from this one intervention.
Whole-protein model accuracy and training-domain applicability are unverified.
Fixed-charge ALPB remains an additive solvent term, not self-consistent
polarization of the vacuum model.

Whole256 MACE+ALPB already exists for LBT and passed its declared directional
derivative tests in job12958. Reuse those qualifications and all completed LBT
searches. The inspected 7CCO hybrid-motion DFT comparison is consumed development
evidence; agreement there alone does not validate a replacement.

## Mex fixed-motion panel

Freeze the exact start and terminal coordinates of job12970's La full-coordinate
search, without selecting on a new energy. Preserve all2,049 physical atoms,
167 waters, the EF3 target and the three Nd spectators at EF1/EF2/EF4. The
80-step endpoint is admitted but not converged. Both metals use the identical
two coordinate arrays. Reuse the four existing matched ALPB receipts from the
completed common pool12981; no new solvent or geometry calls are needed.

Whole-system charge is +2. The old local-core multiplicities La1/Dy6 do not
specify the aggregate state when all three Nd ions enter the MACE graph.
Freeze two explicit conditioning hypotheses, preserving the same hypothesis
before and after the motion:

| Hypothesis | La multiplicity | Dy multiplicity |
| --- | ---: | ---: |
| Aligned counting: three Nd 4f3 moments and target | 10 | 15 |
| One Nd moment reversed in spin counting | 4 | 9 |

These are aggregate model inputs. A single total-spin category does not specify
localized spin populations, identify a reversed ion, establish coupling or prove
that either hypothesis is the physical electronic state. Check actual nuclear
electron parity, model element coverage and embedding-category support before
any forward. Do not select the lower-energy hypothesis after inspection.

Denominator: two geometries × two metals × two spin hypotheses = eight fixed
whole-system MACE energies and gradients, with four reused ALPB calculations.
The first scout is the before geometry, La, multiplicity10. Reuse this exact
scout in the remaining panel. No optimization, additional spin scan, refolding,
new DFT or geometry repair belongs to this Mex panel.

## Execution gates and finite budget

Use the existing pinned MACE-OMOL checkpoint, float64, native adapter/kernel and
software environment. Slurm: one H200 GPU, at most12 CPUs, existing scheduler memory
defaults, no time/memory override. A source-only adapter may generalize the
physical-size/state validation but may not change model arithmetic. Verify
actual graph charge/spin, identities, coordinates, element support and unchanged
parameters; save the full force array and actual peak GPU memory and runtime.

Zero-call preflight checks the complete eight-cell manifest, four solvent
receipts and parent search/pool joins. The scout proceeds only after it passes.
Continue the remaining seven only if the scout returns finite energy/full forces,
all state gates pass, model-forward time is at most600 seconds, and peak allocated
GPU memory is at most90% of the device total. This is a continuation gate, not
a scheduler runtime override. Model-loading time is separately recorded.
On failure or OOM retain the exact attempt; no smaller system, changed precision,
different model or silent energy-only fallback. Root reviews any technical
recovery separately before it is declared. No duplicate completed cell.

The final saved comparison reports per-metal before/after work, Dy−La work
contrast, full force changes and component contributions for both hypotheses.
Use the existing1 kcal/mol reporting resolution for changes in work/contrast
between hypotheses; it is not an error bar. A larger dependence blocks treating
one chosen spin input as a unique prediction. Even agreement does not establish
physical state or whole-protein accuracy. Old hybrid work is a comparator, not
an independent reference. Do not compare absolute energies of different models.

## Independent LBT challenge

Freeze the start and exact terminal of completed7CCN-La whole-MACE search12968,
algorithmically selected before new DFT scoring. Existing four paired MACE
energies and ALPB receipts are reused. A separate versioned specification binds
four whole256 PBE0-D4 vacuum reference cells: before/after × La/Dy, with the same
fixed-FF ALPB added. Keep the established basis/ECP, physical/effective spin,
native SCF admission and all reference limitations. Missing saved Dy MACE forces
remain missing; do not recompute existing energies just to fill them.

Start with one before-La reference scout; review actual native admission and
cost before submitting the remaining three. No optimizer or derivative/grid
sweep is authorized by this declaration. A result that contradicts the new model
must be retained and stops unsupported optimization. Small work differences
below1 kcal/mol are below the existing interpretive resolution, not proof of
equivalence. Report all four cells and failures, never only favorable survivors.

## Decision

This phase implements and challenges a coherent candidate. Successful execution
does not finish the mechanics goal. Natural LanM transfer still requires a
defensible common state/interaction treatment and accuracy evidence. No restart
of the challenged hybrid optimizer or closed Mex derivative experiment follows.
