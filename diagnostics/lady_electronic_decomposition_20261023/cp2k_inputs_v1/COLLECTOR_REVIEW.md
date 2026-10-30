# Independent collector semantic review

**Status: passes source review with the interpretation limits below. No concrete sign, energy-convention or failure-handling correction requested.** This review read the collector, frozen ledger, relevant parser paths and active user attachment. It did not execute/import COLLECT or PARSE, read live pilot results, recompute a scientific result, or submit anything.

## Exact reviewed sources

| Source | SHA256 |
|---|---|
| `COLLECT.py` | `813b6e13a6afbd5d7137bface4e81bfff2bad257b1acb967f8d57c8e5ceb0366` |
| `ENERGY_LEDGER.json` | `24c699cc97c3c55857c56a5ee07dda6213eeaf12e59af69395f23bbd05b732c5` |
| `PARSE.py` | `743ddda796bb469ad21eb5329e0b7df5f02ebd29e34f319b11419f4f3cdf6202` |
| Active attachment `pasted-text-1.txt` | `45b8ff8d178fdabc245d6f55574c55f2902813d07de63cc0272f43423ba93e04` |

## Signs and constrained energies

The fixed coefficient map is correct:

`DD[X] = X(Dy,mixed) − X(La,mixed) − X(Dy,aquo) + X(La,aquo)`.

A negative value favors Dy relative to La on water-to-acetate substitution. The collector checks that this map equals the frozen ledger map and requires all four named cells. It compares the original model using the same convention; the ledger's recorded original model-minus-reference discrepancy is about −6.810868 kcal/mol, not merely the small positive reference contrast.

For each cell, R is the sum of isolated metal and complete-shell electronic references, F the antisymmetrized molecular-guess energy, P the converged block-constrained energy and U the converged fully delocalized energy. The code correctly defines FRZ=F−R, POL=P−F and CT=U−P. Absolute full and CT-suppressed totals are U+D and P+D, with the **same archived full-system D4 scalar D added once**. Therefore their DD difference is DD[CT], with no dispersion reassignment or reference subtraction hidden in that comparison (`COLLECT.py:158–176,234–269`).

The full absolute-energy identity is DD[FRZ+POL+CT+R+D]. The electronic interaction is U−R; it is not silently substituted for U. Observed repeated-reference differences and DD[R] are retained, including actual fragment sums versus the printed native sum. In exact matched calculations, repeated ion references cancel across motifs and each shell reference cancels across metals. The code appropriately reports numerical deviations instead of forcing this cancellation.

## Dispersion, references and basis limits

Missing isolated-shell D4 references prevent a true absolute D4-corrected per-cell interaction energy. Their matched quartet DD cancels under the ledger's same-shell geometry, charge, basis and state conditions. That cancellation does **not** mean the four full-parent D4 corrections cancel: DD[D] is explicitly reported. The direct printed D4 values, final post-COSX electronic energies and archived total tokens remain distinct; their small printing residual is preserved.

This is the declared own-fragment-basis ALMO convention. No ghost-basis counterpoise reference, independent fragment deformation energy, relaxed geometry, or missing electrostatic/Pauli split is invented. Whole-shell internal delocalization remains allowed. Release of the metal–shell AO-block constraint can include finite-basis relaxation/BSSE effects; CT here is the energy lowering under that defined release, not a unique observable or measured transferred electron count. Combined FRZ is the available quantity.

## Numerical interpretation

Four per-cell 1e-4-Eh archive-agreement gates imply a triangle-inequality bound of 4e-4 Eh, or 0.2510037896 kcal/mol, **only on the corresponding full cross-code DD disagreement**. The collector correctly states that this is not a CT/component uncertainty, SCF error bar or physical-accuracy bound. Its `observed_absolute_CT_DD_above_four_cell_comparison_scale` boolean must not be interpreted as a significance test or certification of CT precision. Native residuals, variational ordering, printed precision, telescoping arithmetic, archive density flags and final COSX shifts answer different questions.

The saved outputs support the requested final assessment against the model discrepancy, but do not automate it. Root must judge whether actual convergence/reference residuals permit a material-CT conclusion, explicitly retain any unquantified component uncertainty, and provide the user's one evidence-supported recommendation. A substantial CT term alone cannot identify missing CT as the force-field defect or establish a binding free energy.

## Failure and provenance behavior

The collector requires four successful terminal runner/parser receipts, every mandatory parser gate, exact configuration/source/ledger/state/deck/raw-output joins, distinct parent jobs, and pilot admission for the remaining cells. It checks saved raw evidence lines against the pinned output and reconstructs components from the saved decimal energy tokens. It never invokes the parser or native program. A missing, failed, changed or inconsistent slot yields a nonpassing result and no quartet contrast; partial per-cell evidence remains visible. New output directories are exclusive (`COLLECT.py:97–156,190–234,290–318`).

Source admission here is not actual quartet admission. The collector deliberately leaves cumulative recovery authorization/budget, earlier failed attempts, build costs, scheduler parent accounting and peak RSS to root's retained ledger. Its four-successful-receipt cost is explicitly limited to those receipts. It neither erases prior OOM work nor releases further calculations. No extra approval or scientific experiment is introduced by this review.
