# Mex motion and mechanics diagnostic: frozen selection rule

## Purpose and inputs

Use the two completed broad-coordinate Mex La/Dy searches, their exact restricted
starts, physical maps, Cartesian requests, native MM groups and full gradients.
No molecular evaluations, optimizations, repaired structures or new reference
outcomes enter this analysis. Root runs the saved-data diagnostic on one Slurm
CPU only after both actual search results exist. An unavailable or state-changing
Dy result leaves segment selection unavailable; La alone never selects a segment.

The full system remains 2,049 physical atoms: 1,544 peptide atoms, 501 retained
water atoms, the exchanged EF3 target and three Nd spectators. The current local
electronic region and boundary map are inherited exactly. Source metadata about
obsolete spectator mechanics is not used: the adopted Cartesian state and MM
ledger specify moving Nd monopoles with their explicit LJ interactions.

## Alignment and frozen segment selection

For each metal separately, fit its terminal to its own restricted start by a
proper, equal-weight Kabsch rotation using **all protein heavy atoms**. No
reflection is allowed. Apply that same rigid transform to the complete terminal
physical system. Preserve both original coordinate pins, centroids, rotation,
translation, singular values and raw/aligned per-atom displacement. Alignment is
only a coordinate frame for motion analysis; any future reference uses the raw
physical endpoint and corresponding full environment.

Enumerate chemically contiguous three-residue peptide windows using the saved
native C–N bond graph, preserving chain/residue/insertion identifiers. For each
window, sum squared aligned displacement over its protein heavy atoms **outside
the real electronic region**, first per metal and then La+Dy. Do not normalize by
atom count or use energy/preference signs. Rank descending by this pooled sum;
exact ties use the ordered `(chain, numeric residue ID, insertion code)` triples.
The first nonzero window is the sole candidate for subsequent feasibility review.
If every score is zero, selection is unavailable. All windows and both metal
contributions remain in the output. Report QM overlap, cap-anchor overlap and
whole-window atom inventories; partial QM residues are visible, not silently
excluded from the window.

This rule is frozen before new scaffold-reference outcomes. Existing independent
DFT derivative outcomes at original coordinates do not affect selection.

## Mechanics and localization reporting

- Report each global MM-group work separately: `retained_bonded`, `MM_LJ`,
  `MM_Coulomb`, `QM_MM_LJ`. Reproduce their sum and the saved search MM work.
  These are global energy terms, **not residue energy decompositions**.
- Map the actual group gradients into the same physical atom order. Replay
  their sum against the assembled MM gradient, and all four composite components
  against the full gradient (absolute tolerance 1e-10 kcal/mol/Å for reporting
  summation only). Replay energy sums to 1e-7 kcal/mol. No original admission
  or convergence tolerance is changed.
- Report total and per-component residual-force locations at both endpoints,
  including all atoms and residues, with named top-20 atom listings. Retain full
  search force norms and finite-search termination; no equilibrium claim.
- Separate peptide, external peptide, real QM peptide, water oxygen, all water,
  target and spectators. Report counts, squared displacement, RMS and maximum
  motion in both raw and aligned frames. Water/internal motions and peptide
  motions are not combined into a peptide-only explanation.
- Provide the full residue contribution distribution, top-window share, number
  of residues accounting for 50/80/90% of outside-QM peptide motion, and effective
  participating-residue count `(sum w)^2/sum(w^2)`. Coverage percentages are
  descriptive, not admission thresholds.
- For the selected window, report nearest real-QM heavy-atom distances at both
  endpoints, peptide-edge separation from existing QM residues, and the number
  of peptide-connected components in the proposed QM-residue/window union.
  Report endpoint gradient projections along its aligned displacement; these
  are first-order projections, **not integrated local energies or work**.

No arbitrary localization pass threshold decides whether one region can explain
the global release. The selection nominates one local stiffness check. Its
coverage and connectivity describe representativeness; a defensible expansion
still requires a complete donor/charge-boundary and local hydrogen-bond-network
review. Highly distributed motion can therefore be reported without pretending
that the top three residues explain all global Coulomb/LJ lowering.

## Artifacts and execution

`MOTION.py --stage bind --config INPUTS.json --output CONFIG.json` hashes the
provided La/Dy result paths and plan, and freezes the source/spec. This stage is
metadata only. INPUTS has `plan` and `search_results: {La: path_or_pin,
Dy: path_or_pin_or_null}`. Do not pass an anticipated result as completed evidence.

`MOTION.py --stage run --config CONFIG.json --output FRESH_DIRECTORY` requires
Slurm and NumPy; one CPU suffices. It writes RESULT.json, ALIGNMENTS.json,
ATOM_MOTION.tsv, RESIDUES.json, WINDOWS.json, MECHANICS.json and REPORT.md.
Missing cells remain explicit; source/map inconsistencies fail closed. Root owns
submission and any later feasibility/reference stage. No chemistry is launched.
