# Frozen common8 solvent-guided line probes

Execution agreement: `diagnostics/nikasha_parallel_pilots_20260922/PLAN.md`.
Exact source panel and existing source/model/state pins: its `INPUTS.json`.
Archive replay precedes this design; no new geometry energy has been inspected.

## Question and evidence

On204complete fold pools, native MACE and composite minima differ for93Ca and
29La rows; composite regret exceeds the existing0.1kcal/mol selection scale in
90Ca/24La rows. This justifies testing proposals around the geometry preferred by
the actual composite. It does not demonstrate an accuracy improvement by itself.
The current pipeline already solvent-selects the existing common pool. This pilot
asks whether using that preference to propose additional coordinates finds useful
candidates missed by vacuum-only optimization; it does not claim selection is new.

## Exact geometry rule

For each of the eight fixed sources and each of Ca/La, use the archived current
operational composite winner as the center in the existing four-angle kinematics.
The two metals use the same source graph, active modes, chemical composition and
origin. All unselected physical modes remain fixed.

1. If native and composite operational winners differ, set direction to
   `q_composite − q_native` using their actual physical angular coordinates.
2. If they agree, use `q_composite − q_origin`. This tests extension/contraction
   along the existing accommodation path, not a metal translation or new mode.
3. Normalize by the largest absolute selected-angle component. A direction of
   maximum magnitude≤1e-12rad is explicitly unavailable; no alternate axis.
4. Independently probe each sign at the first admissible maximum-angle step in
   the fixed sequence **0.10,0.05,0.025,0.0125rad**. Backoff uses geometry only;
   no energy is evaluated at rejected steps or used to choose a distance.
5. Retain the existing±0.8rad bounds, maximum0.8Å source-heavy-atom displacement,
   stated numerical tolerances, source bond/cap/fixed-atom/overlap checks. If no
   listed step is admissible, record an unavailable probe with every reason.
   The shared finite adapter's conservative1e-10Å displacement-roundoff allowance
   is also applied during backoff; no point relies on the older1e-7Å allowance.

This direction is a finite secant construction from actual candidate preferences,
not an analytic composite gradient or a claimed minimum. Centers and directions
must replay their archived XYZ within the already established1e-12Å tolerance.
Unsupported coordinate mappings remain explicit; no geometry or chemistry rescue.

## Finite scoring and interpretation

At most32new proposed geometries (four per source). Deduplicate identical
coordinates before execution, retain every proposed ID/reason, and give both
metals the complete common candidate union. Maximum new cross-scoring is
**64native float64 MACE calls and128native GFN2 singlepoints**. No optimization,
new DFT, CPCM, new folds, numerical DFT gradient or surrogate fitting.

Native MACE + native GFN2(ALPB−vacuum), fixed300K/state, unchanged convergence and
qualifiedMaxIter500. Existing exact five-geometry pool and source origin energies
are reused. Root's finite-candidate adapter owns paired scoring and preserves
failures; this branch emits only pinned physical candidates.

Select actual composite row minima with the existing0.1kcal/mol origin-retention
rule; report mathematical minima separately. Do not choose a geometry using class
labels. Report both metal works and native/solvent contributions, ΔR, every
unavailable probe/case, old released/adaptive-band transfer and complete common8
denominator. No new threshold or canonical reference is fitted to this panel.
Changes in energy alone are not an accuracy claim. Broader evaluation remains a
separate parent decision after this contained result.
