# Fixed-field DFT/MACE short-component screen

Declared after the canonical direct scorer failed, before evaluating this new
combination. Covered by Jacob's active MACE goal and discretionary pilot approval.
Production baseline and previous experimental records remain unchanged.

## Question and smallest informative calculation

Can the exact learned local MACE contribution account for the residual
representation dependence of the existing DFT plus permanent-field model?
The previous density study left about 10 kcal/mol dependence when Asp303 moved
across the core boundary. The learned global charge response also failed.
This screen combines already computed, compatible components; it requires
zero new DFT, MACE, charge extraction or solvent calls. Do not rerun a panel.

Use only the consumed 1H4I qm33/qm36 La/Ca states with their exact original
coordinates, caps, assembly, protonation, zero waters and native vacuum
r2SCAN-3c settings. Use both saved analytic MACE checkpoints: medium is primary,
large a separate sensitivity result. Use the four uniformly extracted native
CHELPG charge fits and the existing fixed ff19SB environment. Exact saved-density
coupling is a diagnostic comparator, never a per-state charge-scheme choice.

## Defined candidate and accounting

For partition P and metal M, define

```
V_P,M = E_DFT,vac(core_P,M)
      + C_CHELPG(core_P,M; environment_P)
      + E_MACE,short(full_M) - E_MACE,short(core_P,M)
R_P   = V_P,Ca - V_P,La
D     = R_qm36 - R_qm33
```

`short` means the installed model's `interaction_energy` scalar, after local
readout/scale-shift and before global charge updates. It is not total MACE
energy, total forces, or the field-dependent `electron_energy`. Isolated atom
offsets cancel in the full-minus-core difference between metals. Verify that
both full inputs are identical between partitions and match the archived source.

The DFT core contains all internal core terms; its own MACE short approximation
is subtracted before the full MACE local contribution is added. Explicit direct
Coulomb couples the quantum charge representation only to environment charges;
no environment charge is owned by a source QM atom. Environment-only Coulomb
terms are constant between metals within a partition and cancel from R. They
need not be equal between partitions to cancel from each separate R. The full
MACE short contrast is identical in both partitions and cancels from D.

The original boundary charge redistribution and cap positions are preserved,
including their documented limitations. No new redistribution or neutralization.
CHELPG sums retain their printed precision; require the existing 1e-5 e tolerance.
Compare against exact density using the same recorded atomic-unit Coulomb
convention. Convert DFT Hartree and MACE eV once; saved coupling is kcal/mol.

This is a new approximate energy model. The trained short/field split is not a
unique physical decomposition: its transfer to external charges is unvalidated,
and point-charge replacement may mismatch short-range charge penetration. A
successful partition check alone cannot validate the split or affinity. No
solvent, CPCM, reaction field, self-consistent polarization, mechanical response,
entropy, aquo reference, absolute decision or existing calibration is included.
Solvent is **unavailable**, not assumed numerically zero for a solution score.

## Frozen checks and outputs

- Verify actual execution receipts and hashes, exact core/full mapping, source
  state equality across checkpoints, endpoint charge closure and immutable inputs.
- Recompute direct CHELPG coupling from the saved charges and actual environment;
  require agreement with its archived value within 0.01 kcal/mol. No new fit.
- Retain DFT, coupling and short components, raw R, partition D, and the exact
  density diagnostic. No fitted weights or selection of a favorable checkpoint.
- Audit the algebraic isolated reduction using an actual archived core as both
  full and core with the empty environment. Added terms must vanish within
  1e-8 kcal/mol. This is an algebra test, not an executed solvent identity test.
- Check saved rigid transformations of each short component against the existing
  0.01 kcal/mol tolerance. Coupling rotates analytically with its source charges;
  a newly fitted rotated CHELPG distribution is not tested here.
- Keep the previous 2 kcal/mol absolute partition target as the screen criterion.
  Report all results even if it fails. A pass permits designing a separate
  common-boundary solvent model; it does not authorize calling this a solution
  affinity score. A failure ends this frozen screen, with no parameter tuning.

Report both checkpoints and the approximation's limits. All input cases have
been consumed; no biological classification test is performed. Record local
analysis wall/CPU time and the original cached-job costs separately. There is
no project time/CPU stopping budget and no Slurm submission for this screen.
