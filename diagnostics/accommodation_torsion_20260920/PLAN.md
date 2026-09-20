# Fixed whole-carboxylate profiles

Parent-approved scope is [TORSION_PLAN.md](../accommodation_goal_20260920/TORSION_PLAN.md).
This implementation clarification is frozen before any new energy evaluation.

- Exact archived original sources: 1H4I, 4MAE, and PLM sample0
  PQQSEQ_83440678cbbd658047c9 / PQQSEQ_07ab500e3df76b30d71c.
  PLM biological labels remain unknown; extremes were previously inspected.
- Anchor Glu terminal CG–CD torsion in each case; extra Asp terminal CB–CG
  torsion in all except 1H4I. Role selectors come from the source core manifest.
- Each path uses −0.4, −0.2, 0, +0.2, +0.4 rad. Common origins shared between
  paths: 32 geometries / 64 metal endpoints. Complete physical carboxylate
  rotations; source scaffold, PQQ, water inventory, protonation and charges fixed.
- Exact original context origins and mappings are reused for the two controls.
  Existing complete_expansion prepares PLM contexts from original archived H,
  without normalization, reprotonation or geometry minimization.
- Native float64 OMOL100M scalar energies: 64 evaluations; native ORCA GFN2
  vacuum/ALPB(Water), 300 K electronic smearing/native mixer, published primary
  convergence recipe: 128 singlepoints. This is the established composite;
  analytic pilot tight gradients are a separately quantified reference.
- Native r2SCAN-3c/CPCM(Water)/DefGrid3: 16 frozen validation tasks (extra-Asp
  ±0.2 rad, three cases × two metals; four PLM origins). Exact existing 4MAE
  context DFT origins reused. No finite-difference DFT gradients or optimizer.
- Full mapping/bond/Jacobian check at each point. A newly introduced atom-pair
  distance below 0.55 Å (either H) or 1.0 Å (both heavy) is an explicit severe
  geometry failure; unchanged source contacts remain recorded. This is only a
  catastrophic-overlap screen, not a claim of relaxed physical plausibility.
- Primary control-origin repeat tolerances: 0.01 kcal/mol per OMOL endpoint,
  0.10 per GFN2 solvent transfer, 0.20 per composite metal contrast. Freeze these
  before results; report exact differences. Native DFT finite-work errors are
  reported individually with signs, no fitted scale/weights or label-selected
  favorable angle. Boundary-descending/unstable paths remain visible.
- No probabilities, entropy, harmonic relaxation scalar, broad accuracy or
  production promotion. All failed points remain in the denominator. Finite
  manifests define the experiment; there is no arbitrary aggregate budget cap.

The intervention tests whether this specific affordable solvent-consistent
surface describes real physical donor displacement energies. It cannot prove
that absent minima do not exist along other coordinates or coordination basins.
