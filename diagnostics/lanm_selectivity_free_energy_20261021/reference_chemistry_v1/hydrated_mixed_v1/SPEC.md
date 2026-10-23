# Hydrated donor substitution: physical preparation only

## Question and declared scope

Does a model's metal-dependent donor-versus-water interaction survive first-shell
many-body polarization and water competition? The isolated dimers cannot answer
this. Prepare a small matched explicit first-shell panel before selecting its
electronic or force-field scoring protocol. No electronic calls, optimization,
sampling, fit or submission is part of this preparation.

Use the frozen initial panel's exact eight-water control. Construct three
substitution families from the supplied acetate/amide internal geometries:

1. M(W)8 + acetate -> M(acetate-mono)(W)7 + W.
2. M(W)8 + acetate -> M(acetate-bi)(W)6 + 2W.
3. M(W)8 + acetamide -> M(acetamide)(W)7 + W.

M is La3+ or Dy3+. Acetate complex charge+2; water/amide complex charge+3.
No proton changes. The physical spin metadata remain La1/Dy6. Each family has
central and +/-0.15Angstrom rigid ligand translations along its radial axis,
with all retained water coordinates fixed. Exactly18 new mixed configurations
plus2 reused aquo controls; matching isolated free fragments are included only
for reaction bookkeeping. Metal-difference cycles cancel those free fragments.

## Deterministic construction

Monodentate substitutions remove water0. Map the original coordinating oxygen
to the removed oxygen's radial axis at2.50Angstrom, rotating/translating the
entire ligand rigidly. Bidentate acetate removes adjacent upper-ring waters0/1;
map its O--O midpoint onto the bisector and retain the exact O--O separation,
placing both oxygens at approximately2.50Angstrom. All internal ligand distances
remain unchanged under orthonormal transformations.

For each family consider eight evenly spaced orientations around the radial
axis. Choose the one maximizing the minimum ligand-to-retained-water distance
over all three declared radial positions; exclude the metal. Record all eight
scores and the selected angle. This purely geometric, metal-independent choice
avoids an arbitrary initial overlap; it does not use energy, reference forces or
the desired preference. Both metals receive exactly the same final positions.

Check finite coordinates, preserved ligand internal distances, exact retained
water coordinates, proper rotation determinant, covalent inventory, charge and
electron parity, and balanced atom/charge stoichiometry for all three reactions.
Report the shortest interfragment H-containing and heavy-heavy contacts, including
any contact less than1.2 or1.8Angstrom respectively; reject such fixtures rather
than silently changing the design. These contact limits admit a stress fixture;
they do not prove a minimum or realistic equilibrium hydration.

## Interpretation and next gate

The intended contrast is E(Dy mixed)-E(Dy aquo8)-E(La mixed)+E(La aquo8).
Negative means greater fixed-state Dy stabilization relative to water, under
the selected model only. Report ligand-radial and full forces as well as energy.
This panel contains one representative donor per cluster, not a protein pocket.
Successful preparation supplies inputs, not qualification of either model.
Scoring requires root's review of the first native pilot and its remaining
electronic-state/force limitations. Actual sampled mixed-donor environments and
outer-shell convergence remain necessary later.

The free-water bookkeeping geometry is copied from the constructed aquo8
control, whose water geometry differs from the published dimer water. Free
ligands retain their respective published monodentate/bidentate/amide internal
geometries. Their energies cancel within each paired metal contrast, but the
different acetate internal geometries must remain explicit across donor modes.
These conditional Ln3+ reactions exclude buffer complexation and protonation.
In particular, MES is not assumed inert; matching assay conditions needs an
independent paired La/Dy speciation treatment.
