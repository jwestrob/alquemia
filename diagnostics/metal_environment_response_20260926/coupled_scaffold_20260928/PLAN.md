# Repaired Hans selected scaffold ledger and finite checks

Declared before new energies. Exact repaired Hans8DQ2 full source, native1,887atom ff19SB/TIP3P parent, target EF3La/Dy, actual frozen La/La/Na spectators. Complete physical model1,891atoms;191realQM including target,1,700MM. Four fixed-length source caps exist only in electronic mapping.

Classical energy = retained native bonded (all MM, four cross radial bonds,12crossangles,21crosstorsionentries,4crossCMAP) + MM–MM LJ + shifted-charge MM–MM Coulomb + realQM–MM LJ. Omit whollyQM terms and cap-represented8crossangles/11crosstorsions, including doubly capped supporting-amide torsion. No capFFparticles or classicalQM–MM Coulomb/C4/additional dispersion. External electronic energy with mapped caps/field is not evaluated here and stays explicitly unavailable. This is classical-component/geometry qualification, not complete hybrid force qualification.

Literal pure12–6 Amber IOD/TIP3P La/Dy and HFE/TIP3P Na parameters as declared in COUPLED_SCAFFOLD_DESIGN. Lorentz–Berthelot mixing, native exceptions. The four zero-field-charge MM1CA particles remain mechanical. Use existing source redistribution, no neutralization; update MM exception products using native exclusions/1–4scales. Preserve actual source ion charges and coordinates. Parameter transfer and bulk-solvent adequacy remain unqualified.

Use the source graph N83–CA83 bridge; rotate all814downstreamprotein atoms about it, carrying176realQM/638MM atoms. Three actual configurations only:0,±0.05degrees. All water/metal/upstreampositions fixed. Regenerate core/caps/Jacobians/field coordinates. Two targetelements×three configurations=6classical native evaluations, no optimization or QM calls.

Checks fixed before evaluation: source bond lengths and angles unchanged within1e−10Å/radian; cap Jacobian coordinate finite differences(step1e−5Å) maxresidual1e−8; radial cap J_MM*u norm<1e−12; centralphi energy derivative vs analytic projectedgradient residual≤0.02kcal/mol/radian+2e−5*abs(analytic), each component/metal. Finitephi cap-coordinate derivative residual≤1e−5Å/radian, consistent with O(h²) at±0.05°. Exact native cut-bond radial curvature equals forceconstant; direct±1e−4Å native harmonic expression residual≤1e−6kcal/mol/Å²+1e−8relative. Geometry/finite contacts reported, no electronic acceptance claim.

OpenMM Reference,1sharedCPU,mem0, noGPU. Keep actual failures and costs. Root owns any new electronicforce checks or adoption after inspecting this ledger. No shared runner or production changes.
