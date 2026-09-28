# Smallest coherent coupled scaffold-response candidate

**Implement one explicit finite additive model on repaired Hans8DQ2, then test one real peptide torsion.** The parent parameters and maps exist; the missing implementation is a correctly defined boundary ledger plus physical force assembly. No Hessian, fitted spring or whole-protein optimization is needed to learn whether scaffold motion is energetically meaningful.

This is a design-only audit. No endpoint, structure, parameter or active job changed.

## Fixed-length caps determine the boundary policy

Current caps are `l = rQ + d*(rM−rQ)/|rM−rQ|`, with fixed C–H1.09Å/N–H1.01Å. Moving the omitted MM atom radially along the cut leaves the cap unchanged: `J_M * u = 0`. Therefore electronic cap energy supplies **no restoring stiffness for that real cut-bond length**. Keep the four original native cut-bond stretch terms. Deleting them because a generic QM/MM default deletes link bonds would leave a real mechanical mode unsupported.

For angles/torsions use source-support substitution, not merely a count of MM atoms. A native cross term is represented by the capped electronic region only when every exterior atom in the term can be replaced by its unique actual cap and all remaining atoms are real QM. Angles and dihedrals under that replacement have the same physical angular coordinate, because caps follow the real cut direction. Omit these native terms to avoid applying a second approximate local angular potential. This is a declared boundary approximation, not proof that C–H cap and original C–C/C–N potentials are identical.

The actual Hans topology gives:

| Cross native term | Keep | Omit as represented by caps |
|---|---:|---:|
| Cut-bond stretch |4|0|
| Angles |12|8|
| Proper/improper torsion entries |21|11|
| CMAP |4|0|

The torsion distinction matters: ten omitted entries have one MM atom, but the eleventh is `CA65–C65–N66–CA66`. **Both** exterior CA atoms have actual caps, so this torsion is represented entirely in the doubly capped supporting amide. Blindly retaining all two-MM torsions would duplicate it. The retained cross angles/torsions/CMAP include additional actual exterior atoms absent from the electronic region, such as `N65–CA65–C65` and the upstream backbone CMAP. Wholly-QM native bonded terms remain omitted, and wholly-MM terms remain native.

Verify the exact term IDs against the graph before building the new ledger; the counts above are from the actual1,887-atom XML/source selection. This policy intentionally differs from both the earlier all-cross1H4I ledger and unexamined native ORCA link defaults.

## Complete declared finite energy

For fixed source occupancy and charges:

`E_M(x) = E_PBE0-D4,embedded,M(Q(x),caps(x); qtilde_MM,R_MM(x))`
`         + U_native_MM_bonded(x)`
`         + U_native_cross_bonded,selected(x)`
`         + U_LJ(MM,MM;x) + U_Coulomb(qtilde_MM,qtilde_MM;x)`
`         + U_LJ(realQM,MM;x)`.

Use current literal lcecp1TZ/ECP46La/ECP55Dy, light-atom def2TZVP, matching AuxJ and unchanged source-specific charge states. `DoEQ false` excludes external–external Coulomb from the quantum component; add that term exactly once. No classical QM–MM Coulomb, QM–QM LJ/bonded, cap LJ, C4 induction or extra D4/gCP. Keep native MM–MM/realQM–MM LJ exceptions and original1–4 scales; regenerate MM exception charge products from actual redistributed charges. The zero-field-charge boundary CA atoms stay physical with bonded/LJ terms. Nonperiodic finite conventions are fixed; no hidden continuum or tail correction.

Physical derivative = real-QM electronic gradient + both-anchor cap chain-rule gradients + `.pcgrad` distributed onto the actual MM source positions + analytic derivatives of all declared classical terms. The shifted charges are fixed numbers attached to source atoms, so this is movement of permanent charges, not polarizable MM response. On any collective motion regenerate every cap coordinate/Jacobian and every moved pointcharge coordinate. Never use an intrinsic isolated QM gradient or leave the field fixed while moving its source atoms.

This is a mathematically explicit model whose derivative can be tested. Its cap representation and short-range QM/MM interactions remain approximations; numerical force agreement would not prove biological accuracy.

## Outer ions and available parameters

The complete physical Hans model has targetEF3 plus fixed LaIII EF1/EF2 and NaI EF4. Keep identities and source positions explicit. Spectators may be held fixed for the first torsion test, but their forces/interactions with moving protein atoms cannot be omitted from the energy.

Candidate pure12–6 values are available in installed `frcmod.ions234lm_126_tip3p`: LaIII Rmin/2=1.718Å, epsilon=.15060822kcal/mol; DyIII1.609/.08389240 (IOD sets). `frcmod.ions1lm_126_tip3p` supplies NaI1.475/.03171494 (HFE set). Real protein/water QM atoms already have native LJ. These are actual parameter candidates, **not a validated frozen-f QM/MM cross fit**. Use the same declared mixing rule/native exceptions throughout; no arbitrary zero or reconstructed C4-stripped parameters. A finite model can be tested with these frozen candidates while keeping transferability unqualified.

Only electronic embedding represents polarization of the QM region. Adding12–6–4 induction would risk counting that response twice. Source spectators represented as fixed charges/LJ still cannot reproduce their own responsive electron clouds.

## One source-derived covalent-preserving mode

Use the repaired Hans source and the actual backbone `N83–CA83` bond. Removing that bond in the native source graph gives a genuine bridge. The CA83-side component contains814 protein atoms spanning residues83–133:176 real-QM and638 MM atoms. Rotate that whole connected component about the original N83→CA83 axis, carrying every attached H; leave upstream protein, actual waters and all four metals fixed. This is a physical backbone-phi motion coupled to the scaffold, not independent fragment translation. The rotation preserves covalent bond lengths and local bond angles; it changes the backbone torsion and actual contacts. It is intentionally a small collective perturbation, not a proposed equilibrium folded state.

At the repaired origin, use the complete physical gradient to project onto this mode. Evaluate both metals at ±0.05°: four new electronic endpoints if matching origin energies/gradients already exist, otherwise two origins additionally. Compare central energy differences with the analytic directional derivative of the **whole defined energy**. Pin an acceptance rule before those calls, based on the actually qualified reference numerical precision; do not invent a pass tolerance from the observed chemical signal. If that comparison is numerically ambiguous, the declared diagnostic refinement is ±0.025°, at most four further endpoints, not a geometry search.

Before quantum calls, evaluate the inexpensive classical pieces and source Jacobians on the same mode, checking preserved covalent geometry, finite contacts and the expected nonzero MM/cap force paths. Test the cut-bond radial null property of the cap algebra separately and verify the retained native bond supplies its known stiffness; that algebra/classical test needs no quantum evaluation. Do not mistake the previous Asp85 sidechain rotation (zero exterior motion) for this test.

## Hard limits and next decision

The first implementation is blocked on **building and checking this particular selected-term ledger and force projection**, not missing full-protein mechanics. The source-native H repair supplies sane C-bound H geometry; exchangeable-H strain remains and the source is not an equilibrium state, so retain the affine force. Do not infer thermal covariance or entropy.

The present finite environment has crystal waters only. It is not a bulk-solvent relaxation Hamiltonian. A successful small torsion comparison can establish conditional coupled response and identify whether the exterior changes the local metal contrast; it cannot authorize unconstrained whole-protein relaxation, affinity or state populations. Moving toward that goal requires a separately coherent solvent model and short-range metal-response qualification, not a larger allocation of the same dry calculation.

Recommendation: **one selected-ledger implementation plus this finite two-metal phi response**, after current repaired electronic outcomes are interpreted. It directly tests the missing coupling while keeping chemistry fixed. If this explicit model cannot pass its force accounting, repair that concrete defect; if it passes but adds no useful metal-dependent response, close that candidate rather than escalating immediately to a Hessian or global optimization.
