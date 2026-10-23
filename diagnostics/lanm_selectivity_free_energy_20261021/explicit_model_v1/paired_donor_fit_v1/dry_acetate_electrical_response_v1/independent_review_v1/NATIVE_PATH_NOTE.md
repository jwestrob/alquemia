# Original donor field route: independent source check

The original installed analyze binary supports the proposed static field-energy diagnostic without a build. This note precedes final adapter/preparation review; it is not an execution release. No molecular program was invoked.

`prmkey.f:402` parses `EXTERNAL-FIELD` in MV/cm and divides by `elefield`. `induce.f:215–225` adds that field to **both** permanent-field arrays before the native mutual solves. Analysis uses `epolar3a`, whose pair sum employs the energy-channel permanent scaling; its external-field addition is `−K/2 Σ uind·f` (`epolar3.f:553–565`). `exfield3` adds `−K Σ(q r + dipole_perm)·f` using literal laboratory coordinates (`exfield.f:346–365`). Native analysis therefore includes field coupling exactly once. Binary symbol inspection found the corresponding parser, analysis, field and induction routines.

With the original symmetric mutual operator A, native d/p fields Ed/Ep, and a uniform field vector f repeated over induced sites,

```
u = A⁻¹(Ed+f),       v = A⁻¹(Ep+f)
U_pol = −K/2 (Ed+f)ᵀ A⁻¹(Ep+f)
−∂U_total/∂f = K [mu_permanent + ½ Σ(u+v)] .
```

The donor graph does **not** justify equal u/v. Its methyl group (types401/404) and carboxyl group (402/403) are separate polarization groups. Their C–C bond receives direct scaling1 but energy-channel polarization scaling0. This follows from the selected key, original prm and `dfield0a` scaling at `induce.f:785–825,922–925`. Energy odd differences correctly recover the conjugate dipole without requiring equal fields. Neither the usual u-only saved induced dipoles nor COM-based printed M should replace that derivative without a separate proof.

Use the already qualified native-energy conversion: for field h in atomic units, `F_MV/cm = h * HARTREE * ELEFIELD / (K * BOHR)`; with the pinned constants this is `5142.206758812493*h`. Then `−[E(+h)−E(−h)]/(2*h*HARTREE)` is in e·Bohr, and multiplication by `BOHR*DEBYE` gives Debye. Curvature is `−[E(+h)+E(−h)−2E0]/(h²*HARTREE)` in a0³. The fixed-geometry original model is quadratic in field at exact mutual convergence; step disagreement tests numerical consistency, not a model hyperpolarizability signal.

The exact selected monopoles sum to +2e for both metals. Native field energy uses the fixed laboratory origin, not the mass center. No COM correction applies to an energy derivative. A fixed origin translation changes each dipole by the same charge-dependent constant, which cancels in same-charge metal differences and dμ/dq. Recentring each geometry or rotating the field with q would instead change the observable. Keep the declared field direction fixed in the inherited frame.

At each fixed geometry the bonded and vdW terms are field independent; their odd energy difference cancels. Thus dμ/dq measures the mixed derivative `−∂²E/(∂F ∂q)`, an electrical response along the chosen motion. It does not identify which field-independent energy correction would repair the zero-field force, or uniquely separate permanent multipoles, induction and charge transfer in a quantum reference. Centroid rotation also changes metal–donor distances. Agreement in one projected response cannot establish the full vector, tensor, local-field response or physical accuracy.

Native field terms increment the Atomic Multipoles and Polarization interaction counters once per site; these eight-site inputs therefore change those counts from13 to21. Other component counts and the eight component labels remain unchanged. Field parsing must account for this source-supported inventory change.

For each printed energy with half-unit ε, the absolute odd-moment printing bound is `ε/(h*HARTREE)` e·Bohr and curvature bound is `4ε/(h²*HARTREE)`. Propagate absolute coefficients through both field and geometry differences, and again through Dy-minus-La. Printed rounding bounds are not induction-error bounds. Zero energies must come from exactly the unchanged selected center/rotation bytes; the odd estimator does not itself require E0. No new tolerance is introduced by this note.
