# Dry mono-acetate: electrical derivative of the existing rotation

## Question and finite scope

Determine whether the frozen donor model also misrepresents the dipole and its change along the consumed dry mono-acetate rotation. Both model and reference passed their declared numerical derivative checks, yet their rotational energy slopes have opposite signs. A field-independent short-range energy correction cannot change a dipole or this mixed electrical/geometrical derivative. A disagreement here would therefore establish an additional electrical defect; agreement would leave an orientation-dependent energy defect plausible without identifying its unique cause.

Exactly **40 original-model `analyze system.xyz E` calls**: the existing La/Dy centers and eight existing rotation endpoints, each at four signed fields. Reuse all ten zero-field model energies and all ten full-precision quantum dipoles. No new quantum calculation, geometry, instrumentation, executable, parameter point, force call, or fit belongs to this design. The conditional quantum extension recorded in the earlier inventory is outside this frozen scope.

`SOURCE_INVENTORY.json`, bound by `READY.json`, gives the ten literal row IDs, raw reference/property files, original model XYZ/key/parameter files, and zero-energy authorities. Rotation coordinates are the original normalized displacement coordinate q = ±0.010 and ±0.005 Å. Preserve the original rotation length, axis, centroid, coordinates, atom order and metal at the laboratory origin. These are consumed diagnostic fixtures; no reserved outcome is opened.

## Frozen direction, field convention and energy

Use h = **0.001 au**, half = **0.0005 au**, and signed fields −h, −half, +half, +h along one fixed laboratory unit vector:

```
t = (0.7562812564509392, -0.15594410634600103, -0.6353897204369424)
```

This vector is determined by the original central geometry: normalize n × (O3 − ligand centroid), with n the original rotation axis. O3 denotes zero-based atom 3 in the inventory construction. Keep this same vector for both metals and every endpoint; do not rotate it with the ligand or select it using an observed error.

The exact original `analyze` SHA256 is `f67420e79f8c1793bf491569f77ab4354219926fe5b4c49a2311e64f217f7a28`, identical to the independently admitted 19-call water-field experiment. `FIELD_AUTHORITY.json` binds that actual result/review and the original native field, induction, energy and unit sources. Its zero-identity and bare-ion ±h controls tested field sign, scale and energy factor. Water moment equality does not qualify donor moment equality.

The native `EXTERNAL-FIELD` vector is in MV/cm. For a field s in au render `s * 5142.206758812493 * t`; this conversion uses Hartree = 627.509474063 kcal/mol, Bohr = 0.529177210544 Å, native electric = 332.0637133, dielectric = 1 and native elefield = 1439.96455. Record the actual rendered vector and source constants. Preserve all original active key tokens, including donor parameters, induction settings and frame definitions; only the static field is added, with path relocation documented separately.

Native total energy already contains permanent and induced external-field contributions. Add no further −μ·F term. The field enters both induction right-hand sides. For the symmetric linear operator A, the polarization energy is −K/2 (Ed+f)ᵀ A⁻¹(Ep+f); its field derivative contains half the sum of the two induced solutions. Donor direct/polarization exclusions differ, so neither Ed=Ep nor u=v is assumed. Energy differences measure the conjugate dipole without using `SAVE-UINDUCE` or a water-only moment parser.

## Common origin and reused quantum evidence

All complexes have charge +2. Model field energy uses the literal laboratory coordinates, with metal at zero. Convert the existing quantum COM-origin total dipoles to this same laboratory origin before any comparison:

```
O = (sum_A Z_eff,A r_A - mu_nuclear(O)) / 42
mu_total(0) = mu_total(O) + 2 O
```

Use property-file coordinates and dipoles in atomic units, effective nuclear charges [11,6,6,8,8,1,1,1], 40 electrons, and full-precision property tokens. ECP cores are not added back. Inventory checks already join these values to the native printed origin and mass-COM reconstruction. Printed-origin-only values have a conservative 1.0005e−6 ea0 bound per Cartesian component; project and propagate this separately from full-precision property arithmetic. A moving COM must not be used for the rotation derivative. A common fixed origin shift cancels from the same-charge derivative and metal difference.

The centers and rotation SPs contain dipoles and population analyses, but no quantum static polarizability. The saved projected quantum dμ/dq values at q=.010/.005 Å are La −.718669/−.717914, Dy −.789076/−.788829, and Dy−La −.070407/−.070915 ea0/Å. These are descriptive property differences, not newly qualified quantum susceptibility or a physical accuracy bound. Population charges remain partition-dependent.

## Estimators and prospective numerical checks

Convert model total energies to Hartree. For every geometry and s∈{h,half}, retain

```
mu_E(s)  = -[E(+s)-E(-s)]/(2s)                  # ea0, projection on t
chi_E(s) = -[E(+s)+E(-s)-2E(0)]/s²             # a0³, projection tᵀχt
D(q;s)   = [mu_E(+q;s)-mu_E(-q;s)]/(2q)        # ea0/Å
```

Report all ten moments, their model-minus-quantum differences, and D for both q steps, both field steps, both metals and Dy−La. The mixed derivative obeys ∂F_q/∂F_t = ∂μ_t/∂q with F_q=−∂E/∂q. It measures the field dependence of the generalized rotational force, not a quantum energy-component decomposition.

All 40 slots and ten compatible zeros must be admitted, finite and exactly joined. Reuse the original native energy/induction admission without weakening it. Then apply these frozen numerical resolutions:

- Every geometry: |mu_E(h)−mu_E(half)| ≤ **1e−4 ea0** and |chi_E(h)−chi_E(half)| ≤ **0.5 a0³**.
- Every matched La/Dy geometry pair: the direct two-step discrepancy of chi_Dy−chi_La is ≤ **1 a0³**. Retain the analogous moment discrepancy and its implied 2e−4 ea0 envelope.
- D field-step differences retain the implied endpoint envelope **1e−4/|q| ea0/Å per metal**, twice that for Dy−La, and the actual propagated print bounds. These follow from the individual moment checks and are not independent physical acceptance criteria.
- Compare the two geometrical q steps descriptively. Add no new q-step, dipole-accuracy, torque-accuracy, or model-minus-quantum threshold. Do not use the quantum moment as a numerical consistency gate for the model.

The 1e−4/.5/1 resolutions come from the prior field design, prospectively applied here to energy-only step consistency. They are not variances, confidence intervals or reference error bars. At fixed geometry the original model is linear in induced moments and quadratic in a static field; field-step disagreement tests numerical/input/solver consistency. It does not by itself assign a cause. No failed check is replaced by extrapolation, averaging, new steps or relaxed thresholds.

Retain actual raw decimal tokens and half-last-unit bounds b. Use b_mu=(b_plus+b_minus)/(2|s|) and b_chi=(b_plus+b_minus+2b_zero)/s², then absolute-coefficient propagation through geometry/metal differences. With eight-decimal kcal/mol energies, half-step bounds are approximately 1.59e−8 ea0 for mu and 1.275e−4 a0³ for chi; these bound printing only, not induction or method error.

## Interpretation, staging and cost

If moment or D errors accompany the torque error, a field-independent short-range correction alone cannot repair all observed defects. If they agree in this fixed projection, the result weakens that specific electrical explanation while leaving other field modes, anisotropy and non-electrical energy errors unresolved. This panel does not separate erroneous local driving fields from erroneous response stiffness: no quantum field curvature is available. Model chi is a consistency diagnostic, not a quantum susceptibility comparison. It identifies neither unique charge transfer nor a correspondence between a quantum contribution and an AMOEBA energy component.

Preserve the original PBE0-D4, frozen-f ECP46/55, effective singlet/physical La-singlet and Dy-sextet distinction, basis/Aux-J and final-integration limitations. The central EnGrad and displaced SP outputs retain Energy Check stopping with failed printed density/DIIS/rotation criteria and nonidentical integration counts. Their admitted directional derivatives do not establish method accuracy or a physical electronic state. D4 is field-independent here and supplies no dipole; its rotational energy remains separately recorded.

Use the two central quartets as the first **8 calls within 40**; admit their source/input/energy/induction/step checks before the remaining 32 rotation calls. Any native or numerical failure stops continuation and retains every unavailable slot. Root owns release and execution. Run original serial, one-thread native programs with scheduler defaults, no retry or extra controls. The preceding 16 energy calls took about 7.43 summed native seconds / 8.89 job seconds; linear scaling estimates about 18.6 native seconds / 22.2 single-CPU job seconds for 40, not a promise or cap. Record actual cost. Preparation and this design make zero native calls.
