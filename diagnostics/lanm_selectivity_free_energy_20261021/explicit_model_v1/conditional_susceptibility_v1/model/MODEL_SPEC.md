# Original-baseline model conditional susceptibility

This implements only the model counterpart of the frozen shared reference SPEC (`reference_chemistry_v1/conditional_susceptibility_v1/SPEC.md`). Original 17676 six parameters, original CPU analyze binary, no pair-response keywords, EXCHPOL, new geometry, gradients, fit or production change. The guarded prototype is untouched.

## Finite sequence

One Slurm task/CPU, all native and numerical thread settings one. Exactly 19 maximum `analyze system.xyz EM` programs, sequential, fresh immutable runtime directories:

1. Exact saved W9 far Dy trimer with an explicit zero external field.
2. Exact saved bare La ion at z=+0.001 au.
3. The same bare ion at z=−0.001 au.
4. Only after all three controls pass: W9 far M+water-1+water-3 and its exact M+water-1 dimer, La/Dy, each z=−0.001, −0.0005, +0.0005, +0.001 au. All sixteen slots are retained.

First native/parse/control failure stops further starts. All successfully parsed output and every unattempted slot remain. The final finite-step summary applies every shared numerical gate without fitting or selecting a field point. A failed consistency gate retains the entire actual curve as numerically unqualified; it is not automatically attributed to SCF error rather than nonlinearity. Root owns release/submission. No calls occur during preparation or saved controls.

## Exact reuse and field energy

The five baseline physical inputs and zero-field moments come from the admitted model78 and trimer results. Preparation replays these five raw energies/site dipoles, native receipts and exact parameter/input joins. The four quantum real-atom masks/coordinates must match Singer's immutable shared inventory, with no rotations or fragment recentering. Old results and sources stay unchanged.

Each copied key gains exactly one `EXTERNAL-FIELD 0 0 F` record. F=field_au×hartree×elefield/(electric×bohr), using pinned native constants: bohr=0.529177210544 Å, hartree=627.509474063 kcal/mol, electric=332.0637133, elefield=1439.96455 and Debye/eÅ=4.803206802. Conversion is 5142.206758812493 MV/cm per au, preserving native energy/moment consistency. No extra −mu*F correction is added. Unchanged SAVE-UINDUCE/POLAR-PRINT report both native induction solves; every residual must be finite, nonnegative and ≤1e−8.

Original complete-water groups have equal direct/polarization exclusion scaling; this scoped equality permits using native uind for the energy-derivative dipole. It does not generalize to arbitrary protein graphs or enabled pair-response kernels.

## Precision and moments

Six-decimal site induced dipoles are Debye. Full laboratory moments are reconstructed as the fixed permanent moment plus their sum. The permanent contribution uses only the already admitted local-frame signatures: q*r, O Bisector and H Z-then-X rotations, with native bohr conversion for input dipoles; zero metal permanent dipole. Exact native kmpole/rotpole/moments sources are pinned. Preparation revalidates the frame/atom/type records and cross-checks all five saved outputs against printed total M after the explicit +Debye*Q*COM correction. Runtime repeats that cross-check.

The native total M has a ±0.0005 Debye/component print bound. Reconstructed total moments use N×0.5e−6 Debye/component plus a prospective 2e−9 Debye arithmetic allowance. The low-precision M is retained as a check; it does not set the finer 1e−4 ea0 odd-slope precision. Energy rounding is ±5e−9 kcal/mol per value. Bounds propagate with absolute cycle coefficients, field divisions and native unit conversion; solver residuals alone are not asserted as observable error bounds.

## Initial controls and final gates

Before matched fields, require zero-field identity |dE|≤1e−8 kcal/mol and maximum site induced-component change≤1e−6 Debye; retain exact-equality flags. For bare La, alpha=1.6597 Å³, position at the origin and intrinsic moment zero give native expected induced dipole `Debye*alpha*F/elefield` and energy `−electric*alpha*(F/elefield)^2/2`. Require each signed energy error≤1e−7 kcal/mol and each induced-component error≤1e−6 Debye. These nonzero expected responses detect a missing/no-op field parser, wrong sign, conversion or half-energy factor. Zero residual is valid. The ordinary parser's zero-field ion proof remains unchanged; only explicitly declared field-ion calls replace it with this analytic test.

Retain energies, full dipole vectors, all site moments and print bounds. Compute shared chi_E, chi_mu and odd mu_E at both steps for all four states, each trimer−dimer increment and Dy−La increment difference. Apply exactly the shared 0.5/1 a0³ consistency and 1e−4 ea0 per-state odd-slope gates. There is no physical response/sign target and no alpha3 or full-shell B claim.

## Runtime/collector contract

`RUN.py prepare --out <fresh preparation>` writes MANIFEST.json and CONFIG_RUN_TEMPLATE.json; root supplies a one-way release pin in the config. Release binds source, prepared manifest, shared SPEC, nineteen calls, one worker and fresh execution path. It does not hash a config that hashes itself.

`RUN.py run --config <bound config> --out <fresh execution>` writes full19 RESULT rows, per-cell START/RESULT/raw analyze/induced files, and a durable call counter. `admit_result(config_dict, manifest_dict, result_dict)` independently reparses saved outputs with the same unchanged native readers and returns rows plus summary without native calls.
