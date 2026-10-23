# Matched chemistry comparison

Question: on the exact shared geometries, how do current MACE and the two
declared PBE0 recipes differ in balanced water-to-donor exchange contrasts,
radial deformation work, and the change in forces upon La-to-Dy substitution?

Inputs are the pinned reference panel, the complete native MACE result, and
available admitted native reference results. Missing reference cells remain
missing. No molecular calls, new fit, geometry change or free-energy claim.
Run saved-data arithmetic on one local CPU within root's allowance. Each
collection uses a fresh output directory and records its inputs and elapsed time.

The optional completed native polarizable panel adds the same 28 physical
configurations, retaining its explicit unqualified Dy parameter-transfer label.
Its four free fragments are not included as extra reference configurations.
Both models receive identical balanced-contrast and paired-force comparisons.

The balanced fixed-configuration contrast is E(DyL)-E(DyW)-E(LaL)+E(LaW).
Its negative sign means greater Dy stabilization relative to water on these
configurations. Free-fragment and isolated-metal terms cancel. The 8-to-9-water
contrast is balanced in the same manner. Report all three radial configurations,
not only the most favorable. Full-force paired error compares (F_Dy-F_La)
between methods, alongside individual-metal errors. This avoids mistaking large
shared force errors for equally large metal-discrimination errors.

Raw DFT recipe admission is not electronic-state or accuracy qualification.
Do not compare uncorrected dimer totals to published counterpoise interaction
energies. No acceptance threshold is fitted to these results.
