# Saved-data checks and numerical comparison

Use `collect.sh` in a one-CPU Slurm allocation after the relevant jobs finish.
It launches no model or quantum calculations. Preserve raw arrays and verify
their hashes before analysis. Assert exact endpoint coordinate correspondence,
the two intended element changes, float64 forces, finite arrays, charge sum and
the recorded charge/spin/head checks. Reproduce the distant-water extrema and
the existing directional finite-difference check on real preparations.

Report the declared periodic translation and each remote-water deletion
separately, including metal force-vector changes. Retain the original 0.5
kcal/mol robustness criterion. Removing a water changes real electrostatics;
an exceedance requires diagnosis before being assigned to a model defect.

Before viewing GPU values, require CPU/GPU agreement on both completed start-0
endpoints within 1e-5 eV in total energy and 1e-5 eV/A in each force component.
These are numerical implementation checks, not chemical accuracy tolerances.
Retain every initial fixed-configuration energy gap but explicitly mark it as
not a free energy. No comparison to the conditional experimental contrast is
performed on those gaps.

Bibliography checks follow CODE1's one-based first-record pointer to the native
ligand identifier and use CODE1's separate logical identifier for LMBIB. Verify
the four predeclared ligands only. This recovers provenance without broadening
the target search or redistributing the native compilation.
