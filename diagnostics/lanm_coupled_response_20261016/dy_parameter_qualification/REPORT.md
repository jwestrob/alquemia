# Corrected La/Dy native parameters pass the real-core comparison

Slurm12488 completed all eight declared48-atom7CCO controls in19s on4CPUs.
Actual application time18.138s;76allocatedCPU-s. All eight endpoint admissions,
four absolute-energy/atomic-population comparisons and two solvent-component
comparisons pass their frozen criteria. See RESULT.json for every residual.

New full profile SHA256:
`635e111bac31e7c5755bf612c8e148ce73b401c110ff31887fd5e31b590fa85f`.
Only Dy kcn/shpoly/lgam change from the previously corrected-La full profile;
the policy lists both La and Dy corrections relative to the original export.
The source shell labels, every other field and all historical files are intact.

The exact finite plan is ../DY_QUALIFICATION_PLAN.md. This is originalGFN2,
effective multiplicity1/uhf0, charge−1, vacuum or ALPBwater, on one real peptide
core. Physical Dy(III) sextet metadata remains distinct. Nothing here qualifies
spin energetics, biological affinity or every element/charge/structure.

Native endpoint acceptance uses the shared parameter preflight/postflight,
actual imported-file evidence and all printed SCF criteria. Standalone checks
its actual method/version/state/settings, normal convergence, full populations
and analytic energy consistency. Neither normal termination alone nor array
matching alone was sufficient for admission. No arbitrary offset, numerical
tolerance change, retry or fabricated scientific result was used.

This eight-cell qualification uses the corrected profile only; it is not a new
bidirectional erroneous-versus-corrected causal test. Source schema evidence
identifies the misplaced arrays, and the real computation demonstrates recovery
of published-parameter standalone behavior under the tested conditions.

The full corrected file and policy are preserved here for new, explicitly named
paths. No production runner, La-only approval registry, PQQ calibration or cached
historical result has been replaced. SOURCES.json links exact raw inputs,
runtime/implementation pins, task-level observations and the original result.
