# Hans native mechanics: accounting passes; hydrogen strain is substantial

**The exact-source classical partition is numerically correct, but the normalized source is far from mechanical equilibrium, mostly because of hydrogen angles.** Address that preparation issue before interpreting a whole-scaffold relaxation or a harmonic compliance model. This result does not establish a La/Dy preference or invalidate the separate paired electronic response by itself.

Job1220314 evaluated only the existing1,887-atom ff19SB/TIP3P protein/water System, at the two already prepared A/B configurations. No ions, synthetic caps, electronic energies, new coordinates, optimization or parameters were added. The unmodified native parent and independently reconstructed QM/cross/MM groups agree to2.31e−10kcal/mol (fixed acceptance~1.53e−6). All native exceptions and CMAP tables were preserved.

## What moves, and what does not

A→B native parent work is−0.667753254kcal/mol:

| Component change | kcal/mol |
|---|---:|
| Wholly-QM nonbonded |−0.954550114|
| Cross QM/MM nonbonded |+0.289270060|
| Wholly-QM torsions |−0.002473200|
| MM-only terms |0|
| Cross bonded terms |0|
| Other QM bonded terms |roundoff|

This confirms that the chosen intact Asp85 motion changes local contacts without deforming covalent bonds/angles or the exterior. It does **not** test scaffold compliance. The native cross nonbonded term includes classical Coulomb as well as LJ; it must not be added wholesale to embedded QM, which already represents electrostatic coupling. Native parent charges also differ from the electronic boundary redistribution. No hybrid correction is provided.

## Source strain and immediate implication

Native A energy is5331.793845kcal/mol; its zero has no binding interpretation. Total harmonic-angle energy is4922.935452kcal/mol. Direct source/XML angle decomposition places4846.850908kcal/mol in H-containing angles versus76.084544kcal/mol in all-heavy angles. Largest actual examples:

| Source angle | Actual degrees | Native equilibrium degrees | Energy kcal/mol |
|---|---:|---:|---:|
| CA79–CB79–HB2 |172.539|109.5|60.526|
| N124–CA124–HA2 |168.768|109.5|53.501|
| HG3(110)–CG110–CD110 |167.659|109.5|51.518|
| CA83–CB83–HB2 |166.711|109.5|49.853|

The maximum parent force is215.159kcal/mol/Å at residue73 HB3. The largest individual bond strain is only0.471kcal/mol. Thus corrected radial H bond lengths do not establish plausible H directions. Exact topology/ID ordering was verified in the preceding feasibility audit; this is not evidence of a missing ion parameter, and the dominant angle contribution does not require ions.

The source is not stationary under its native parent model. A response approximation must retain affine forces; blindly interpreting its curvature as equilibrium compliance would be misleading. The next useful preparation investigation is a source-mapped hydrogen-direction audit/repair that preserves actual heavy atoms, atom identities, proton inventory and water occupancy, evaluated separately from old frozen inputs. Do not rewrite current electronic artifacts or use the classical residual to tune a desired metal sign. This test does not uniquely diagnose prior electronic or whole-protein failures.

## Execution and reproducibility

- OpenMM8.5.1 Reference platform, two configurations,32 energy-state queries and2 parent-force queries; no integration steps.
- Actual scheduler:4wall seconds,1CPU,4allocatedCPU-seconds. Script measured1.567s and95,700KiB process peak RSS. Slurm step MaxRSS reported0, so that value is not treated as measured zero memory.
- **Resource-policy deviation:** this tiny job requested2GiB, contrary to the standing user mem=0 preference. It completed before correction; no cancellation/repetition occurred. The request and receipt remain preserved. Any future submission must use mem=0; the archived batch is evidence, not a future submission template.
- Full result: `workspaces/metal_environment_response_20260926/scaffold_decomposition_v1/RESULT.json`; angle audit alongside `ANGLE_AUDIT.json`; scheduler evidence `ACCOUNTING.txt`.
- `RUN.py` and `PLAN.md` are pinned by the actual result. `ANGLE_AUDIT.py` is a direct algebraic inspection of the same source/parameters, no Context or additional configuration. It verifies the dominant angle term independently.

No rerun is required. Read the saved result/angle audit. A future altered scientific model or new coordinates require a distinct version; this native parent is not a complete hybrid Hamiltonian and its missing ion terms remain unavailable rather than zero.
