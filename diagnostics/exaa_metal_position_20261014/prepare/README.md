# Reciprocal metal-position preparation

The runnable snapshot is frozen in `workspaces/exaa_metal_position_20261014/preparation_v1/implementation/`. Its manifest SHA256 is `28aebabb3f4a99f00b80b708d446e3ace61f0adf4c836415beec3a4d7957ce39`. Preparation has not run when this record is written. Metadata generation and Python syntax parsing made no scientific imports or calculations.

Root submission command:

```bash
sbatch --nodelist=node-64-768g-5 \
 /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/exaa_metal_position_20261014/preparation_v1/implementation/RUN.sbatch \
 /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/exaa_metal_position_20261014/preparation_v1/manifest.json \
 28aebabb3f4a99f00b80b708d446e3ace61f0adf4c836415beec3a4d7957ce39
```

The wrapper requests one standard-partition CPU with no memory or time override. Runtime checks require an actual one-node, one-CPU allocation and one numerical thread before importing any scientific code. Root owns submission and completion wakes.

## Exact construction and admission

- Original192-atom two-Ca samples0/3, exact inventory SHA `757a13523549697fd52cf9e94b3ba2ffdd01c463f3969fd8a5d01c1459d28343`.
- One unweighted proper-rotation fit of all named PQQ heavy atoms, `C1–C14`, `N1–N2`, `O1–O8`, from source0 to source3; the inverse supplies the reverse transfer. Planar rank2 geometry is allowed; collinear/undefined alignment fails. Fit residuals are reported, not tuned to score outcomes.
- The resulting point for each host is identical across origin/adaptive_La and Ca/La. The writer changes only the metal coordinate line and diagnostic comment; every nonmetal coordinate line is retained exactly.
- Exact frozen `adaptive_metal_proposals.final_geometry` is reused through the original source snapshot. Three metal axes are added to the original four selected donor coordinates. The saved donor q is unchanged. **All displacement and overlap checks remain relative to the original host origin**, including0.8Å physical bound with its existing tolerance, component bounds, source/cap/bond/Jacobian checks and new-overlap guards1.0Å heavy/0.55Å H-involving. No0.45Å packing-only substitute or rebasing to the moved donor endpoint.
- Same source maps, state signature, protons, dry inventory and original electronic charges/spins. The small rendered-versus-saved tolerance verifies kinematic representation; it never authorizes altering saved nonmetal coordinates.

## Downstream contract

Output `execution/PREPARATION.json` has two rows, each with `source_case_id`, `sample`, `role`, unchanged maps/state/context pins, a `geometry_admission` receipt and `matrix[Ca|La][origin|adaptive_La]`. Prepared cells expose XYZ/state/atom count and their original XYZ/MACE/solvent receipts. `GEOMETRY_ADMISSION.json` has `pass:true` only after all four source cells pass, including paired coordinates. A failed host retains four unavailable cells; technical global failure retains the complete two-source/eight-cell denominator and an execution traceback. There are no energy calls, optimization, source repairs or silent retries.

The result is a modeled coordinate perturbation, not an observed structure or equilibrium sample. Root must inspect actual preparation before binding the separate eight-call MACE batch. Independent static review found no blocker; execution remains unqualified until the Slurm result exists.
