# Reciprocal metal-position diagnostic: independent review

## Evidence motivating the test

The completed eight-cell ExaA vacuum reference reproduces the source dependence of the saved donor response. The actual Slurm comparison reports:

| Ca-minus-La donor work, kcal/mol | MACE vacuum | r2SCAN-3c vacuum |
|---|---:|---:|
| Control sample 0 | +5.261497 | +5.377921 |
| Target sample 3 | +16.403481 | +17.135466 |
| Target minus control | +11.141984 | +11.757545 |

The stored discrepancy in the last row is +0.615561 kcal/mol. These are common origin→adaptive-La paths, not operational candidate minima or binding affinities. Agreement weakens a MACE-specific explanation for the large target response and makes a bounded coordinate intervention useful.

The reference remains conditional. Five endpoints finish with admitted TRAH residuals. Control Ca origin instead explicitly stops on a gradient check at 1.4799e−6; its final energy-change summary exceeds TolE. Control La origin and target Ca proposal stop on energy checks with final gradients above 1e−5. Consequently, all eight satisfy the frozen solver-specific diagnostic rule, but not a universal all-criteria convergence standard. Late gap/level-shift warnings, untested state stability and the dry oxidized-PQQ model remain limitations. No new numerical failure was found in the completed control proposals.

Reference record: `workspaces/exaa_donor_response_20261014/full_comparison_v1/execution/RESULT.json`, SHA256 `3a1502e41a6c5a349f5d6aa5e9f817a1bab7ecdd355e9f1ecd5cfd32f0b6eb80`. It pins all eight actual collections; this review only read their existing results and raw endings.

## Causal question and limits

Align the other source's 24 named PQQ heavy atoms to a host's original PQQ by a proper rotation, then transfer only that source's metal point. Apply the resulting single host-frame metal point to both host geometries and both metal identities. This asks whether PQQ-relative metal placement changes the metal dependence of the **same host donor motion**. It does not transfer the donor coordinates or reoptimize them.

If the control-derived point reduces the target response and the reciprocal target-derived point increases the control response, placement is implicated within this fixed local model. Either direction, an asymmetric response or no response is informative; no desired sign is an admission criterion. Changes can reflect metal interactions with every retained host group, not only PQQ. A successful reciprocal effect still cannot establish which original placement is physically correct or prove that either fold is wrong.

Preserve each host's proton inventory, state, charge, caps, nonmetal coordinates and spectator treatment. No solvent computation is part of the proposed eight MACE endpoints; cavity equivalence or ALPB accuracy is therefore not established. The existing frozen solvent classifier remains unchanged.

No molecular preparation, model call, scientific algebra or submission has been run by this reviewer.

## Scoring adapter review

Static review of `scoring/SCORE.py` SHA256 `7e815ca00bae318116b58832a0c63c4fe8e546439441d69f24e2acaa8667c9ad` found no execution or energy-sign blocker. It reuses the unchanged, byte-pinned `CPU_BATCH.WarmCPU` and `cell` interfaces with one seven-thread model worker. Actual packing qualification schemas match the adapter's assumptions. Source, state, checkpoint and runtime identities are verified; the eight archived baseline energies are admitted from their original request/result records rather than reevaluated.

The new XYZ checks enforce unchanged nonmetal coordinates relative to each host candidate, paired Ca/La coordinates and one identical transferred metal point across both host candidates. The collected individual works, Ca-minus-La differences, per-host `delta_D` and target-minus-control contrasts follow the declared signs. No descriptive ±1 kcal/mol hypothesis rule alters scientific admission. A failed model cell has a null energy and blocks the dependent contrast; no partial minimum, hidden retry or old-energy substitution is introduced.

The current metadata binder requires two prepared rows before scoring. Thus a preparation failure must remain a visible unavailable two-source experiment in the preparation/root report; the binder must not be used to quietly select one successful host. Scoring launcher SHA256: `d6a054a210140fcd0645517fd7c13d3671a048fe6da918011fbce36de5c5c14e`.

## Frozen preparation review

**No concrete blocker to the one-CPU Slurm preparation was found.** Actual metadata schemas match the adapter's accesses: both sources have 192 atoms/14 caps, the original admitted four-angle `full_q`, and three single-metal Cartesian translation modes with the expected axes. The output row/cell schema matches the scoring adapter.

The named 24-atom PQQ correspondence checks uniqueness, elements and whole-cofactor membership. One row-vector proper Kabsch rotation (`U D Vt`) and its analytic inverse generate the reciprocal points. The second singular-value condition rejects a collinear alignment while allowing planar PQQ. The fit uses original PQQ coordinates, and their equality along the host donor motion is checked; the code does not refit around each candidate or metal.

Original host `Kinematics` and the saved donor coordinates are retained. The existing frozen `adaptive_metal_proposals.final_geometry` is called for both old and transplanted arrangements. It enforces original-origin heavy displacement and coordinate bounds, cap/source/covalent consistency, and the original `geometry_check` rejection of **new <1.0 Å non-H pairs or <0.55 Å H-involving pairs**. This is not the packing branch's weaker 0.45 Å duplicate-coordinate check. The frozen 0.8 Å physical bound and its existing numerical tolerance are retained; none is adjusted to admit a transplant.

Only the metal coordinate line and explanatory XYZ comment are rewritten. Nonmetal coordinate lines remain byte-identical, physical reconstruction must agree within the frozen replay tolerance, and both metals share the resulting point. Failure within a host invalidates all four host cells, retaining partial artifacts and the full two-host/eight-cell denominator. No geometry failure is converted into a valid unchanged-origin substitute.

This is **static implementation clearance**, not a passed physical-preparation result. Numerical construction and gates have not been executed by this reviewer. Root must inspect the actual Slurm preparation before binding scoring.

Pins:

- Manifest: `workspaces/exaa_metal_position_20261014/preparation_v1/manifest.json`, SHA256 `28aebabb3f4a99f00b80b708d446e3ace61f0adf4c836415beec3a4d7957ce39`.
- `prepare/PREPARE.py`: `6e33e953be4f2d4fcad57797e491a9c82bb5c77b8435861e6f7642f3d588d87d`.
- `prepare/RUN.sbatch`: `ce3fd1064d37dea673bb148dee250e6e5bfcb0d42b087d2600fd54d0d7edd6d9`.
