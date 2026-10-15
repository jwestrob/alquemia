# ExaA PQQ residual-force feasibility

**Saved MACE projection is possible; matched saved DFT-force comparison is unavailable.** All eight r2SCAN-3c calculations were single-point energies. Their inputs request `r2SCAN-3c NoAutostart DefGrid3 TightSCF`, with no nuclear-gradient keyword; no `.engrad`/`.grad` file or Cartesian-gradient section exists in any endpoint directory/output. Each retained a GBW. Orbital-gradient convergence diagnostics are not nuclear forces. The completed finite donor-work comparison therefore supports that finite energy path, not uncomputed PQQ forces.

[full inventory pointer](INVENTORY_POINTER.json) pins all eight exact endpoints from [the completed DFT result](../../exaa_donor_response_20261014/FINAL_RESULT.json), the corresponding MACE receipts/requests/force arrays, and source maps. All eight MACE force-file hashes and coordinate/state joins passed metadata verification. No force array or coordinate array was numerically loaded. DFT electronic stability and gradient qualification remain unestablished.

## Exact cofactor identity

Both sources have the same27-atom intact PQQ block, zero-based indices165–191:24 heavy atoms plus `H_N1`, `H_C3`, `H_C8`. State is the existing modeled oxidized PQQ3−; this task changes no chemical state. Every current donor-mode moving-index list excludes this block, so these degrees of freedom are absent from the existing search.

The generic names follow the pinned `protenix_generic_pqq_v1` graph. The existing cofactor audit maps:

| Generic atoms | Indices | CCD atoms | Role |
|---|---|---|---|
| C5–O1 |169–170| C5–O5 | Metal-facing quinone C–O |
| C6–O2 |171–172| C4–O4 | Other quinone C–O |

`H_N1` belongs to N1; `H_C3` to generic C10; `H_C8` to generic C1. Source naming must not be replaced by CCD labels without the explicit map. The heavy covalent graph/name map is pinned to `workspaces/pede_solvent_response_20261013/cofactor_review_v2/manifest.json`; the actual source/core preparation records supply the three H parents. The context treats PQQ as an opaque fixed cofactor, so its protein-only context bond list is not its complete covalent graph.

## Frozen finite saved-force analysis

Root has frozen the following definitions and numerical tolerances in [PLAN.md](PLAN.md) before any projection. Analyze precisely the eight saved cells: samples0/3 × origin/adaptive_La × Ca/La. One shared Slurm CPU is sufficient; zero model/energy/gradient calls are required.

1. For each quinone C–O, define u from C toward O. Hold C and all other atoms fixed; displace only terminal O along u. The generalized force is `F_O·u`, positive toward elongation; its negative is dE/dlength. Report actual bond length only inside the declared Slurm analysis.
2. Treat all27 PQQ atoms, including their three existing hydrogens, as one rigid object. Report total force and torque about the centroid of its24 heavy atoms. Use a chemically fixed right-handed frame: e1=C5→O1; e2 is C5→C6 orthogonalized against e1; e3=e1×e2. Normalize each direction by its heavy-atom displacement norm, retaining raw force/torque and raw norms separately. Include H in the rigid motions and force sums. No fitted axis or favorable direction is selected from observed forces.
3. Preserve the full per-cell values. Compare Ca−La at the same geometry, their arithmetic common load, proposal−origin for each metal, and the resulting source contrasts. Check exact paired geometry, fixed PQQ coordinates across donor motion, complete source/name/element correspondence, finite192×3 force arrays, and the stored negative-energy-gradient convention. Require frame/norm residuals below1e−10 and projection invariance under a fixed joint coordinate/force transform below1e−9eV/Å.

This would diagnose residual MACE force outside the available donor modes. A force magnitude or directional change alone cannot establish an erroneous fold, a unique strain mechanism, beneficial relaxation, a redox change, or a classification remedy. It would not create missing DFT forces or qualify the solvent gradient. A real DFT force comparison would require a separately declared analytic-gradient evaluation; none is prepared or authorized by this inventory.

## Audit scope

Only authoritative JSON, input/output text, file existence and hashes were inspected. No molecular calculations, numerical force projections, imports of scientific backends, submissions, structure edits or broad archive scan occurred. Inventory SHA256: `6bacbe1c78d38061cd11801e992f745ece05bb151f68e2a6bf0fcd769523d1af`.
