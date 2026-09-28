# Matching classical origins are ready

The reviewed finite classical ledger transfers to both remaining repaired sources without missing ion parameters or source remapping. Four source/target origins completed. The classical Dy-minus-La contributions are **+0.073347255864 kcal/mol for Hans8FNR** and **+0.066701362053 for Mex8FNS**. Their Hans-minus-Mex difference is only **+0.006645893811 kcal/mol**. This is a mechanical component, not a preference prediction or a solvent-complete binding model.

At fixed coordinates the bonded, MM-LJ and MM-Coulomb terms cancel exactly between targets. Only the target metal's cross-region LJ changes. Do not attribute specificity to the magnitude of those cancelling terms. Root can append these values to compatible electronic results, keeping electronic-only contrasts visible.

## Actual components (kcal/mol)

| Source | Target | Retained bonded | MM LJ | MM Coulomb | QM–MM LJ |
|---|---|---:|---:|---:|---:|
| Hans8FNR | La | 1530.276715820 | 109.668774820 | -4703.874110063 | -63.190777857 |
| Hans8FNR | Dy | 1530.276715820 | 109.668774820 | -4703.874110063 | -63.117430601 |
| Mex8FNS | La | 1405.369441224 | 153.217153695 | -5030.427048178 | -43.902913253 |
| Mex8FNS | Dy | 1405.369441224 | 153.217153695 | -5030.427048178 | -43.836211891 |

## Exact expression and exclusions

`E_classical = E_native,retained_bonded + E_MM,LJ + E_MM,Coulomb(q_shifted) + E_QM–MM,LJ`.
Native topology and parameters come from each source's own ff19SB/TIP3P XML, not Hans8DQ2 indices. Every real QM charge is excluded from classical Coulomb; omitted boundary CAs are physical atoms with zero MM charge. Native LJ exceptions and their 1–4 scales are retained; shifted-charge MM Coulomb exceptions are regenerated with verified native graph scaling. No QM–MM classical Coulomb, C4 induction, synthetic cap LJ or metal bonded constants are added.

The fixed-cap graph policy retains cut radial bonds and exterior terms, while omitting cross angular terms reproducible by the actual caps and their retained anchors. All four cut bonds are present in both sources. Source-specific decisions and counts are recorded in each `cross_term_decisions.json` and `LEDGER.json`.

The actual spectators are three Dy3+ in Hans8FNR and three Nd3+ in Mex8FNS. Their field identities remain unchanged. This classical component adds documented pure12–6 TIP3P IOD LJ to those same spectators. Exact parameter lines/file hashes are in the ledgers: La 1.718/.15060822, Dy 1.609/.08389240, Nd 1.681/.12564307 (Rmin/2 Å, epsilon kcal/mol), Lorentz–Berthelot mixing. These are available parameter candidates, not a validated cross fit to the frozen-f electronic Hamiltonian.

## Checks and saved gradients

Both real-source artifact suites pass: 2206 physical atoms for Hans and2049 for Mex; native ordering/charges, every source field row, selected QM membership, source A coordinates and caps, actual spectator identities, field closure (+9/+2), native exceptions and four cut bonds. Each configuration retains all four full Cartesian component gradients plus their sum. Energies/arrays were serialized immediately after evaluation, before aggregate checks. The first three components' energies and gradients are exactly identical across target metals.

An initial reporting test incorrectly demanded bitwise equality when re-summing gradients after JSON alphabetically reordered dictionary keys. The roundoff-only summation check now uses1e−10 kcal/mol/Å. No molecular calculation or scientific tolerance was changed or rerun. This is not new finite-difference qualification; the earlier Hans quarter-step result and original failed coarser checks remain separate.

Force-assembly interface: `source_ids`, `coordinates_A`, `gradient_kcal_mol_A` (gradient, not force; kcal/mol/angstrom), `component_gradients_kcal_mol_A`, `components_kcal_mol`, `ledger`, `inputs`, and explicit false QM_MM_Coulomb/C4 flags. Actual source EF3 IDs A/203//DY and A/203//ND label the physical target even when the evaluated target element is La. No synthetic caps enter these arrays.

## Cost, limitations and handoff

Slurm1220338 COMPLETED:5 elapsed seconds ×1 allocatedCPU = **5CPU-s**, mem=0, zeroGPU/quantum/optimization. Measured script wall4.679391s and process peakRSS190280KiB; scheduler MaxRSS reported0 and is not used as measured memory. Exactly16 component energy/force queries.

These are monomer, actual-occupancy, finite-crystal-water components. Source proton inventories/compositions and spectator species differ; no absolute cross-source minima are compared. Bulk solvation, validated cross parameterization and relaxed conformational populations remain absent. Full hybrid qualification remains false. Root owns combining matched electronic/physical gradients, interpretation and any further work.

Results: `workspaces/metal_environment_response_20260926/scaffold_origin_transfer_v1/RESULT.json`; per-source `{La,Dy}_A.json` contains full results and pins. No production changes.

Recheck without molecular work:
```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/metal_environment_response_20260926/scaffold_origin_transfer_20260928/TESTS.py
```
The accompanying batch command is an executed record; do not rerun the completed immutable output directory.
