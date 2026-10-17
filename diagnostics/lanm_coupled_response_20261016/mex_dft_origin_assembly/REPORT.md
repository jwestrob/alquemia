# Mex saved DFT origin assembly — complete

Job **12849** completed both La and Dy states with **zero molecular calls, new
geometries or optimization steps**. The saved native PBE0 electronic terms now
have an explicit whole-scaffold energy/gradient ledger:

`embedded PBE0(Q209; 1840 point charges) + existing MM + fixed-FF-charge ALPB`.

Both records retain all 2,049 physical atoms, four cap derivatives, waters and
three fixed Nd spectators. The electronic regions have charge 0 and effective
singlets, with physical La 1/Dy 6 metadata and ECP 46/55. Global physical charge is +2.
MACE and the GFN embedding correction are absent from this new sum; no additional
classical QM–MM Coulomb or C4 term is included.

## What passed, and what remains unresolved

| Saved replay | La | Dy |
|---|---:|---:|
| DFT output versus gradient-file energy residual, Hartree | 0 | 0 |
| MM component energy sum residual, kcal/mol | 0 | 0 |
| ALPB component energy sum residual, Hartree | 0 | 0 |
| Total gradient reconstruction residual, kcal/mol/Å | 2.84217e−14 | 1.42109e−14 |
| All printed SCF residuals within their printed tolerances | false | false |

The native calculations reported convergence under `Total+1el-Energy`,
`ConvForced 0`. Density, DIIS and orbital-rotation residuals exceed their printed
tolerances in both states; small-gap and symmetry warnings remain in the saved
evidence. The all-residuals flag is an observation, not a replacement admission
rule. No SCF rerun or changed setting was used.

Observed total translation-gradient sums, in kcal/mol/Å, are La
`[−2.03368e−6, −1.47398e−6, −3.98792e−6]` and Dy
`[−2.05740e−6, −1.54987e−6, −4.10413e−6]`. These are observations, not a new force
qualification. Reconstruction agreement verifies reuse and arithmetic; it does
not establish finite-difference consistency, electronic stability or physical
accuracy. Both `total_derivative_FD_qualified` and `physical_state_qualified`
remain false. Absolute energies for different metals are not affinity scores.

The copied `spectator_occupancy` is original preparation metadata and retains its
historical monopole-only wording. The **actual reused classical ledger** includes
Nd pure 12–6 parameters and the corresponding MM interactions; that ledger is
authoritative for this assembled model. No spectator term was removed or added.

## Actual cost and preserved evidence

Slurm: **COMPLETED, exit 0:0, 5 seconds, one CPU, 5 allocated CPU-seconds**.
Measured TotalCPU is **3.733 seconds**; batch MaxRSS is **450108K**. Application
elapsed/allocated CPU time is 3.914018 seconds. No GPU was allocated. Previously
completed molecular costs are not counted again in this saved-only job.

The exact result, both full assembled records, deployed source/config, execution
report, submission and completion snapshot are copied here. `PRESERVATION.json`
records their original paths and hashes; `SACCT.tsv` is the actual scheduler query.

| Artifact | SHA256 |
|---|---|
| `RESULT.json` | `b5322a5fd2666569d4c9dded8f99980ca977baf8c3c584e2870f64c93eb0618c` |
| Deployed source | `c19de37d038f6259c447049015247fc9fb51c2f2b76c0500c049a48d93f4a2fb` |
| Deployed config | `9c428c87a96323818ea9cad2a8262cdab4a7ce3c47edb8ffe0b8c48f2702fbac` |
| Declared plan | `9dfabb813214393c0748f9ce2ac99498960814cc9bb5f10b19ca9aac622254ab` |

This delivers traceable DFT/MM/solvent origins for a separately declared future
test within the active research goal. It establishes saved assembly correctness,
not a qualified relaxation surface, preference classifier or production default.
