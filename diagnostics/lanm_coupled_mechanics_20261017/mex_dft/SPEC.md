# Mex embedded-DFT total derivative qualification

Declared before new energies, under Module 4 of [the mechanics plan](../PLAN.md). Root owns all submissions and continuation decisions.

## Question and finite matrix

Does the analytic physical derivative of `PBE0_embedded(Q;MM) + MM_ledger(R) + ALPB(fixed_FF_charges;R)` agree with centered energy differences for an external-charge displacement and a cap-boundary displacement? PBE0 replaces both local MACE and the GFN2 embedding difference. No classical QM–MM Coulomb or full MM internal-Q energy is added.

Reuse the two admitted 2,049-atom origin assemblies from job12849 (`mex_dft_origin_assembly_v1/execution/RESULT.json`) and their original Q209 / 1,840-point-charge PBE0 gradients. The 2,045 native atoms, EF3 target and three Nd(III) spectators, original fixed charge shifts, exclusions, cap rules and unshifted whole-system FF/formal3 solvent charges stay unchanged. Whole charge +2, electronic region charge0. This is EF3 exchange at fixed Nd spectators, not pure-metal protein affinity.

Use the two source-defined directions already declared for job12709, without using forces to select them:

1. `MM_O97`: move A/97//O (physical1019, external row832) along the source C97→O97 unit vector. Electronic coordinates stay exactly unchanged; the point-charge position, MM geometry and solvent geometry change together.
2. `boundary_CA83`: move omitted cap anchor A/83//CA (physical801), along CA83→N83 projected perpendicular to C83→CA83. Reconstruct cap205, including both-anchor Jacobians. The external field stays byte-identical; MM and solvent use the same displaced physical coordinates.

For both La and Dy use identical physical displacements ±0.001 Å: eight displaced singlepoints total. No optimization, new coordinates selected from outcomes, or automatic retries. All eight slots remain in the denominator. No central-gradient calls are needed with the compatible saved origins; a specific compatibility blocker would have to be reported before any of the plan's at-most-two conditional central calls.

The first scout is the La MM_O97 ± pair, reused within the eight. Remaining six require actual pair admission, its derivative gate and root's measured tractability review. A failure stops dependent reference work; it does not authorize changed SCF settings or another step size.

## Electronic state and numerical recipe

Copy the actual admitted Mex input deck, removing only `EnGrad` for a singlepoint. Preserve PBE0-D4, def2-TZVP, def2/J, RIJCOSX, NoAutostart, DefGrid3, VeryTightSCF, PModel guess, `%method DoEQ false`, every custom basis/ECP/auxiliary byte and the finite external charge potential. La uses ECP46 and physical multiplicity1; Dy uses ECP55 and physical multiplicity6. Both execute the same spin-free frozen-f valence representation, effective multiplicity1, 832 explicit electrons (416α/416β). No SOC or physical open-shell-state qualification is inferred.

Use ORCA6.1.1 and the existing pinned native task runner, memory renderer and SCF evidence parser. One allocation runs one named endpoint with 12 actual MPI slots; MaxCore derives from that allocation's CPU/memory grant, with the established reserve. No Slurm memory/time overrides. Fresh directories and NoAutostart prevent density/restart transfer. An execution manifest contains only the endpoint actually run; the prepared inventory retains all eight.

Admission remains native normal termination, native SCF convergence and exact Hamiltonian/state/source/receipt checks. Record actual ConvCheckMode, ConvForced, all printed residual values/tolerances, warning/gap lines and separate all-residuals flag. The old origins converge under Total+1el-Energy / ConvForced0 while several density/DIIS criteria remain false. Those flags are not silently converted into a stricter or looser SCF gate.

## Frozen derivative gate and limits

Compute `(Eplus−Eminus)/(2h)` in each component's native units before conversion or summation. Compare PBE0, MM, fixed-FF ALPB and total against the archived physical gradients. Each residual must be at most `0.005 + 0.0001*max(abs(FD),abs(analytic))` kcal/mol/Å. Retain separate MM-group comparisons and both-metal/common/differential values; the gate applies to each own-metal component and total. The 0.005 absolute target is below the optimizer's 0.03 RMS and 0.10 maximum-force targets.

The fixed 0.001 Å step reduces energy-subtraction sensitivity relative to the earlier approximate-model 0.0001 Å fixtures. This eight-cell design has one step size: truncation error and step-halving stability remain unqualified. Report printed-energy quantization and binary64 subtraction scales, and a clearly labeled sensitivity proxy derived from printed SCF energy tolerances/residuals. SCF stopping thresholds are not rigorous energy-error bounds. A derivative pass cannot overrule unresolved SCF/density/state uncertainty or establish physical accuracy, full-gradient accuracy outside these two directions, or stable electronic states.

## Implementation and actual costs

`QUALIFY.py prepare` replays saved origins and writes eight current physical/core/field maps and decks, with zero energies. `singlepoint --task ID` reuses the existing one-cell MM primitive and native fixed-charge solvent executable, then runs one native DFT endpoint. `collect` admits saved endpoint receipts and reports either the scout pair or all eight; it launches no molecules. Configs pin actual preparation/results, every helper, source and runtime asset. No new engine or production default.

Record attempted/complete DFT, MM configurations/group queries and ALPB calls; process/allocation costs remain distinct from Slurm accounting. The old Mex La/Dy origin outputs report 1h28m53.155s / 1h39m57.052s molecular wall time under historical allocations. These are context, not predictions for 12-slot displaced singlepoints. Root uses the measured first pair before the remaining six.
