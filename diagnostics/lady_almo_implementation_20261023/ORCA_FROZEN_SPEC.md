# Frozen ORCA density: RIJCOSX diagnosis

## Question and finite scope

The converged La-aquo PySCF baseline disagrees with the archived ORCA reference.
Measure the effect of native RIJCOSX versus direct Coulomb/exchange on the exact
archived ORCA orbital set. Exactly two independent native starts are declared:
one RIJCOSX control and one NOCOSX/NORI counterpart. Each computes one Fock
matrix with `MOREAD CALCGUESSENERGY NOITER`, zero SCF iterations and no gradient.
These are fixed-density diagnostics, not new converged-SCF reference energies.

`orca_frozen_diagnosis.py prepare` pins the original La-aquo input/output/GBW,
all geometry copies, complete basis/ECP and auxiliary assets, original result,
ledger, native executable and relevant executable modules/shared libraries.
It copies the GBW and geometry without changing the originals and preserves the
complete literal `%basis` block. All pins are checked before and after a run.

Preserved settings: PBE0, D4, def2-TZVP light basis, supplied La orbital basis
and ECP46, auxiliary basis, exact coordinates/order, charge +3, effective
singlet, DefGrid3, VeryTightSCF and MaxCore3000. The original four MPI ranks
become eight. The PModel/MaxIter300 block becomes MORead/MaxIter0. EnGrad becomes
SP. The direct member replaces RIJCOSX with NOCOSX NORI and `%method RI off`.
No functional constants, integral thresholds, grid settings or scientific
admission thresholds change.

## Resources and admission

Root owns submissions. Each native job uses standard_expand, eight allocated
CPUs, one numerical thread per MPI rank, 32 GiB maximum allocation and one hour
maximum job wall time. Native timeout is 900 seconds, with termination of its
process group and no retry. The archived complete four-rank energy/gradient
calculation took 374.269 seconds; this informs the finite cap without promising
direct-J/K cost. A manifest-level exclusive attempt receipt prevents repeating
either variant. Fresh output directories and basenames protect saved orbitals.

Admission requires normal native termination, exact source/input pins, one
affirmative one-Fock message, explicit MaxIter0/SCF-skipped messages, identical
native geometry/basis, imported occupations, no removed overlap modes, and
unchanged 88-electron/391-AO/ECP46 inventory. Unexpected SCF iteration tables,
convergence banners, gradients or basis projection fail closed. Same-basis
native renormalization/reorthogonalization is retained and reported; bitwise
density equality is not claimed. Original and emitted GBW files are preserved.

## COSX phase and interpretation

DefGrid3 is unchanged. Capture every printed DFT/GRIDX section and any final
reintegration/phase-switch messages. A generated grid list alone does not prove
which COSX grid the one-Fock operation used. Report actual control-minus-archive
electronic and D4-inclusive residuals, together with the pre-existing 1e-4 Eh
archive gate. Passing that gate means agreement within that tolerance, not exact
identity or proof of final-COSX-grid reproduction. A failed control comparison
remains visible and prevents interpreting its paired difference as the whole
archived discrepancy. No phase override, grid sweep or automatic follow-up.

The paired difference measures the native energy-operator change at imported
orbitals under NOITER conventions. A converged direct SCF would also include
orbital relaxation and is outside this scope. This diagnostic cannot admit the
failed baseline, change ALMO gates, establish physical state/accuracy, or prove
that all remaining ORCA/PySCF differences have disappeared.

## Evidence

- [ORCA 6.1 RI controls](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/RI.html#keywords): NOCOSX retains RIJONX unless RI is also disabled.
- [ORCA input-density energy](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/initialguess.html#calculating-only-the-energy-of-an-input-density): CalcGuessEnergy/NoIter evaluates one Fock without orbital updates.
- Existing same-version native precedent: `workspaces/mace_omol_20260917/chelpg_sampling_qualification_v1/GGR_2FVY_La_original/endpoint.out`, lines914–944, supplies the affirmative log phrases used for admission.

No native calculation is performed by preparation or validation.

## Observed MPI startup failure and scheduling-only recovery

22432/22433 failed in orca_startup_mpi before a Fock evaluation: one Slurm
task with eight threads provides only one MPI slot. Preserve both failed
attempts. Retry the identical two native inputs with eight tasks, one CPU each,
32 GiB, one hour and the same native 900-second cap; standard_expand. A fresh
manifest and recovery record bind the originals. No oversubscription, changed
method, new geometry, extra successful energy evaluation or gate change.
