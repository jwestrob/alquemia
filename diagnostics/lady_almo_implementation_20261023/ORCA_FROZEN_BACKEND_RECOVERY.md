# Saved-only correction of the auxiliary-metric parser veto

Native jobs22436/22437 completed the declared one-Fock evaluations normally.
The direct receipt failed only `requested_jk_backend`: the parser treated any
startup `Calculating RI/J V-Matrix` message as active RI in the Fock build.
Preserve that failed receipt and both immutable source/input bundles.

Independent evidence supports correcting this predicate for these exact saved
artifacts, without another native call:

- Direct output714 computes the auxiliary metric during SHARK startup, before
  SCF SETTINGS762 and one-Fock execution883. Its Hamiltonian omits RI/COSX
  activation, its exact input requests NOCOSX/NORI and `%method RI off`, and
  no COSX grid is generated.
- Native `frozen_direct.property.txt:456` records a positive `J` timing and no
  `SPLITRIJ` or `COSX`. The corresponding control property file records positive
  `SPLITRIJ`/`COSX` timings and no `J`. These independently identify the actual
  SCF work, beyond input keywords and startup preparation.
- Generated method bibliographies corroborate the distinction: RI Coulomb
  DOI10.1002/jcc.10318 and COSX DOIs10.1016/j.chemphys.2008.10.036 and
  10.1063/5.0058766 appear in the control and are absent from direct.
- The [ORCA RI manual](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/RI.html#ri-j)
  defines V as the auxiliary-basis metric, describes its startup Cholesky
  preparation, and documents `%method RI off` as disabling RI-J. Preparing
  this metric alone is not evidence that the density-dependent Fock uses it.

`recover_orca_frozen_backend.py` is a one-CPU, Slurm-only saved collector. It
pins the exact source, manifest and both result/output/error/property/bibliography
sets, verifies all original artifact pins, and calls the original parser.
The sole override is `requested_jk_backend` for job22437, conditional on the
independent native evidence above and every other unchanged admission check.
It rejects any additional failed check. It makes no native call and never
edits old artifacts. Common printed DFT grid and density observables must match.

Report paired electronic/full/D4 differences, the control-minus-archive
residual, actual grid sections and the unchanged 1e-4 Eh comparison gate.
Do not call the NOITER result converged SCF or assert exact final-COSX-phase
reproduction. The failed PySCF baseline and all ALMO gates remain unchanged.
