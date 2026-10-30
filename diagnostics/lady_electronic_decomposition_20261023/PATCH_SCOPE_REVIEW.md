# Independent review: CP2K Dy core55 reader support

Source-only review, 23 October 2026. No third-party source edited, installed, compiled or executed. Reviewed CP2K v2026.2 commit `67b5da876dd6a76b8b021d5a04d1c81ba79a4c50`, the archived Dy ECP, and the quartet input audit.

## Finding

**A correctly scoped core55 default-occupation case can preserve the intended semilocal ECP Hamiltonian.** The present failure occurs after coefficients and `zion = Z − Ncore` have been read, in the default occupation switch. A new case need only assign the missing occupation bookkeeping. It must leave charge, coefficient arrays, angular channels, integral routines, and all existing cases unchanged. This is a credible compatibility repair, not a validated executable or a new electronic result. [Reader source, lines2894–3068](https://github.com/cp2k/cp2k/blob/67b5da876dd6a76b8b021d5a04d1c81ba79a4c50/src/atom_types.F#L2894).

For the archived Dy(III) ECP55, the core is `[Kr] 4d10 4f9`: angular totals `(s,p,d,f)=(8,18,20,9)`. The stable periodic table gives neutral Dy `(12,24,20,10)`, leaving the default explicit occupation `(4,6,0,1)`, total11. Subtracting Xe plus one electron would describe the wrong core shells. The molecular/fragment charge then gives eight electrons for M3+; it must not be implemented by changing the ECP's effective nuclear charge11 or replacing the neutral default with an eight-electron core definition. [Periodic table, Dy](https://github.com/cp2k/cp2k/blob/67b5da876dd6a76b8b021d5a04d1c81ba79a4c50/src/common/periodic_table.F#L864).

## Downstream consistency and limits

`ELEC_CONF` is parsed earlier, but the ECP reader must return before that override can replace its default. After return, `qs_kind_types.F` passes `zion` and the exact local/semilocal coefficient arrays to the potential, then adopts the default only when no explicit occupation exists. Its existing `set_pseudo_state` already handles lanthanide Z58–71/core47–60 as `[Kr]+4d10+4f(Ncore−46)`. Thus core55/4f9 requires no second core-shell algorithm change. [Override and potential setup](https://github.com/cp2k/cp2k/blob/67b5da876dd6a76b8b021d5a04d1c81ba79a4c50/src/qs_kind_types.F#L2504), [existing lanthanide handling](https://github.com/cp2k/cp2k/blob/67b5da876dd6a76b8b021d5a04d1c81ba79a4c50/src/qs_kind_types.F#L3880).

Occupation bookkeeping feeds atomic initialization and shell setup; it is not cosmetic text. A different initial density can change convergence or the attained SCF solution even with an unchanged formal Hamiltonian. The default `(4,6,0,1)` is neutral pseudoatom bookkeeping, not a claim that the molecular Dy3+ state contains a valence 4f electron. Retain the declared effective singlets and molecular electron counts separately. General core55 support for arbitrary elements is not established; scope the proposed compatibility case to this Dy/core55 convention, retaining rejection of unsupported combinations.

## Minimal validation within the existing quartet

1. Inspect the exact inverse diff: only the missing default-occupation case, preferably guarded to Dy, with the existing nonnegative-occupation assertion retained. Show core count55, effective charge11, default sum11, and unchanged La46 branch. Do not alter LIBGRPP, coefficient/radial-power parsing, core compensation, ALMO equations or thresholds.
2. In the already declared quartet input preparation, verify every original exponent/coefficient/angular channel and Gaussian contraction through the format translation. Keep `ECP_SEMI_LOCAL` true. The no-change source argument does not validate a mistaken ORCA-to-CP2K conversion.
3. In the declared quartet calls, retain actual electron/spin/fragment counts: aquo8 88 electrons at +3; mixed acetate 110 at +2; M3+ fragment8, shells80/102. Report actual convergence/residuals and the full, polarized and CT-enabled states. Apply the predeclared original-energy agreement and decomposition-closure gates, with D4 kept separately consistent. No new atoms, geometries, parameter tuning or convergence sweep is needed for this patch review.

The two La members exercise an unchanged reader branch; the two Dy members exercise the added case. Agreement with archived ORCA full energies remains essential because code-port differences in XC/quadrature, exchange approximation, basis conventions and SCF state are not settled by this repair. A quartet pass would qualify this bounded use, not general CP2K core55 support or physical Dy spin-state accuracy.
