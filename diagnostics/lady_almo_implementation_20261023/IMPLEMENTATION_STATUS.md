# ALMO implementation status

Updated after terminal jobs22344/22439 on24October2026. The reusable solver is
implemented and passes its core numerical tests. The archived La/Dy application
remains unvalidated: ordinary-reference agreement is unresolved, and both
independent native benchmark attempts stopped before ALMO. No quartet
charge-transfer result is admitted. No jobs remain in this validation round.

## Implemented mathematics and interface

[almo.py](almo.py) minimizes the full PySCF restricted HF/DFT energy while each
fragment's occupied orbitals retain support on that fragment's AOs. With real
occupied coefficients `C`, AO overlap `S`, and `M = C.T S C`, it uses

```text
D = 2 C M^-1 C.T
dE/dC = 4 (F C M^-1 - S C M^-1 C.T F C M^-1).
```

The density includes the interfragment occupied-overlap contribution. Local AO
whitening, polar retraction, and projected L-BFGS preserve the coefficient
constraint. Admission requires electron count, density idempotency, support and
rank checks, plus both orbital-gradient and energy-change convergence. Failed
searches return the last accepted state with trial history; invalid initial
inputs fail explicitly. [THEORY.md](THEORY.md) explains the derivation and scope.

Public API:

```python
problem = ALMOProblem(mf, fragment_aos, fragment_nocc)
evaluation = problem.evaluate(C)  # energy, density, Fock, dE/dC, invariants
result = solve_almo(
    mf, fragment_aos, fragment_nocc,
    initial_coefficients=C,       # alternatively: fragment_mfs=[...]
    max_iterations=300, max_evaluations=1500,
    gradient_tolerance=1e-6, energy_tolerance=1e-10,
)
released = release_constraints(mf, result)  # ordinary SCF from result.density
```

`fragment_aos` partitions the complete AO list; `fragment_nocc` counts doubly
occupied spatial orbitals. `result` includes coefficients, density, energy,
convergence/status, gradient, histories, invariants and evaluation counts.
`release_constraints` explicitly updates the supplied mean-field object.

## Actual saved evidence

- **PBE0 core tests, job 22332:** all checks passed. Four orbital directions at
  two finite-difference steps agree with the analytic derivative; largest error
  is `1.94e-8` Hartree. Density reconstruction, occupied-block gauge invariance,
  invalid-rank/support rejection, variational ordering and constraint release
  pass. The constrained solve used 25 iterations/27 evaluations and reached
  gradient norm `4.996e-9`. One-fragment recovery differs from ordinary SCF by
  `1.71e-13` Hartree.
  [Saved result](../../workspaces/lady_almo_implementation_20261023/core_tests_v1/RESULT.json)
- **PySCF PBE benchmark side, job 22333:** passed, with gradient norm
  `4.844e-9` after 19 iterations/20 evaluations. Independent CP2K agreement
  remains unvalidated because neither native run reached its ALMO stages.
  [Saved result](../../workspaces/lady_almo_implementation_20261023/pyscf_benchmark_v1/RESULT.json)
- **All four unique isolated fragments pass native SCF admission:**
  [La, 22334](../../workspaces/lady_almo_implementation_20261023/fragment_La_v1/result/RESULT.json),
  [Dy, 22411](../../workspaces/lady_almo_implementation_20261023/fragment_Dy_1e_v2/result/RESULT.json),
  [aquo shell, 22336](../../workspaces/lady_almo_implementation_20261023/fragment_aquo_shell_v1/result/RESULT.json),
  [mixed shell, 22337](../../workspaces/lady_almo_implementation_20261023/fragment_mixed_shell_v1/result/RESULT.json).
  Earlier Dy initialization failures remain preserved.

Three ordinary molecular references converge and pass their density/gradient
checks, but each fails the unchanged `1e-4` Hartree archive-agreement gate:

| Reference | Job | PySCF minus archived electronic energy, Hartree |
|---|---:|---:|
| [La aquo](../../workspaces/lady_almo_implementation_20261023/baseline_aquo_La_v1/result/RESULT.json) | 22327 | +0.000825580291 |
| [Dy aquo](../../workspaces/lady_almo_implementation_20261023/baseline_aquo_Dy_1e_v2/result/RESULT.json) | 22409 | +0.000956700399 |
| [La mixed](../../workspaces/lady_almo_implementation_20261023/baseline_mixed_La_v1/result/RESULT.json) | 22329 | +0.000893668327 |

Mixed-Dy22410, native benchmarks22343/22344 and saved comparison22439 are
terminal. Their failure details follow below; no result is admitted merely
because a backend reported convergence or its energy appeared stationary.

## Remaining work and limits

1. Independent native validation remains unavailable after both bounded attempts
   stopped before ALMO. A longer unchanged run is not supported by the observed
   residual plateau. Preserve this distinction from failure of the ALMO solver.
2. Resolve the ordinary-reference mismatch through the separately recorded
   backend diagnosis; preserve failed comparisons and the original tolerance.
3. Once the required gates pass, run the declared constrained quartet using
   admitted reusable fragment/reference states, then collect its decomposition.

Convergence establishes constrained first-order stationarity, not a global
minimum or physical electronic-state accuracy. The effective singlet Dy model,
finite basis, and fragment definition remain scientific assumptions. The present
software evidence does not establish La/Dy charge transfer or exchange free energy.

## Subsequent backend diagnosis

Saved ORCA one-Fock comparison22441 passes after correcting a startup-metric
parser error. Direct J/K minus RIJCOSX is +0.461420371kcal/mol at the same
imported orbitals; the RIJCOSX control reproduces the archive within the
existing tolerance. This identifies a large numerical-backend contribution.
Documented PBE constants separately account for8.5034e-5Eh on the fixed
PySCF density. Neither observation alone admits the original failed baseline
or supplies a charge-transfer result. See [CROSS_CODE_AUDIT.md](CROSS_CODE_AUDIT.md).

## Native400Ry benchmark terminal update

Job22343 failed in the first isolated-water SCF at300iterations before ALMO.
Printed residual1e-8 exceeds the original1e-10 criterion although the final
energy change is4.97e-13Eh. Integrated net charge is-0.0306942524, a separate
numerical discrepancy with no unique cause established. No native benchmark
energy or ALMO comparison is admitted. The600Ry and strict collector are now
terminal as recorded below.

## Mixed-Dy terminal update

Job22410 is terminal. PySCF reports convergence, but its final returned state
has gradient norm1.0670606603907637e-7 versus the required1e-7. All other
native density checks pass; strict admission remains failed, and the archived
energy comparison was not reached. The final orbitals/checkpoint are saved.
No thresholdchange or automatic rerun was launched.

## Native600Ry and strict comparison terminal update

Job22344 hit its declared3500-second native cap during first-water SCF.
The last completed step246 had residual7e-8 versus1e-10, energy
-76.5599956682Eh and change+8.38e-13Eh. Native elapsed time was3501.126seconds,
return code15; Slurm records3507seconds on8CPUs. No second fragment, ALMO
states or terminal integrated-charge summary were obtained. This is a timeout,
not the400Ry iteration-limit abort.

Collector22439 correctly rejected native return codes1/15 before comparing
energies, using zero new native calls and7seconds on1CPU. See its
[result](../../workspaces/lady_almo_implementation_20261023/cp2k_saved_comparison_v1/RESULT.json).
The benchmark pair and collector used48,175 allocated CPU-seconds. Neither
agreement nor disagreement between native and custom ALMO has been established.

Recommendation: preserve the implemented solver and passing core tests, and
resolve the independent-reference setup before applying the solver to infer
La/Dy charge transfer. Do not extend these identical native runs or admit the
quartet from energy stability alone.
