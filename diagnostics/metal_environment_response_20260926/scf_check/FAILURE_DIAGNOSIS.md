# Stricter stopping did not remove the rotation discrepancy

**Recommendation: close this SCF remedy; do not launch another SCF algorithm
comparison now.** The additional calculations changed energies much less than
the persistent rigid-transform discrepancy. They do not support spending on
NoSOSCF/DIIS merely to satisfy the density-summary parser. All four cells remain
invalid under their frozen admission rules, and the original failed gates stay
failed. No historical collection or executable was changed for this diagnosis.

## What the actual strict run did

Job 1219450 completed all four native outputs. Contrary to an interpretation
based only on their final summary, none ended in ordinary SOSCF: **all four
switched automatically to TRAH** and completed its macro/micro iterations.
The actual switch message reports insufficient decrease of the maximum gradient
error over ten iterations. It is not a final stop caused by ten unchanged
energies. Ca endpoints report 57 total cycles; La A reports 68 and La rigid 65.

The final-summary density/rotation values cannot be assigned to the final TRAH
stopping criterion from these outputs alone. They also do not match the last
printed SOSCF density steps, so there is no verified exact stale-value mapping.
For example, La A's last printed SOSCF iteration 52 has RMSDP 8.35e−7 and MaxDP
6.59e−6, while its final summary prints 1.3874e−3 and 1.0769e−2. The final TRAH
macro orbital residual is instead 5.528534e−7. The origin of those summary
density values remains unverified; it would be wrong to relabel them as direct
measurements of the final TRAH density change or declare the physical
wavefunction invalid solely because of that interpretation.

Actual final TRAH macro orbital residuals:

| Endpoint | Final macro residual |
|---|---:|
| Ca A | 4.690025e−6 |
| Ca rigid | 2.760470e−6 |
| La A | 5.528534e−7 |
| La rigid | 9.941081e−7 |

The parser did exactly what the frozen diagnostic required: it rejected the
printed density/orbital-summary criteria. This remains the recorded result;
reinterpreting a solver-specific summary after inspection cannot retroactively
turn this diagnostic into a passed qualification.

The [ORCA 6.1 SCF manual](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html)
documents automatic TRAH activation when orbital-gradient convergence slows,
and identifies TolG as its gradient-norm convergence threshold. Mode 0 can
allow overachievement exceptions. NoSOSCF alone would not isolate DIIS because
AutoTRAH is a separate default mechanism; a true DIIS-only experiment would
also need to disable that switch. These facts explain why the keyword change
was not equivalent to a separately measured final density-difference guarantee.

## Raw numerical diagnostics — explicitly unqualified

These values are extracted directly from the actual completed outputs/engrad/
pcgrad files. They are **not admitted endpoint results**, do not overwrite
FINAL_COLLECTION.json, and are used only to evaluate whether the attempted
numerical remedy moved the known error.

| Quantity | Ca | La |
|---|---:|---:|
| Raw rigid-minus-A energy, Eh | +1.77566862476e−5 | +1.43655070133e−5 |
| Raw rigid-minus-A energy, kcal/mol | +0.01114248885 | +0.00901449175 |
| Maximum transformed core-gradient residual, Eh/bohr | 4.27911358329e−5 | 5.83932420000e−5 |
| Maximum transformed external-gradient residual, Eh/bohr | 4.33293765392e−7 | 6.01190357398e−7 |
| A energy change versus preceding refined run, kcal/mol | −8.91259488e−6 | −1.31287633e−5 |
| Rigid energy change versus preceding refined run, kcal/mol | −5.47525475e−5 | −9.90521362e−6 |

Both raw rigid energy differences still exceed 1e−5 Eh. The gradient residuals
are below 1e−4 Eh/bohr, but this is not sufficient to pass the failed combined
gate. The preceding refined-grid rigid differences were 1.7829736862e−5 Eh
(Ca) and 1.4360369960e−5 Eh (La); the extra SCF work barely changed them.

This rejects the attempted stopping-policy change as an effective remedy for
the observed rotational energy error. It does **not** uniquely distinguish
remaining quadrature, screened-integral, grid-construction, or other numerical
approximations. The small changes are evidence against the original early-stop
hypothesis being the dominant error at this precision, not proof of universal
SCF convergence or a license to relax tolerances.

## Cheapest useful next evidence

The cheapest decisive evidence is already available: compare the above raw
shifts with the previous refined run and terminate this numerical branch.
No new NoSOSCF/DIIS jobs are justified by its result. A future independently
motivated region-sensitivity diagnostic may report coarse changes beside the
known numerical failures; it must not inherit force qualification, authorize
optimization, or interpret effects comparable to the unresolved numerical error.
No B endpoint was run under this SCF policy and no strict-policy environmental
response or affinity is available.

## Exact source artifacts

All paths below are relative to repository root:

- `workspaces/metal_environment_response_20260926/scf_check_v1/FINAL_COLLECTION.json`
  retains the failed frozen qualification and original receipt/artifact hashes.
- Its `Ca_A_strict`, `Ca_rigid_strict`, `La_A_strict`, and `La_rigid_strict`
  directories contain the actual `endpoint.out`, `endpoint.engrad`, and
  `endpoint.runtime.pcgrad` files used here.
- `workspaces/metal_environment_response_20260926/grid_check_v1/FINAL_COLLECTION.json`
  supplies the preceding refined endpoints and actual energies.

This investigation ran only file parsing/algebra and read the official manual.
No new molecular calculation, parser amendment, or submission occurred.
