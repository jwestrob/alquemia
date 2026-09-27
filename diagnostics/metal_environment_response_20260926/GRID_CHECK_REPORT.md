# Denser quadrature preserves response, but does not pass rigid qualification

Eight/eight endpoints completed, job1219437,583s on344 allocated CPUs
(200552 allocated CPU-s, zero GPU). Collector1219438 completed6s/1CPU.
The automatic completion event resumed root successfully. Four real-artifact
implementation tests passed before execution. Actual ORCA headers confirm
IntAcc7.000, AngularGrid7 Lebedev770, GridPruning0(no pruning), with explicit
HGridReduced false. No parser repair or tolerance change was needed.

Both original-grid translation checks pass. Both refined joint-rotation checks
still fail the frozen energy tolerance1e-5Eh: Ca1.7829736862e-5Eh and
La1.4360369960e-5Eh. Both refined core/external gradient residuals pass1e-4Eh/bohr.
Thus this refinement does not resolve the original failed qualification.
It does not uniquely identify SCF convergence, quadrature, or another approximation
as the remaining cause. No further grid round is automatically launched.

The useful response is insensitive to this refinement. Ca B-A=-3.32846906155,
La B-A=-2.97582879859 kcal/mol; delta_env_el=+0.352640262954kcal/mol,
a change of -0.000155554501 from the original. Both individual response changes
also pass the predeclared0.05kcal/mol target. This is a modeled hydroxyl
perturbation's embedded electronic response, not affinity or a class decision.

Original force qualification remains failed; refined directional finite
differences were not run. Exact ML weights remain unavailable; missing hybrid
cross interactions remain unsupported. Do not extend to LanM optimization on
these results. Retain this finite test and stop unchanged numerical spending.
