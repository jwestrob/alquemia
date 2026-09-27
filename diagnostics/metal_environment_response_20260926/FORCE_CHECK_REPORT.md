# Native force checks: local derivatives pass; rigid gate fails

All 20 real endpoints completed in job1219319 (455s,344 allocated CPUs,
156520 allocated CPU-s; no GPU). Collector1219320 took9s on1CPU. The completion
watcher automatically resumed root, acknowledged in force_checks_v1/ROOT_ACK.json.

14/16 checks pass: both metals pass two step sizes of external hydroxyl and
physical cap-boundary directional derivatives, step convergence, and exact repeats.
Both joint rigid checks fail the frozen gate: Ca energy shift0.012791905 and La
0.011360295 kcal/mol exceed1e-5Eh; La core-gradient residual0.000104952Eh/bohr
also exceeds1e-4. Original qualification remains not_qualified.

Printed D4 and gCP terms are unchanged; the SCF component accounts for the
rigid energy shifts. This suggests quadrature error but does not prove its cause.
A separate eight-cell grid diagnostic separates translation from rotation and
checks a denser unpruned quadrature without changing the chemical Hamiltonian.
Selected force checks qualify only this embedded electronic component; missing
full hybrid cross interactions and inaccessible exact ML weights remain blockers.
See RESULT_1219319.json and the immutable full collection for actual residuals.
