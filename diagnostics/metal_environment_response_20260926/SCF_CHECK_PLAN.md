# Density-convergence diagnosis — declared 27 September 2026

Jacob requested continued work. The denser-grid branch is complete; this new
four-endpoint test targets a distinct, observed stopping-policy issue.
All four refined A/rigid runs stopped with RMS density changes approximately
1.85e-4 to2.61e-4 and maximum changes0.00398 to0.00579, despite TightSCF density
thresholds5e-9 and1e-7. Energy/orbital convergence is not by itself proof that
this caused rotation sensitivity; test that hypothesis directly.

Use actual refined Ca/La A and rigid coordinates and point charges unchanged.
Same r2SCAN-3c, native ECP, singlets, charges, analytic gradients, no continuum;
same AngularGrid7 IntAcc7 unpruned HGridReducedfalse. Keep TightSCF thresholds
and MaxIter125. Add explicit ConvCheckMode0 and ConvForced1, recording all
existing density, energy and orbital thresholds. Official ORCA6.1 SCF manual
explains mode2 tests total/one-electron energy, mode0 checks all criteria, with
possible overachievement exceptions. Therefore parse actual terminal residuals;
do not infer strict density convergence from the success banner alone.
https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html

Four new cells, Ca/La x A/rigid; no B or new protein yet. Rigid tolerances remain
1e-5Eh and1e-4Eh/bohr. Compare electronic energy and gradients to refined originals,
retain failed/unavailable cells and costs. No new scientific result or label may
alter tolerances or state. This test neither repeats a grid sweep nor establishes
refined response/finite-difference qualification. Follow-on requires interpreting
this diagnosis; no automatic increase of125iterations.

Four workers x86MPI ranks=344 on exclusive node, --mem=0, verified full RealMemory,
existing25percent runtime headroom, normal queue priority and no PQQ dependencies.
Actual collector+Codex queue watcher; Monday shutdown pre-start cutoff retained.
