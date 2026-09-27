# Three real cores: environmental response measured, ML qualification unavailable

All18declared primary reference endpoints completed: Ca/La x A/B/isolated for
1H4I,4MAE,1F6S. This is a working embedded-electronic response harness spanning
PQQ and a consumed non-PQQ source, not evidence of improved classification.

| Conditional model | CaB-A,Eh | LaB-A,Eh | delta=(LaB-LaA)-(CaB-CaA),kcal/mol |
|---|---:|---:|---:|
| 1H4I54 | -0.0053045715499138 | -0.0047423556279682 | 0.352795817455 |
| 4MAE80dry | -0.005706942110919 | -0.0053007094361419 | 0.254914852071 |
| 1F6S52normalized | -0.00017750889787749 | 0.00025118617509179 | 0.269010219745 |

All three selected positive10degreehydroxyl motions have positive differential
responses; no class sign was prescribed.1F6S individually lowersCa electronic
energy and raisesLa electronicenergy. These are distinct sites/motions, not
replicated perturbations and not comparable absoluteprotein preferences.
The controls' evidence types remain separate; this table adds no biologicallabels.

1H4I fullThr159 promotion changesdelta+.352795817 to+.520767257, same direction.
The original rigid-motion qualification failed; denserquadrature and stricterSCF
stopping did not fix it. Original selecteddirectionFD/repeat checks passed.
4MAE/1F6S derivative/rigid/refinement checks are unrun. No numericalerror bound is
inherited by a newcore. ExactMACEPOL-EFweights/interface remain unavailable, so
there are zeroactualMLcomparisons and no validatedsurrogateaccuracy claim.

4MAE conditionally excludes historicalwater/15Padduct;1F6S uses explicitlydistinct
normalized52state with2waters, not40atomrepair.1H4I/4MAE retain longsourceOHbonds;
1F6S OH0.96A and more remote4.64009A hydroxylselection. No perturbation was changed
to rescue a response. No newH/water/minimization sweep or libraryrun follows.

Latest1219489 completed389s/344CPUs=133816allocatedCPU-s,zeroGPU;collector2s/1CPU.
Cumulative molecular workers54outputs/1319928allocatedCPU-s,zeroGPU, including
all numericaldiagnostics. ActualCPUutilization and aggregateMPIpeakRSS unavailable;
these costs do not establish productioninference affordability.

Next: finite32configuration classicalcomponent energy/force checks of exact1H4I
additiveledger, with no newQM/MD/optimization. This addresses assembly and force
accounting. It cannot pass the unresolvedelectronicgate or supply missing solvent.
Complete hybrid/solventrelaxation and predictivebenefit remain unsupported.
