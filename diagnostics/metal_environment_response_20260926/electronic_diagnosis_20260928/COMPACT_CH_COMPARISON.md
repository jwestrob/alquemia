# C–H repair removes strain, not the compact source reversal

**All six repaired endpoints completed. The opposing Hans source contrasts remain.** Hydrogen repair is a real geometry improvement, but this comparison demonstrates no repaired discrimination or robust within-series preference. No further compact calculation follows from this report.

The same frozen-core PBE0-D4 treatment evaluated original and repaired compact geometries. Only admitted real C-bound hydrogens changed; heavy atoms, metals, synthetic caps, exchangeable H and waters stayed fixed. The two source structures are retained separately, not averaged.

## Conditional balanced exchange

`D = (E_Hans,Dy − E_Hans,La) − (E_Mex,Dy − E_Mex,La)`. Positive is more La-favoring in Hans relative to Mex under this fixed-state isolated model. This is not binding free energy or a site-specific experimental affinity label.

| Hans source | Original D | Repaired D | Change |
|---|---:|---:|---:|
| Hans8DQ2 | +11.330886345 | +11.090948176 | -0.239938169 |
| Hans8FNR | -18.468182860 | -18.612513541 | -0.144330681 |

Units kcal/mol. Both signs persist. The original source disagreement is approximately29.80kcal/mol; repaired disagreement remains approximately29.70kcal/mol. No threshold, reference offset or source selection was changed.

## Each metal benefits energetically, but similarly

| Source | La repair work | Dy repair work | Dy−La repair work |
|---|---:|---:|---:|
| Hans8DQ2 | -219.273422 | -219.589175 | -0.315753 |
| Hans8FNR | -168.539870 | -168.760015 | -0.220146 |
| Mex8FNS | -52.259952 | -52.335767 | -0.075815 |

Repair work is repaired minus original electronic energy, kcal/mol. The large decreases support the interpretation that the old C–H geometry was strained. Their near cancellation between metals explains why the selectivity contrast barely changes. These are endpoint energy changes from a classical preparation proposal, not quantum relaxation minima or a free-energy correction.

## Metal loads remain large and source dependent

Directions were selected from heavy geometry before results and are exactly unchanged between old/repaired analyses. Positive gradient means moving the metal toward its nearest real oxygen raises energy; force is the negative.

| Source / metal | Original gradient load | Repaired gradient load |
|---|---:|---:|
| Hans8DQ2_La | 27.930113 | 28.554868 |
| Hans8DQ2_Dy | 19.004141 | 19.008727 |
| Hans8FNR_La | 64.320034 | 63.599668 |
| Hans8FNR_Dy | 50.263668 | 49.343762 |
| Mex8FNS_La | -56.664720 | -56.940824 |
| Mex8FNS_Dy | -56.637209 | -56.981514 |

Units kcal/mol/angstrom. The load change is small beside each endpoint load. Hans favors displacement away from its selected nearest oxygen along this direction, whereas Mex favors moving toward its selected oxygen; these are different source-specific directions, not a common trajectory between proteins. They indicate remaining local strain/surface response, not a preferred metal or an instruction to optimize the contrast.

All-atom capped-system translation residuals and metal-gradient norms are retained in the JSON; no residual was subtracted from forces. No ligand-to-protein force mapping is asserted. The earlier Dy single-direction qualification does not automatically qualify every La/source direction or the full physical electronic state.

## Cost, coverage and next decision

Six of six original and six of six repaired endpoints are available. Worker1220316 completed in1773s on48allocatedCPUs:85,104allocatedCPU-seconds. Collector1220317 adds2; total85,106 for this electronic batch/collection, zeroGPU. The separately executed upstream classical H repair cost32allocatedCPU-seconds, giving85,138 including that preparation. Batch-step MaxRSS was9,705,744KiB; this is not an asserted total over MPI children. Per-endpoint receipt times are in JSON.

No affinity, equilibrium occupancy/populations, broad transfer qualification, or full hybrid force validation is available. Old exchangeable-H/water geometry remains and is explicitly outside this targeted test. The absence of a selectivity improvement does not make the geometry repair unnecessary; it rules out this particular C–H defect as a sufficient explanation for the source reversal.

**Recommendation:** retain the corrected preparation for separately declared future research, preserve both negative discrimination results, and stop this compact repair branch here. The independently running embedded work asks a different question. Production and historical calibrations are unchanged.
