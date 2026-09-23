# Three-source motion envelope: mixed pilot result

**The4.3 Å envelope retains all25 canonical and three crystal calls and repairs both earlier A8 probe abstentions, but introduces one A0A3 Ca-sample3 abstention.** Its six probes are5correct/0wrong/1inconclusive, versus6/0/0 for matched strict ten-fold scoring. This supports wider qualification of the practical three-source preparation; it does not establish equivalent fidelity or justify promotion.

All34 pools are available. Expected classes and source membership are unchanged. These are already consumed development references and fold probes, not independent proteins or direct affinity validation.

## Separate preparation and search results

| Method, with its own canonical-only reference | Canonical | Crystals | Six probes |
|---|---:|---:|---|
| Strict ten-fold origins |25/25|3/3|3correct,1wrong,2inconclusive|
|4.3 Å envelope origins|25/25|3/3|2correct,1wrong,3inconclusive|
| Strict ten-fold adaptive pool |25/25|3/3|6correct|
|4.3 Å envelope adaptive pool|25/25|3/3|5correct,1inconclusive|

The envelope search repairs three of its own origin probe calls. Preparation alone does not solve this panel. Its adaptive class gap is7.071396249113604 modelkcal/mol; the separate origin-only gap is10.117019388068002. A larger origin gap is not an accuracy gain. Mathematical and operational adaptive choices give the same calls here. Through the unchanged strict ten-fold adaptive bands, envelope scores give31correct/3inconclusive, including two canonical abstentions: the new representation requires its own reference.

Only the25 designated canonical geometries fit either reference, using the existing extrema/minimum-gap rule. Crystals and probes are excluded. The separately authorized origin ablation was added after origin completion and remains distinct from the adaptive calibration.

## Prespecified probes

| Source | Expected | Envelope adaptive | Margin to own class edge, kcal/mol | Origin call |
|---|---|---|---:|---|
|a0a3f2yly8 Ca sample1|La|correct|4.305488|inconclusive|
|a0a3f2yly8 Ca sample3|La|inconclusive|-0.569343|inconclusive|
|a0acd6b9f2 Ca sample4|La|correct|2.867545|wrong|
|a0acd6b9f2 La sample4|La|correct|16.864168|inconclusive|
|a8r3s4 La sample1|Ca|correct|4.700190|correct|
|a8r3s4 La sample3|Ca|correct|5.272498|correct|

## Why A0A3 Ca-sample3 becomes inconclusive

The matched original core coordinates agree exactly, and both contexts retain Ca/La charges−3/−2. The envelope has202 atoms versus193. Its origin contrast shifts−2.960946848kcal/mol: nativeMACE−5.287715162, solvent subtraction+2.326768314. Differential accommodation improves+1.399701620, leaving a final−1.561245228 shift from the strict ten-fold result. The resulting score is0.569342826 below its new La edge. This is mainly a starting-context offset, partly repaired by search; it is not an excess accommodation penalty. Composition and continuum cavity change together.

The unchanged four-mode selector keeps197χ1,197χ3 and329χ2, but replaces197χ2 with327χ1. Neither endpoint is boundary limited: new maximum heavy displacement is0.741058Å(Ca)/0.638445Å(La), versus0.765690/0.674366Å before. No Ca3-specific optimization, radius change, threshold adjustment or member replacement was performed.

Probe-pair spreads are mixed. A0A3 changes4.752626→4.874830kcal/mol; A0AC changes15.032615→13.996623; A8 changes0.421140→0.572308. This is not a uniform structural-robustness improvement.

## Actual scope, integrity and cost

Preparation retains94complete/6unavailable declared triples and prepares135/135 unique source-context pairs. Only the frozen34-source pilot was scored. All68 searches succeed;25 final proposals touch a bound, with no unconstrained-minimum claim. The existing SLSQP kernel also records175 evaluated trial points outside final admissibility; none is selected as a final proposal. The4.3 Å reachability argument applies to admitted anchor movements against omitted fixed atoms, not every optimizer trial or continuum-solvent convergence.

Actual new work:68 MACE origins +998 search evaluations +68 cross evaluations =1,134 MACE calls;136 origin +272 candidate nativeGFN2 cells =408, all accepted. Scalar inputs use the separately qualified nativeTolE1e−10,300K,MaxIter500,freshNoAutostart,one-rank profile. No loose scalar result is reused. No DFT, folding, new protonation/waters or production change.

The five outer allocations used67,040 CPU-seconds and1,250 requestedGPU-seconds. Their elapsed times are20s(originMACE),1,178s(search),629s(originGFN2),52s(crossMACE),216s(candidateGFN2). They overlap; their sum is not end-to-end latency.

Search executor wall is1,090.311s; summed per-search wall862.759s; summed actual MACE-evaluation wall246.640s; one-time model-load time1.556s. The executor minus recorded GPU-worker interval is210.446s, retained as unattributed overhead. Model construction occurs once per batch; remaining initialization/collection/setup work is not assigned an unmeasured cause or assumed fully reusable. The origin/candidate scalar runners report440.028/170.007s inside their allocations. Local preparation took73.108s; local preflight/report work is additional. This34-source batch does not measure isolated three-source latency.

## Reproducibility and limitations

The reference writer was archived byte-for-byte before correcting a report-only component-key alias. REFERENCE_v1 remains intact; REFERENCE_PINNED_v1 preserves every scientific field and its frozen UTC while pointing to that exact archived implementation. Singleton crystals explicitly have no declared triple ID. An earlier local terminal-file watcher exited before the completed collection became visible; direct preparation used the final artifact. Neither reporting repair nor watcher retry made molecular calls.

All18 focused actual-fixture tests pass across preparation, execution preflight, source joining, pool integrity and final algebra. Raw matrices, gradients, traces, scalar diagnostics and failed preliminary report logs remain in the workspace. The full100 transfer has not been launched; a separate finite inventory is the next step. Keep the Ca-conditioned probe regression visible in any later La-conditioned transfer result.

Compact artifacts: [RESULT](RESULT.json), [probe table](PROBES.csv), [costs](COSTS_v1.json), [commands](COMMANDS.md). Full numerical artifacts are under `workspaces/motion_envelope_20260923/`.
