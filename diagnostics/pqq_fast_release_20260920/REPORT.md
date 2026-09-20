# Fast PQQ standard mode released — 2026-09-20

**The source-to-score release test passed: 25/25 canonical PQQ calls and 3/3
consumed crystal controls.** All 28 raw sources regenerated exactly the archived
complete-context inputs, including every H coordinate, atom order and charge.
All energies were recomputed. A separate literal standard-CLI run also passed.

The promoted mode is `fast_PQQ_OMOL_GFN2_ALPB_v1`, selected by default for
compatible explicit PQQ source requests through `affordable_workflow.py standard`.
The explicit `--mode dft-reference` and existing `baseline` workflow remain
available. Other chemistry retains its current water-prepared DFT workflow.
No raw inbox watcher, existing job, historical score or PLM result was changed.

## What was actually tested

| Check | Result |
|---|---|
| Fresh source normalization + seeded protonation + fixed core + complete polar context | 28/28 exact source/context reproductions |
| Fixed PQQ calibration, unchanged bands | 25/25, 14 Ca-associated and 11 La-associated |
| Previously consumed crystal controls | 1H4I Ca, 4MAE La, 1KB0 Ca: 3/3 |
| Fresh energy-only MACE versus archived native scalar | All 56 endpoints identical |
| Fresh composite versus archived composite | Maximum difference 0.0004057501 model-kcal/mol, Q60AR6 |
| Frozen composite execution repeat target | Passed: maximum below original 0.02 model-kcal/mol target |
| Exact band edges | Both edge references retain their calls; no padding or refit |
| Literal public standard entrypoint | Fresh 1H4I source → Ca-supported; exact context and repeated score |
| Explicit DFT dry route | 56 original-core XYZ/native input files preserved; no new DFT execution |
| Focused source, dispatch, baseline, composite, qualification and graph tests | 57 passed, zero skips |

The 28 structures represent 26 recorded accession groups; crystal replicas are
not extra independent biological observations. All cases were already consumed.
These are PQQ functional-association labels, not interchangeable La/Ca affinity
measurements. This release establishes an operational faster route retaining
reference fidelity, not broader biological validation or improved reference
accuracy over the DFT method that already separated the panel.

## Method and scope

The preparation reuses the actual versioned source normalizers, fixed seeded
standard-only protonator, canonical PQQ(3−) core carver and complete 3.5 Å polar
neighbor graph. Direct-donor checks retain the established 3.3 Å rule. Assembly,
source selectors and five homologous core roles are explicit. Dry water inventory,
protonation, heavy coordinates, caps and charge policy remain fixed. No geometry
optimization or label-based source selection occurs. New source protocol ID:
`fixed_core_PQQ_source_to_complete_context_v1`.

The scoring expression and prior calibrated bands are unchanged:

```
E_comp,M = E_native_OMOL,vac,M + E_native_GFN2,ALPBWater,M − E_native_GFN2,vac,M
R = E_comp,Ca − E_comp,La
Ca-supported: R <= -405464.18774828 model-kcal/mol
La-supported: R >= -405459.1113819997 model-kcal/mol
otherwise: inconclusive
```

The eV and Hartree components are each converted once. There is no compatible
aquo reference for this composite, so no composite S or zero-based affinity
interpretation is supplied. Native OMOL uses the exact prior checkpoint and
float64; energy-only evaluation matches the archived scalar with forces. Native
GFN2 uses the frozen primary ORCA 6.1.1 native mixer, 300 K, NoAutostart and matched
vacuum/ALPBWater inputs. This is a compact-context solvent approximation, not
whole-protein polarization or CPCM added to DFT.

Compatible unseen sources can be scored but explicitly retain
`unvalidated_input_domain`; high folding confidence alone does not certify local
coordination geometry. There is no new geometric cutoff and no prediction-derived
label. The source request must identify the canonical core roles and physical
state. Unsupported donors/cofactors/source gaps, mismatched preparation, changed
method pins and failed components remain visible; no success substitution occurs.

The separate native TightSCF check passed six consumed contexts with maximum
score-correction change 0.025545 kcal/mol under its predeclared 0.20 tolerance.
Independent ordinary-SCF paired agreement remains unavailable. Those are separate
qualification dimensions; this release changes neither solver nor bands. See
[qualification report](../compact_qualification_20260920/REPORT.md).

## Actual cost

| Job | Scope | Wall | Allocation |
|---|---|---:|---:|
| 1203463 | 28 fresh source preparations +56 MACE +112 GFN2 | 1185 s | 37920 core-s; 1185 GPU-s |
| 1203729 | Literal CLI 1H4I: one fresh preparation +2 MACE +4 GFN2 | 44 s | 1408 core-s; 44 GPU-s |
| Total release | 29 preparations, 58 MACE, 116 GFN2; zero DFT | 1229 s | **39328 core-s; 1229 GPU-s** |

Both used one H200, 32 CPUs/MPI slots and 200000 MiB requested host memory.
The full panel spent 305.472 s in source preparation and 874.000 s in prepared
scoring/collection. Median per-site source-plus-endpoint path was 41.104 s;
the full allocation averaged 42.321 s/site. Literal standard CLI timing was
42.386 s from source preparation through score execution, excluding batch startup
and final separate collection. Folding is excluded from every figure.

Peak observed MACE allocated CUDA memory was 4071723008 bytes (3.792 GiB);
peak worker RSS was 1824616 KiB (1.740 GiB). Neither that worker RSS nor batch RSS
is a simultaneous whole-job memory maximum. Allocation core/GPU seconds include
idle allocation and source work; they are not measured hardware utilization.
No matched fresh DFT timing experiment was run here, so no new exact DFT speedup
factor is claimed.

## Artifacts, checks and use

- [Commands and real source example](COMMANDS.md); [release parameters](../../params/pqq_fast_v1.json).
- [Actual results and cost receipts](RESULT.json), [full panel](RELEASE_PANEL.json),
  [finite source manifest](SOURCES.json), [plan](PLAN.md), [literal CLI addendum](CLI_CHECK_PLAN.md).
- Fresh molecular products: `workspaces/pqq_fast_release_20260920/full_v1/`
  and `standard_1H4I_v1/`; old artifacts were not overwritten.
- [57-test log](COMBINED_TESTS_v1.txt). Parser/dispatch/graph assertions reuse real
  outputs; the two actual scientific executions are counted separately above.
- The first core-only implementation replay exposed precision mismatch with the
  old canonical six-decimal XYZ writer; using that exact writer fixed all 28
  before any fresh scientific execution. Both replay records remain saved.
- Initial local test log `RELEASE_TESTS_v1.txt` contains one failed assertion that
  assumed 1e−7 composite equality. The method never promised that tolerance.
  The corrected assertion uses the previously frozen 0.02 repeat target;
  the real Q60AR6 discrepancy remains reported and the bands never changed.

**Recommendation:** use the named fast PQQ standard mode for supported inputs,
retain DFT as an explicit reference, and assess structural robustness of new folds
separately. The source reconstruction and classification release work is complete;
no further calculation is needed to use the documented route.
