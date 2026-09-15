# GPT-6 Pro scientific review brief: La/Ca energetic discriminator

## Assignment

Audit whether the current frozen-geometry DFT strategy can demonstrate
lanthanide-versus-calcium information beyond the canonical PQQ-MDH D+2-Asp,
fragment-charge, and atom-count confound. Focus on scientific validity and
experimental identifiability, not software style or cluster engineering.

No energies from the 2026-09-15 benchmark expansion had been inspected when
its design and inputs were frozen. Do not use a desired outcome to revise a
label, select a site, add a water, alter a carve, or move a threshold.

## Read first, in this order

1. `diagnostics/pqq_pmdh_fixed_core_calibration_20260914/EXPERIMENT.md`
2. `diagnostics/pqq_pmdh_fixed_core_calibration_20260914/RESULT.md`
3. `diagnostics/pqq_pmdh_fixed_core_calibration_20260914/reserved_crystal_holdout/result/HOLDOUT_RESULT.md`
4. `diagnostics/pqq_boundary_pair_20260914/README.md`
5. `diagnostics/pqq_boundary_pair_20260914/RESULT.md`
6. `diagnostics/pqq_balanced_embedding_20260914/EXPERIMENT.md`
7. `diagnostics/pqq_balanced_embedding_20260914/RESULT.md`
8. `diagnostics/laca_benchmark_expansion_20260915/MASTER_EXPERIMENT.md`
9. `diagnostics/laca_benchmark_expansion_20260915/benchmark_manifest.tsv`
10. `diagnostics/pqq_q46444_1kb0_external_validation_20260915/EXPERIMENT.md`
11. `diagnostics/pqq_q46444_1kb0_external_validation_20260915/PREPARATION_AUDIT.md`
12. `diagnostics/nonpqq_direct_site_benchmark_20260915/PREREGISTRATION.md`

Consult the lane manifests and preparation reports only when needed to audit a
specific assertion. The implementation is useful evidence of what was actually
done, but the requested judgment is scientific.

## Facts that must remain distinct

- The 25-member canonical PQQ panel separates perfectly, but biological class,
  D+2 Asp, fixed-core charge, and fixed-core atom count are perfectly coupled.
- The preregistered crystal transfer passed: Ca-supported 1H4I scored 7.52 and
  Ln-supported 4MAE scored 38.09 kcal/mol under fixed-core-v3.
- Expanding a distance carve from 3.3 to 3.6 A moved MxaF in the La direction;
  it did not solve the confound.
- Charge-conserving whole-chain point-charge embedding produced a roughly
  64--68 kcal/mol discontinuity when Asp(-) crossed the QM/MM representation
  boundary. That embedding model was rejected, not rescued post hoc.
- The new campaign separates canonical PQQ transfer from an architecture-neutral
  direct-site lane. The latter has one strong same-assay falsifier at present:
  GGR/MglB 1GLG is about 29-fold Ca-selective and must score below zero.
- Aequorin and parvalbumin are supporting/site-vector evidence, not independent
  scalar calibration labels. Calbindin D9k remains an important planned
  opposite-site-direction test but is not in the currently prepared six-site
  panel.
- The score is a frozen electronic compatibility contrast relative to an aquo
  gauge. It is not a binding free energy, Kd, physiological assignment, or
  within-lanthanide discriminator.

## Questions to answer

1. What, exactly, do the completed calibration and crystal transfer establish,
   and what do they not establish because of label/core-charge leakage?
2. Can the frozen 2026-09-15 external campaign falsify the claim that DFT adds
   information beyond motif and net charge? Identify any remaining circularity,
   pseudoreplication, label ambiguity, or architecture-dependent quantity that
   prevents a clean answer.
3. Is the vertical La(III)/Ca(II) swap against a symmetric CN8 aquo gauge a
   defensible directionality test for these frozen sites? Analyze the unequal
   ionic charge, fixed geometry, coordination number, proton linkage, missing
   entropy/reorganization, CPCM cavity, fragment net charge, and explicit-water
   assumptions. Separate fatal problems from bounded limitations.
4. Is rejecting the balanced point-charge embedding scientifically justified?
   Diagnose the large boundary discontinuity and say whether any electrostatic
   embedding or QM/MM formulation could repair it without making the benchmark
   adaptive after results are seen.
5. Are Q46444/1KB0 and GGR/1GLG genuinely informative external controls under
   their frozen preparation rules? Are the omitted TFB/water choices in 1KB0
   and the dry deposited Ca geometry in 1GLG defensible, or do they change the
   question being tested?
6. Before expensive execution, what is the smallest set of preregistered changes
   or added proteins needed for a difficult-to-narrativize yes/no assessment?
   Prioritize independent same-assay La/Ca direction labels and matched opposing
   outcomes. Do not propose second-shell mutants or within-lanthanide work as a
   substitute for first validating La/Ca discrimination.
7. If the panel succeeds, what claim is supportable? If GGR or the new Ca-class
   PQQ controls fail, what specific methodological conclusion follows and what
   single next model should be tested?

## Requested output

Give a blunt verdict under three headings:

1. **Current evidence** -- strongest supportable claim now.
2. **Benchmark adequacy** -- whether the frozen expansion can answer the stated
   question, with fatal flaws called out explicitly.
3. **Recommended decision** -- run unchanged, amend before running, or abandon
   this discriminator; list only the highest-value amendments in priority order.

Treat a method that merely rediscovers donor count, D+2 Asp, or fragment charge
as a failure to demonstrate incremental energetic information, even if its
classification accuracy is perfect.
