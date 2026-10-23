# The larger water-shell error is in the collective response

21 October 2026. Fixed-geometry diagnostic; independent actual-output review passes 402 checks.

**The model underestimates the additional whole-shell response by 5.60 kcal/mol
on this La/Dy contrast. Its isolated-pair error is smaller and has the opposite
sign. Changing pair interactions alone could mask the collective error.**

## Result

We compare the change from the saved eight-water arrangement to the saved
nine-water arrangement, for Dy minus La. Each isolated metal–water fragment keeps
its literal cluster coordinates. The original cluster energies are reused.

| Contribution | ECP-PBE0-D4 reference | Frozen polarizable model | Model minus reference |
|---|---:|---:|---:|
| Sum of individual metal–water interactions, P | −8.933375 | −7.670058 | +1.263316 |
| Additional collective response, B | +8.399168 | +2.797568 | −5.601600 |
| Full-shell contrast, C = P + B | −0.534207 | −4.872491 | −4.338284 |

All values are kcal/mol. The pair and collective errors partly cancel. Both
remain material relative to the declared 1 kcal/mol energy resolution.

The calculation includes 34 literal reference dimers and two own-basis ions,
34 matching model dimers, and four reused clusters per method. The ion energy
difference is subtracted once from the difference of dimer sums; it cannot be
discarded because the coordination counts differ by one. Water-only terms cancel
between metals at each size. Independent decimal arithmetic from raw energy
tokens reproduces all 17 matched water contributions and the complete identity.

## Mechanical interpretation

Within this model, the collective remainder is polarization at printed precision.
Its additive van der Waals contribution is identical in the cluster and extracted
pair sums. The reference remainder consists of +8.476002 kcal/mol electronic
response and −0.076834 kcal/mol D4 dispersion nonadditivity.

This localizes the larger discrepancy to the difference between isolated pairs
and the assembled shell. It supports investigating collective electrical response
before another radius or epsilon fit. It does not uniquely identify a missing
polarization, charge-transfer or short-range term: the quantum remainder also
contains finite-basis and other electronic-model effects.

The separate metal–whole-shell counterpoise test changed the total error from
−4.338284 to −4.663006 kcal/mol. That particular basis-sharing correction does
not explain the failure. It did not establish the isolated-pair basis limit, so
it is not inserted into this own-basis decomposition. Original fitting targets,
parameters and withheld validation remain unchanged.

## Limits and admission

- These two arrangements differ in more than the presence of a ninth water.
  This is a coordination-and-arrangement contrast, not an isolated ninth-water
  binding energy or a bulk exchange free energy.
- All 36 new reference jobs stop under the native energy-check route. Their
  printed maximum/RMS density and DIIS criteria remain false. They pass the
  declared collection protocol, not a claim that every printed SCF criterion
  passes. Numerical, frozen-core and electronic-state limitations remain.
- This energy-only test does not qualify forces, sampling, protein transfer or
  independent physical accuracy. No optimization or parameter change occurred.

## Evidence and cost

The [collected result](../../../../workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/matched_pair_decomposition_v1/collection/run_v1/result/RESULT.json)
retains raw endpoint provenance, electronic/D4 splits and every native stopping
record. The [independent review](../../../../workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/matched_pair_decomposition_v1/independent_review_v1/RESULT_REVIEW.json)
passes 402 checks. An initial review-script metadata-key failure is preserved;
correcting that lookup changed no scientific source or result.

Scheduler cost is 4,948 allocated CPU-seconds for the 36 reference jobs, 36 for
model job 17819 and 7 for collector 17822: **4,991 allocated CPU-seconds** total.
The reference native wall-time-times-rank proxy is 4,922.4005 CPU-seconds and is
kept separate from scheduler allocation cost and measured utilization. The
independent review made zero molecular calls. Preparation and review receipts
retain their own local costs.

This finite diagnostic is complete. The next decision is a targeted check of the
collective-response discrepancy and reference uncertainty; no new molecular
experiment or model correction is released by this report.
