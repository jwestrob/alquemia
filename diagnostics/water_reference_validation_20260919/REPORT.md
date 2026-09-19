# Contextual water preparation transfers to parvalbumin

**Complete: two MACE proposals and two original-recipe DFT single points; no
failed calculations or retries.** The ordered carp parvalbumin 4CPV CD→EF vector
is preserved. Dry CD is an exact identity. Preparing the actual EF water A166
changes `R = E_Ca − E_La` by **−3.3215921675 kcal/mol**, toward Ca on the existing
raw contrast. This demonstrates transfer to another real water-bearing protein;
the evidence does not support an accuracy claim.

## Actual results

| Site | Source core waters | Baseline R | Prepared R | Change in R |
|---|---|---:|---:|---:|
| CD | none | −405428.555612006 | −405428.555612006 | 0 exactly |
| EF | A166 | −405401.4468577985 | −405404.768449966 | −3.3215921675 |

All contrasts are kcal/mol. The large common elemental energy offset is retained;
it cancels in the preparation difference. No new aquo reference or decision band
was fitted or borrowed. Calibrated decisions and occupancy probabilities remain
explicitly null.

| EF endpoint | Original energy / hartree | Prepared energy / hartree | Change / kcal mol−1 |
|---|---:|---:|---:|
| Ca | −1877.263970774381 | −1877.331181292529 | −42.1752368903 |
| La | −1231.215642175284 | −1231.277559399378 | −38.8536447228 |

The inherited preparation **includes water internal-geometry normalization** before
rigid orientation: source O–H distances 1.1863966453 and 1.1755024458 Å become the
established reference distances 0.9571628483 Å. These large energy drops combine
that repair and orientation; they are not an isolated measurement of orientation,
entropy, or increased predictive accuracy. No extra normalization was invented for
parvalbumin.

## Method and invariants

Protocol: `independent_reference_contextual_water_H_transfer_v1`. Reuse the exact
native OMOL100M optimizer from the demonstrated contextual-water manifest, source
and radial-away starts, float64, BFGS with gradient tolerance 0.001 eV/rad and
maximum 200 iterations. Both starts converged for both endpoints; source was
selected in both. There were 20 Ca and 23 La recorded objective evaluations.

The 44-atom source-backed context contains variable water A166 and frozen source
outer water A257. Nearby supported polar groups include Asp92, Asp94, Lys96
backbone oxygen and Glu101. The 3.5 Å context rule, source geometry, charge ledger,
water geometry and optimizer settings were frozen before scoring. The other site
and its remote metal are outside this local model; this is not a whole-protein
calculation. No nearby cofactor was silently omitted.

Only the two H coordinates of A166 transfer to the original amide-v3 EF core.
Water oxygen, every nonwater atom, caps, atom order, water count, protonation,
charges and multiplicities are unchanged. DFT uses the actual original input body:

```text
! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3
```

ORCA 6.1.1; no additional composite correction, TightSCF, gradient or optimization.
Original Ca/La endpoints and preparation records remain intact. The experiment
does not modify the default workflow; the separately authorized integration is
owned by the parent agent.

## Evidence and interpretation

The existing [evidence ledger](../benchmark_set_20260915/EVIDENCE.json) marks both
`CARP_PV_CD` and `CARP_PV_EF` directions **unresolved**, with La affinities and
cross-study Ca support. They are one protein group, retained as an ordered vector,
not two independent matched-assay La/Ca labels. The frozen source is
DOI 10.1021/bi00294a030 / PMID 6661415. No new labels were inferred here.

Numerical implementation and state preservation pass the checks below. Practical
transfer to a second water-bearing protein is demonstrated, and preparation does
not systematically push every score toward La. Biological improvement is still
unresolved for this protein. Lower endpoint energy alone does not establish a
better discriminator.

Aequorin was inventoried but not run: its existing whole-protein prepared objects
include a chromophore, background metals and waters from multiple sites. They
require an explicit coherent multisite/cofactor adapter; changing a provenance
field cannot make them equivalent to a source-backed single-site repair manifest.
No artificial waters or favorable site selection were used.

## Execution, cost and tests

| Job | Work | Status | Scheduler wall / s | Allocated CPUs | Allocated core-s | GPU-s |
|---|---|---|---:|---:|---:|---:|
| 1202475 | two native MACE water-H searches, H200 | completed | 24 | 16 | 384 | 24 |
| 1202476 | two concurrent 16-rank native DFT SPs | completed | 70 | 64 | 4480 | 0 |
| Total | | | | | **4864** | **24** |

Slurm rounded the DFT request of 32 CPUs to **64 allocated CPUs**; the table counts
the real allocation. The DFT runner reports 68.5828063 s and 4389.2996037 core-s
inside the allocation; scheduler accounting includes setup. DFT batch MaxRSS was
2,845,676 KiB. Peak MACE host RSS was 1,710,808 KiB and peak CUDA allocation
1,229,083,136 bytes. Local preparation, collection and tests were not metered.
Different archived hardware prevents a controlled production speed comparison.

**Six real-fixture replay tests pass, zero skips, 0.289 s.** They check ordered
site identity, source waters and frozen context, paired state/coordinate invariants,
exact original DFT recipe, actual MACE/DFT receipt replay, sign and unit algebra,
and the preserved pre-execution missing-result status. These are distinct from
the four scientific tasks actually executed above.

## Reproducible artifacts

- [Frozen config](CONFIG.json), [agreement](AGREEMENT.md), [commands](COMMANDS.md).
- Proposal manifest:
  `workspaces/water_reference_validation_20260919/proposals_v1/manifest.json`,
  SHA256 `21ad25f652d4e6179d3e495f0b257249ae0444ef94e67fdf9feec9efc54a6f70`.
- Actual proposals: `proposals_v1/collection_job_1202475.json` in that workspace.
- DFT manifest: `workspaces/water_reference_validation_20260919/dft_v1/manifest.json`,
  SHA256 `aa408a2b69f3f606dce2eb78f366b9c48607a36f19440969f1ccdcd9c9ec37bd`.
- Actual final paired results: `dft_v1/collection_1202476.json`.
- [Scheduler accounting](COSTS_sacct.tsv), [final tests](TESTS_FINAL.txt),
  [MACE submission](MACE_SUBMISSION.json), [DFT submission](DFT_SUBMISSION.json).

Recommendation: use this as a successful independent **preparation transfer**
check for the authorized integration, preserving original/prepared outputs and
explicit unsupported cases. Do not claim a new parvalbumin accuracy result or
occupancy model. No further parvalbumin calculations are needed for this test.
