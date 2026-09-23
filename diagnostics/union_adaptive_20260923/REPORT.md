# Union membership and adaptive accommodation: completed common8 pilot

The combined method restores the remaining wrong raw ordering in this consumed
development panel. The minimum La-source minus Ca-source contrast changes from
−4.290011 with the original adaptive pool to **+8.238885 model kcal/mol**.
All15 La/Ca source pairs are correctly ordered, compared with14/15 for adaptive
alone,13/15 for static union and12/15 for released static. These are5La×3Ca
structural comparisons across6biological groups, not15 independent observations.

No reference was fitted to these eight sources. Under the unchanged adaptive
bands, calls change7correct/1wrong→8correct. Under released bands the combined
result is7correct/1inconclusive; under union bands it is6correct/2inconclusive.
Those are explicitly **transferred-band diagnostics**, not a calibrated new
classifier or broad affinity validation. Production/defaults remain unchanged.

## Actual paired outcomes

Raw R = E_Ca − E_La. Values below add a common +405000 display offset for
readability; machine artifacts retain unrounded raw components/values. No aquo
reference or physical free-energy zero is supplied. Larger values are La-like.

| Source | Known class | Released static | Union static | Original adaptive | Union + adaptive |
|---|---|---:|---:|---:|---:|
| 1H4I | Ca | -506.647982 | -506.647982 | -505.420051 | -505.420107 |
| 4MAE | La | -441.666967 | -441.666967 | -444.304425 | -444.304425 |
| Q9Z4J7 | Ca | -464.187748 | -467.751090 | -459.015508 | -464.235878 |
| Q88JH5 | Ca | -469.110051 | -468.588303 | -464.122696 | -463.709057 |
| A0A3 Ca1 | La | -464.654353 | -455.914149 | -463.305519 | -450.716602 |
| A0A3 Ca3 | La | -461.993448 | -462.450584 | -455.103872 | -455.470172 |
| A0AC Ca4 | La | -472.672136 | -472.672136 | -454.844127 | -454.844132 |
| A0AC La4 | La | -461.470592 | -461.470592 | -439.808385 | -439.808385 |

A0A3F2YLY8 Ca1 benefits from both ingredients: its raw contrast rises8.740205
from released to static union, then5.197547 from union accommodation. Its pair
with Ca3 narrows from8.201646 in original adaptive to4.753570. The A0AC Ca4/La4
spread stays15.035748 versus15.035742 previously; adaptive already fixed that
hard source, and the exact context is unchanged here. Reduced spread alone is
not the criterion: released A0A3 has a smaller spread but wrong classifications.

1H4I and4MAE retain their chemistry/context and replay prior adaptive contrasts
within0.000057kcal/mol. Q9/Q88 use their already prepared union contexts and
remain Ca-side controls. Five of16 optimizer endpoints reach a trust-domain
boundary. Proposals are finite constrained candidates, not proven unconstrained
or composite-energy minima.

## Endpoint work and component accounting

The common pool is {origin,adaptive_Ca,adaptive_La}, scored by both metals:

E_M = E_OMOL,vac,M + E_GFN2,ALPBwater,M − E_GFN2,vac,M.

The operational per-metal selection follows the existing mathematical minimum
and0.1kcal origin tolerance; no cross-source absolute energy minima. Native eV
and GFN Hartree convert exactly once. Selected work satisfies
ΔR = work_Ca − work_La. All scalar components and selected candidates are retained
in the comparison artifact. For A0A3 Ca1, work_Ca=−5.402320 and work_La=−10.599867;
the corresponding native/solvent components are−10.282027/+4.879707 and
−13.960416/+3.360549kcal/mol. No extra scalar, entropy, solvent-force optimization
or scaffold FF contribution is added.

## Execution and cost

All16 searches yielded admitted proposals, all16 new cross-metal MACE calls and
all64 new GFN2 single points completed; all8 pools are available. The16 exact
union q0 force arrays and q0 components were reused. No new DFT or scientific
retry occurred. Two solver/collector requests were cancelled while pending
(zero allocation/calls) solely to move unchanged tasks to an available host.

| Stage | Job | Allocated wall s | CPUs | New work |
|---|---:|---:|---:|---|
| Adaptive proposals |1210402|107|32|295 MACE energy/force evaluations|
| Cross-metal MACE |1210423|20|32|16 MACE single points|
| Native GFN2 |1210431|172|64|64 single points|
| Final collection |1210432|3|1|No chemistry|

Total new allocation: **15075 core-seconds and 127 requested GPU-seconds**.
CPU stages requested128GiB (solver) and8GiB (collector); GPU stages requested
200000MiB host memory/oneH200. Proposal worker wall90.279s, peak CUDA6.056GB and
reported host RSS1,850,704KiB. Cross-worker wall13.951s. These are development
allocations, not a folding/preparation-inclusive production latency. Historical
reused costs and local unmetered preparation/testing are additional.

## Scope, limitations and next decision

The fragment union was fixed from all ten saved folds per reference protein,
with source graph closure and unchanged states. This has not established an
equivalent practical union policy for future three-fold scans. Context membership
changes composition and solvent cavity together; the result does not identify a
unique second-shell cause. Labels were not used for membership, mode selection
or endpoint geometry choice. All sources were consumed development material.

The improvement justifies the separately authorized canonical25 plus three
consumed-crystal check using this exact rule and a distinct calibration. That
reference must be frozen before any further transfer evaluation. No225 expansion
or production promotion follows automatically.

## Reproducibility

- Frozen scientific scope: [PLAN.md](PLAN.md).
- Runnable validators/comparison: [COMMANDS.md](COMMANDS.md).
- Machine summary: [RESULT.json](RESULT.json).
- Full coordinates/components: `workspaces/union_adaptive_20260923/pool_v1/collection_final.json`.
- Full paired analysis: `workspaces/union_adaptive_20260923/COMPARISON_v1.json`.
- Actual allocations: [SACCT_v1.txt](SACCT_v1.txt).
- Focused real-artifact tests: [TESTS_v3.txt](TESTS_v3.txt); no fabricated scientific results.
