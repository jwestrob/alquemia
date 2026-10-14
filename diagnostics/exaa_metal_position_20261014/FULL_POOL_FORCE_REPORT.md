# Full transfer adds no scoring benefit; one interior probe is justified

**Both full historical scores remain exactly unchanged.** Collector4133 validates
all20 saved cells and reproduces the actual three-candidate scorer with84 exact
checks (maximum residual0). Adding the two full-transferred-ion geometries never
changes either selected endpoint, under either mathematical or operational rules.
Both sources retain adaptive_Ca for Ca and adaptive_La for La. Thus the apparent
1.141070kcal/mol target gain against the earlier incomplete two-candidate pool
is not an improvement over the real scorer. No class/calibration/default changes.

Target Ca's actual adaptive_Ca work is−9.698819kcal/mol, already lower than the
transferred-origin candidate's−9.039922. The large reciprocal donor-response
change remains real within the tested model, but is separate from scoring utility.
Close this endpoint-transfer proposal rather than calling it an improvement.

## New physical information from saved forces

Force analysis4134 admits20/20 saved metal vectors and cap-aware donor projections.
Coordinate and own-metal saved-gradient replays pass. These are MACE vacuum
proposal forces, not gradients of the composite solvent descriptor.

| Host, saved adaptive_La arrangement | Metal | F along transfer at original | F at full transfer |
|---|---|---:|---:|
| Target3 | Ca | +16.324543 | −10.448756 |
| Target3 | La | +32.503637 | −35.991480 |
| Control0 | Ca | +5.410100 | −29.668637 |
| Control0 | La | +21.994408 | −63.549962 |

Units kcal/mol/Å; positive force points toward the transferred position. Under
smooth local energy this brackets an interior stationary point along that line,
not necessarily a unique minimum and not a solvent-composite stationary point.
On the actual selected adaptive_Ca geometry, original Ca projections are+20.461715
for target and+2.003816 for control. Full endpoint positions can miss useful
interior positions, even when the endpoint pool has no benefit.

The resulting next experiment is one fixed halfway point on the same vectors,
with both saved adaptive_Ca/adaptive_La donor arrangements and both metals in
both hosts. Eight new cells, unchanged gates, no optimizer or automatic finer
spacing. See ../exaa_metal_midpoint_20261014/PLAN.md, declared before construction
or scoring. This is a new development probe, not rescue of a confirmatory test.
Even lower selected energies will not alone demonstrate classifier improvement.

## Actual cost and records

These two analyses made zero molecular calls.4133:3scheduler seconds/1CPU,
TotalCPU1.916s,MaxRSS20004K.4134:1second/1CPU,TotalCPU0.818s,MaxRSS49372K.
Full unrounded outputs are in the respective full_pool_check_v1/execution/
FINAL_COLLECTION.json and force_review_v1/execution/RESULT.json under the
matching workspace. FULL_POOL_FORCE_RESULT.json indexes exact hashes/results.
All stages in this namespace are complete; delayed callbacks require ACKs only.
