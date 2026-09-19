# Coupled water response: useful relaxation, occupancy still unsupported

**The DFT-anchored MACE model successfully guides water relaxation on this pilot.**
All11 proposed steps lower actual native DFT energy and pass the frozen energy
checks. Three of four endpoints reach the prescribed water-coordinate minimum
criteria. One calcium endpoint remains nonstationary. The local harmonic model
is not qualified for entropy or occupancy. Baseline/default/PQQ are unchanged.

## What was already available; what changed

We reused the actual12-pattern/24-endpoint occupancy table, native analytic
DFT gradients, contextual MACE water orientations, liquid-water reference, and
successful four-state radial check. No source panel was rerun or relabeled.

New code supplies six physical coordinates per variable water: three COM
translations and three rotations about its COM. It retains all cross-water
curvature, actual masses/inertia, the rotation chain rule and exact fixed-atom
coordinates. The Cartesian anchored potential is

`V(q) = E_DFT(0) + E_MACE(x(q)) - E_MACE(x(0)) + (g_DFT - g_MACE)·(x(q)-x(0))`.

MACE is the native unmasked OMOL0-100M checkpoint; DFT remains native
r2SCAN-3c/CPCM(Water)/DefGrid3/TightSCF with analytic EnGrad. This is an approximate
local response model, not self-consistent MACE/CPCM or a global protein potential.
The rotational contribution of the Cartesian anchor is retained in curvature.
No numerical DFT gradients, full DFT Hessian, new force field, spring fit or MD
was used; the few native displacements check projected response only.

The same existing runners now support preparation, dry-run, execution, collection,
comparison, two bounded recentering rounds and reports. Protocol IDs:

- `native_r2scan3c_cpcm_cartesian_anchored_mace_rigid_water_basin_v1`
- `native_r2scan3c_cpcm_cartesian_anchored_mace_coupled_check_v1`

## Actual results

Stage A completed all20 nonempty metal/water states. All initial and proposed
minimum numerical curvature checks pass, with positive curvature. Four cheap
minima are interior;16 hit the translation boundary. These20 predictions are
not20 native-validated minima.

Native validation used the already consumed1F6S11 and6IP9110 arrangements,
bothCa/La:16 soft/response directional endpoints and4 proposed steps. All four
response directions pass; three of four soft directions pass. The1F6SCa soft
check fails: actual even energy0.012884 versus predicted0.003137kcal/mol,
error−0.009748 against the frozen0.005 tolerance. This failure remains in every
report. Geometry refinement was explicitly separated from entropy qualification
in [the development decision](SEPARATE_RESPONSE_AND_ENTROPY.md).

Two recentering rounds added4 then3 native endpoints, using actual new DFT and
MACE gradients at each center.6IP9La reached the criteria after the first round
and was explicitly reused. No extra calculation was done for that endpoint.

| Endpoint | Total actual energy reduction (kcal/mol) | Final native projected-gradient maximum | Meets minimum criteria |
|---|---:|---:|---|
| 1F6S11Ca | 4.60929 | 0.70245 | no |
| 1F6S11La | 2.63244 | 0.28150 | yes |
| 6IP9110Ca | 2.32673 | 0.17365 | yes |
| 6IP9110La | 1.86610 | 0.26241 | yes |

Gradient components use kcal/mol/Å and kcal/mol/radian, with the declared
1Å/radian numerical scale; the frozen maximum threshold is0.4. These minima
are conditional on fixed protein/metal/outer waters and rigid internal H2O
geometry. The1F6SCa residual barely decreases in the last round, so another
unchanged iteration is not justified as a route to qualified thermochemistry.
Its soft-mode failure and stalled residual are consistent with inadequate local
curvature, but do not uniquely identify the physical source of the error.

All11 actual trial energies decrease. Maximum proposal-energy error is
0.070295kcal/mol; recenter errors are0.000977–0.022960. The first four steps
already capture most of the observed reduction. The complete unrounded values,
gradients and source receipts are in [tables](export_final_v1/TABLES.md) and
[the numerical record](export_final_v1/result.json).

All fixed atoms remain exact and all water shapes/inventories remain unchanged.
Final metal–waterO distances span2.410–2.702Å; the nearest waterO–other-heavy-atom
contact is2.582Å. Maximum total atomic movement is0.421Å after successive bounded
steps. These checks support continuity of the intended water arrangements;
they do not demonstrate a complete basin search.

## What this contributes to discrimination

`R = E_Ca - E_La`; the measured rigid-water response contributes
`DeltaE_relax,Ca - DeltaE_relax,La`. It changes the contrast by−1.97685kcal/mol
for1F6S11 and−0.46063 for6IP9110, relative to their frozen starting geometries.
Thus these motions favorCa relative toLa in both sampled arrangements. The aquo
constant cancels from these changes. No old decision bands are applied to them.

This establishes an actual geometry-response component and checks its cheap
prediction against DFT. It is not new biological accuracy validation. The cases
are one consumed biological comparison;6IP9110 is a selected two-water arrangement,
not the full-water production preparation. Do not translate these numbers into a
new alpha/GGR success rate or into a PQQ classification change. The earlier
water-H preparation improvement remains a separate, preserved result.

## Why occupancy is still unavailable

The largest per-water translation RMS within each endpoint is0.241–0.350Å in
the conditional harmonic model at the latest geometries, beyond the0.20Å local
domain; some rotation RMS values also
exceed0.35rad. These are diagnostics of the cheap harmonic assumption, not
measured physical fluctuations. The existing local checks therefore cannot
justify an unrestricted harmonic basin integral. See [the declared test](THERMAL_EXTENT.md).
We did not clamp modes, enlarge the validated domain, refit a threshold or insert
an entropy correction to obtain probabilities.

Bound internal-water relaxation/vibrations, non-electrostatic solvent terms and
competing orientation basins also remain unresolved. No missing contribution is
zero. [The energy ledger](THERMODYNAMIC_ACCOUNTING_v2.md) subtracts the water
reservoir exactly once and distinguishes rigid response from internal relaxation.
No harmonic free energy, occupancy probability or new calibrated score is released.

## Execution and cost

All nine new allocations completed:27 native analytic DFT endpoints,
27 coupled-response tasks,27 matched MACE endpoint calls;2828 recorded MACE
objective evaluations in total. Twenty native source-center gradients were reused.
No scientific calculation was retried. Four initial collection records failed
strict regenerated-coordinate equality and were recovered from existing outputs;
[details](COLLECTION_RECOVERY.md) preserve the original failure record.

| Work | Jobs | Allocation wall seconds | Allocated core-seconds | GPU-seconds |
|---|---|---:|---:|---:|
| Initial coupled MACE | 1202050 | 408 | 6528 | 408 |
| Native directional/proposal checks | 1202052 | 1861 | 119104 | 0 |
| Matched MACE checks | 1202053 | 154 | 2464 | 154 |
| First recenter MACE | 1202056 | 118 | 1888 | 118 |
| First recenter native checks | 1202058 | 389 | 24896 | 0 |
| First recenter matched MACE | 1202059 | 34 | 544 | 34 |
| Second recenter MACE | 1202063 | 100 | 1600 | 100 |
| Second recenter native checks | 1202064 | 329 | 21056 | 0 |
| Second recenter matched MACE | 1202065 | 26 | 416 | 26 |
| **Total** | | | **178496** | **840** |

This is49.582 allocated core-hours and14 GPU-minutes of development validation.
Peak reported MACE CUDA allocation is2.03GiB; peak worker host RSS is2.05GiB.
All GPU jobs used oneH200/16CPUs/200000MiB. Native jobs used node-64-768g-31,
16MPI ranks per endpoint, four concurrent tasks initially and three in the last
round. Slurm allocated64CPUs even for the final48CPU request; all64 are charged.
[Cost receipts](COSTS.json) retain scheduler and worker measurements, including
unobserved short-job MaxRSS. Login preparation/tests were not resource-metered.
The older live orientation comparator1201825 is outside this cohort and remains
running; it has not been interrupted or counted as a completed new result.

The cheap local response is affordable for this development test. Production
cost/accuracy qualification is not established: repeated native anchoring adds
endpoints, and this is not a matched speed comparison against the baseline.
The cost does not authorize a production rescore.

## Validation and recommendation

71 distinct real-fixture tests pass:15 new coordinate/gradient/curvature/result
checks,33 previous water tests and23 baseline/preparation tests, zero skips.
These parser/algebra checks are distinct from the27 actually executed native
DFT endpoints and the completed native GPU calculations. Frozen analyzers,
input, energy and gradient receipts and exact source mappings remain
under workspaces; compact records are here. No fabricated scientific outputs.

**Recommendation: pursue MACE-assisted water preparation/response; do not release
this local harmonic occupancy model.** The useful next modeling step is a basin
representation that covers the waters' wider configurations, with native checks
and consistent internal-water/solvent accounting. More identical local iterations
or a fitted entropy offset would not resolve the demonstrated limitations.
Baseline/PQQ remain the default and unchanged. No push or production promotion.

Runnable preparation, dry-run, execution, collection, comparison and read-only
report replay are in [COMMANDS.md](COMMANDS.md). The next runnable operation is
the final report replay there; completed compute jobs should not be resubmitted.
