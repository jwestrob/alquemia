# Complete-context PQQ reference test

**The candidate retains 25/25 canonical PQQ classifications and 3/3 consumed crystal transfers.** Its calibration gap contracts from 79.030608 to 2.315459 model-kcal. Combined with the separately completed alpha/GGR direction gain (2/6→6/6), this is useful limited evidence, with a substantially weaker PQQ margin. The research scorer is not promoted; production/default remain unchanged.

## Frozen model and decisions

`native_OMOL_complete_polar_context_disulfide_v2` uses the existing native 100M OMOL checkpoint, float64 vacuum, actual endpoint charge/singlet, exact source coordinates and original PQQ/water/protonation. The original 3.3 Å donor-anchor / 3.5 Å complete polar-neighbor rule is unchanged. The uniform source-disulfide extension supports the two previously unsupported Cys contexts; all other 23 preparations remain exactly identical. No source heavy atom or retained original cap was moved. Complete overlaps replace original cut-bond caps through source connectivity.

The raw descriptor is `R = (E_Ca − E_La) × eV_to_model_kcal`, converted exactly once. It is a vacuum learned-energy contrast, not the CPCM DFT scale, affinity, or a complete binding free energy. No old aquo reference, CPCM contribution, environmental scalar, fitted combination or baseline threshold was added.

The predeclared class-extrema procedure uses only the 25 calibration cases. Frozen bands are:

- Ca-supported: `R ≤ −405347.73205159395` model-kcal.
- La-supported: `R ≥ −405345.416592485` model-kcal.
- The open interval between these values is inconclusive.

All 25 are scored, with zero unsupported or failed cases. The minimum required calibration gap was 0.02 model-kcal before execution. Transfers do not set the bands. All cases were already consumed; 1H4I/P16027 and 4MAE/I0JWN7 share sequence accessions with calibration structures. This is structural transfer, not new independent blind validation. Composition and motif already separate the canonical panel.

| Transfer | Expected region | Context R (model-kcal) | Decision |
|---|---|---:|---|
| 1H4I | Ca | -405428.473351463 | Ca-supported |
| 4MAE | La | -405338.648146804 | La-supported |
| 1KB0 | Ca | -405410.248491153 | Ca-supported |

## Why the margin contracts: existing-output evidence

The upper Ca edge is Q9Z4J7; the lower La edge is A0A3F2YLY8. All four contexts gaining formal charge −1 shift R upward. Every charge-neutral expansion shifts it downward. These are descriptive strata, not label-conditioned scoring rules; no case is excluded or corrected.

| Calibration stratum | Structures | Mean context−core R | Range |
|---|---:|---:|---:|
| Ca, added charge -1 | 3 | 49.848406 | 49.128656 to 51.025181 |
| Ca, added charge 0 | 11 | -13.470841 | -19.475292 to -1.020130 |
| La, added charge -1 | 1 | 19.776419 | 19.776419 to 19.776419 |
| La, added charge 0 | 10 | -30.238487 | -38.366649 to -26.571630 |

The Ca class width grows 15.819355→73.516161 model-kcal; La grows 16.629285→69.135909. This explains why perfect calibration fidelity alone is insufficient. The three charged Ca contexts now have the same total charge as neutral-expanded La contexts. Geometry of the original core did not change; adding real anionic chemistry can legitimately favor La. Native OMOL also broadcasts its categorical total-charge embedding to every atom. The current endpoints have no saved atomic-energy or charge/density decomposition, so the association cannot identify how much is physical electrostatics versus charge-conditioning/representation effects. No empirical subtraction or total-charge change is justified by these data.

The follow-on physical check is recorded separately: unchanged native CPCM DFT on all four charge −1 contexts, with their actual archived original-core endpoints and existing neutral-context controls. It cannot restore blind status or serve to tune a new threshold.

## Actual execution and cost

- Job1202474 completed all 52 new native OMOL endpoints. Four exact context endpoints (1H4I/4MAE) and all original-core scores were reused only after model/input/state checks. Zero new DFT, retries, training or folding.
- 70 allocated GPU-seconds; 16 CPUs ×70s = **1,120 allocated CPU-seconds**. Worker67.546548s; summed endpoint evaluation11.594540s, median0.170356s. Repeated model loading dominates this small job.
- Peak CUDA allocation9,092,231,168bytes; peak worker host RSS2,666,152KiB; scheduler batch MaxRSS2,004,600KiB. Metrics measure different processes/scopes.
- Local preparation, testing and reporting are additional. No controlled matched-hardware end-to-end DFT speedup was measured. Cached folding/protonation are excluded from this inference timing.

Five real-source preflight tests passed (zero skips,16.860s). They cover actual disulfide topology, exact preservation of the other23 contexts, charge/paired-coordinate/reuse invariants, and explicitly corrupted real-source rejection. Scientific calculations above are distinct from those tests. Source-scoped actual-output regression is recorded in TESTS.txt.

Exact outputs and all 28 row values are in RESULT.json, which pins the original result/manifest/receipts. COSTS.json and COMMANDS.md record cost and runnable operations. The original immutable job artifacts remain under workspaces/environment_pqq_20260919/prepared_v1/.

**Recommendation:** retain this as a promising compact research candidate and test its charge/composition robustness against native DFT. Keep the established baseline available; do not promote this context scorer from canonical fidelity alone.
