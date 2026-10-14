# Independent review: conditional 1FLG crystal result

**Chain A produces a complete, numerically admitted shared pool and a conditional Ca call under the existing accommodated calibration. This does not demonstrate that accommodation physically favors Ca.** Chain B remains unavailable. No parameter, threshold, source selection or result was changed in this review.

## Energy response and classification

Both metals select the La-generated geometry under both operational and mathematical policies. The saved works relative to their own origins are:

| Metal | MACE vacuum | Frozen-charge solvent | Total selected work |
|---|---:|---:|---:|
| Ca | −6.061703 | +0.891899 | −5.169804 |
| La | −8.069710 | +0.969909 | −7.099801 |

Units are kcal/mol. The reported Ca-minus-La contrast therefore increases by **+1.929997 kcal/mol**. La gains more energy from the admitted motion. This is a shift toward La on the raw model axis, not toward Ca. Ca selecting the other metal's proposal also confirms that the actual common pool matters here.

| Policy | Reported contrast | Frozen Ca maximum | Frozen La minimum | Conditional call |
|---|---:|---:|---:|---|
| Origin | −405468.157453 | −405468.826318 | −405457.367232 | Inconclusive |
| Operational / mathematical | −405466.227456 | −405461.963102 | −405453.512779 | Ca |

These are different, previously frozen canonical calibrations. The transition from inconclusive to Ca reflects the score together with that policy's calibration. It cannot be presented as a Ca-directed accommodation mechanism. Neither the large raw totals nor their differences are binding affinities.

## Actual acceptance evidence

All six metal/geometry solvent cells are available. Their stored xTB/tblite component parity passes; the largest listed component error magnitude is 2.378312e−8 hartree and population discrepancy is 7.700865e−7 electrons. All four native water-diagonal component checks per cell pass, and every vacuum-charge precision bound is informative under the frozen policy. All 24 stored backend endpoints are available. The older absolute-energy parity flag remains false in all 12 point records; it has not been relabeled as passing. Its retained absolute offsets do not replace the separately qualified solvent-component gate.

The two actual MACE fixture checks pass with exact stored energies. Both optimizer results report success, physical feasibility and no angular boundary; Ca used 15 iterations and La 14. Fixed atoms remain fixed and no new severe overlaps are reported. These are bounded vacuum-Hamiltonian proposals selected by endpoint energies, not demonstrated stationary minima of the composite energy. The solvent stage reused the six MACE values without new MACE calls.

## What the source comparison establishes

The saved comparison retains all ten modeled Q9/ExaA sources plus both crystal chains. It separately displays the five historical one-Ca preparation failures rather than counting them as additional current observations. Existing inconclusives remain: La-conditioned samples 2, 3 and 4, and two-Ca sample 3. The admitted crystal chain provides another conditional Ca-compatible source; it does not resolve all modeled-source sensitivity or prove that their folds are wrong. All concern one already-consumed protein, not independent biological validation.

The crystal preparation is **dry, with no explicit waters, and assumes oxidized PQQ3−**. Its redox state is not experimentally qualified. The admitted 192-atom local state matches the Q9 state signature; Ca/La charges are −3/−2 and both multiplicities are one. The mapped mature sequence is exact over 582 residues, but the full constructs differ. The separately admitted remote terminal OXT completion remains documented. Chain B retains the unchanged exclusion: mapped catalytic Asp partner at 3.511 Å exceeds the 3.5 Å gate. No score was supplied for B.

## Measured scope and records

Native job4095 records 33 actual model forwards and 427.601902 seconds through export (2,993.213313 allocated CPU-seconds). Solvent job4096 records 128.447023 seconds (899.129162 allocated CPU-seconds), with 12 vacuum/water inputs evaluated by both backends and 12 native ALPB component evaluations. These stage timings exclude earlier preparation, queueing and historical work; they are not a GPU-throughput estimate. Saved-output comparison4106 reports 0.080950 seconds and zero molecular calls.

`RESULT_REVIEW_PINS.json` pins the actual preparation, native collection, fixture checks, solvent/component collections and comparison, plus all six native replay results and their twelve point records. Principal result hashes:

- Solvent collection: `c3711fcb0be689a977082f29bcc9af52ccfd69f59d639b5c158a05f60283854b`.
- Source comparison: `fd5c7b663814be6783ca2e5cd160db38e3b115d3c66aaf52ebc20ab789dae351`.
- Frozen canonical reference: `ee00fa03aba64c23e6a737fb2b3b4a9f40dc139e58bdb9225e8ed10c9af2b58d`.

This review read completed artifacts and verified metadata hashes. It performed no new molecular evaluations, scientific recomputation, refitting or production promotion.
