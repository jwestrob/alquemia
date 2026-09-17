# MACE-OMOL: fast PQQ transfer succeeds; broader robustness fails

**The new vacuum descriptor separates all25 canonical calibration cases and
correctly classifies all three consumed transfer cases with its own frozen
calibration rule. It does not resolve the broader affinity failures.** Baseline
scoring, references and defaults remain unchanged.

## Predictive results, kept separate

The25 calibration scores yield a79.03060796106001kcal/mol gap. Their own raw-R
bands are Ca-supported <=-405391.7169957054 and La-supported >=-405312.6863877443.
The open interval is inconclusive. No transfer case set those bands, and no
baseline aquo value, old threshold, solvent energy or fitted weight was added.
A larger gap on this different energy model is not evidence of better accuracy.

| Transfer structure | Expected class | Raw R, kcal/mol | Decision |
|---|---|---:|---|
| 1H4I | Ca | -405405.6231521251 | Ca-supported |
| 4MAE | La | -405303.4009723768 | La-supported |
| 1KB0 | Ca | -405395.5839287484 | Ca-supported |

The separate non-PQQ test requires positive alpha-minus-GGR contrasts. It gives:

| Comparison | OMOL delta R, kcal/mol | Matched DFT baseline delta R | Result |
|---|---:|---:|---|
| 1F6S minus GGR extended58, primary | -18.471538477 | -14.918692215 | Fail |
| 6IP9 minus GGR extended58, primary | -19.428069552 | -17.309144509 | Fail |
| 1F6S minus GGR connected111, robustness | +7.720422072 | -7.575191341 | Pass |
| 6IP9 minus GGR connected111, robustness | +6.763890998 | -9.965643635 | Pass |

The two GGR representations differ by **-26.191960549680516kcal/mol** in OMOL,
versus -7.34350087327 for these DFT references. Selecting the larger core after
seeing the direction would conceal the failed primary test. Both are retained;
primary and representation-robustness acceptance both fail. This is evidence of
sensitivity to preparation, not a unique diagnosis of its cause.

All cases were already consumed. The two alpha structures are one biological
group, and the two GGR cores are one structure/group. Alpha evidence retains
its condition/construct qualifications. Both1H4I/P16027 and4MAE/I0JWN7 overlap
calibration accessions. PQQ functional class differs from La/Ca affinity; motif
and charge already separate that calibration. No absent model-training overlap,
independent28-protein accuracy, broad affinity or incremental-physics claim.

## What changed and what actually ran

New protocol: `mace_omol_0_100m_vacuum_descriptor_v1`. Installed MACE0.3.16 and
torch2.8.0 were reused without upgrades. The official100M checkpoint was
SHA256-verified and inspected:83elements including Ca/La,52,365,482parameters,
float64,6Angstrom cutoff, `omol` head and categorical total-charge/spin inputs.
The input batch actually receives each recorded formal charge and multiplicity1.
No POLAR field/charge adapter is applied. The model does not produce an atomic
charge distribution in this path, so missing density is explicit, not a failed
charge sum silently accepted. Parameters remain unchanged during inference.

The [maintainer's release](https://github.com/ACEsuit/mace-foundations/releases/tag/mace_omol_0)
and [OMol dataset definition](https://fair-chem.github.io/omol25/) describe the
pretrained molecular model and its vacuum wB97M-V/def2-TZVPD target. This is a
separate descriptor, not a native r2SCAN-3c replacement with an inherited gauge.
The checkpoint's ASL license is recorded; the work is an academic pilot, with
no commercial distribution or production promotion.

Job1200797 completed10numerical-qualification calls; all9checks pass. Largest
rotation energy error1.24e-6kcal/mol and force error5.71e-8eV/Angstrom are below
frozen0.01/0.001 tolerances. Repeat/translation also pass. Job1200799 completed
60newbenchmark calls plus4exact qualified endpoint reuses. **70newOMOL calls,
zeroDFT/solvent calls, no failed scientific attempt.** Native source coordinates,
caps, protonation, PQQ state, explicit waters and source mappings are unchanged.

Existing runner dispatch now supports typed OMOL receipts, finite qualification
and benchmark preparation, exact caches, partial-job recovery and read-only
reports. The score is E_Ca-E_La, eV converted once. Solvent, S, aquo reference,
relaxation and entropy remain unavailable. Numerical feasibility and prediction
are reported independently. Every original baseline comparator remains visible.

## Measured affordability and verification

Total recorded new allocation: **670GPU-seconds,10,720allocatedcore-seconds,
753.691reportedactualCPU-seconds**. OneA5000,16CPUs,64,474MiB host allocation per
job. Median model inference is1.454023641seconds per La/Ca pair across32core
preparations (range1.3320–2.0395s). Whole jobs take longer because they include
fresh worker imports, model loading, validation and collection. Peak GPU tensor
allocation3,170,962,432bytes, reserved4,162,846,720bytes; worker host RSS peaks at
1,543,488KiB. No matched-hardware DFT speedup or persistent-process throughput
is claimed. [Actual costs](costs.json) retain allocations and per-case timings.

Checkpoint download took7.013s; metadata inspection14.565s. Preparation and local
checks were not fully profiled and are unavailable, not zero. The first report
read preceded visibility of the completed collection; the same finished output
was then read successfully without rerunning any scientific calculation.

Six real-fixture OMOL tests cover state/charge/checkpoint corruption, exact
inventories/reuse, real energy units/missing results and separation of successful
PQQ calibration from failed primary non-PQQ evidence. Existing baseline-runner
regression:10passed with2explicit environment-dependent skips. Actual70GPU
calculations are distinct from these parser/runner checks; none is fabricated.

## Judgment and next use

- Numerical consistency at the declared precision: passes the real core tests.
- Predictive usefulness: demonstrated narrowly on the retrospective canonical
  PQQ class transfer; broader affinity robustness fails.
- Affordable: measured small-core inference is seconds and fits comfortably on
  one standard GPU. End-to-end production throughput is not yet optimized.

**Recommendation: pursue OMOL further as an opt-in candidate while retaining the
baseline.** The immediate scientific issue is the representation-sensitive
GGR/alpha comparison, not another threshold adjustment. This result is useful
progress toward the MACE goal; it does not complete broad discrimination or
justify automatic deployment. Full endpoints, source pins, all cases and
separate gates are in [result.json](result.json).

[COMMANDS.md](COMMANDS.md) provides runnable preparation, dry-run, finite execute,
collection and report commands. Both OMOL jobs are complete. A next read-only
operation is the report replay listed there, which launches no inference.
