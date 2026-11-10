# The gamble: protein fingerprints in natural rare-earth patterns — report (2026-11-10)

Lead agent: Claude Fable. Mac, CPU, public data. Plan frozen in `PLAN.md` before any natural data was read;
sample rule for the plume frozen in `ADDENDUM_test1_rule.md` before deep-water values were examined.

## Verdict: dead, under the rule frozen in the plan

The two marks that separate protein sites from small chelators in pure solution fail as a geochemical pointer
to lanthanide-binding proteins. Three measured reasons:

1. **Biology in nature acts at the other end of the series.** The two natural biological patterns available
   (methane plume, symbiont-bearing mussel gills) involve La, Ce, Pr and Nd and stop there. The marks sit at
   Gd, Tb, Ho and Y, where these pools hold no biological signal.
2. **Iron oxide, the main abiotic scavenger, rejects yttrium twice as strongly as lanmodulin does** (Y − Ho
   −0.41 against −0.21 log10). The plan assumed abiotic scavengers would look like small chelators.
3. **Carbonate in the water shifts every reading.** In carbonate solution iron oxide lands inside the range the
   plan predicted for proteins (4 of 5 measured runs), and a lanmodulin in seawater lands outside it (6 of 6
   computed settings).

One by-product is worth keeping: the natural biological sink prefers the largest ion (La over Nd about 10-fold
in the water, about 4-fold after correcting for free-ion availability), a direction of selectivity that none of
the 353 lanmodulins and one of the 72 chelators in the project tables shows.

Figure: `figures/biosignature_gamble.png`.

## Step 0: the marks are readable on single patterns in pure solution

- Fingerprint of one log10 pattern: cubic polynomial in nine-coordinate radius plus an additive step for Tb–Lu;
  step in element spacings u; yttrium offset. Variant without Ce and La.
- Chelators against lanmodulins, one pattern at a time (AUC): step 0.81; |step| 0.88; yttrium offset 0.84
  (oxygen-donor chelators 0.98). Stop rule passed.
- The plan's "protein-type" range (|u| below 1 and Y 0.1–0.3 log10 below Ho) holds 339 of 353 lanmodulins and
  4 of 50 oxygen-donor chelators.

## Test 1: Deepwater Horizon methane plume — "uptake-like"

Data: Shiller's per-sample dissolved Y and lanthanides (GRIIDC, 145 samples) joined to oxygen for the same
bottles. May 2010, 800 m or deeper: 43 samples.

**R1, the article's Figure 2 rebuilt.** Five samples with the lowest La/Yb (1,060–1,170 m; oxygen deficit +7 to
+37 µM) against the other 38.

| Element | Removal (95% interval) | Rate relative to La, selective part |
|---|---|---|
| La | 54.7% [40.9, 68.0] | 1 |
| Ce | 18.5% [−2.9, 35.8] | — |
| Pr | 20.0% [9.6, 30.1] | 0.26 [0.21, 0.31] |
| Nd | 8.3% [0.3, 16.4] | 0.09 [0.04, 0.14] |
| Sm | 1.8% [−3.5, 7.1] | 0.00 [−0.02, 0.02] |
| Gd | 1.8% [−3.6, 6.9] | 0.00 [−0.02, 0.02] |
| Er–Lu (mean) | 1.9% [−2.7, 6.4] | 0 (defines the uniform part) |
| Y | 0.5% [−6.0, 6.9] | −0.02 [−0.05, 0.01] |

- "Selective part" = removal beyond the uniform 1.9% seen in Er–Lu.
- The same five samples are selected by lowest La and with a 1,000 m depth floor.
- Sample by sample, La removal runs from 37% to 74% and the ratios hold: Pr 0.16–0.35, Nd 0.02–0.18.

**R3, selector independent of the rare earths.** Slope of ln(concentration) on oxygen deficit, 39 samples at
900–1,400 m: La 16.1% per 10 µM [5.7, 23.7], Pr 5.7%, Nd 3.1%, Sm 1.5%, heavy elements 1.8% (uniform). Selective
part relative to La: Pr 0.27, Nd 0.09, Sm −0.02. R1 and R3 agree.

**Frozen reading:** Nd removal ≥ La removal: false. Sm–Gd removal detectable: false in R1 and R3. Removal falls
La > Pr > Nd: true. Reading: **uptake-like**.

**Against the binder predictions (rate relative to La):**

| | Pr | Nd | Sm | Gd | Yb |
|---|---|---|---|---|---|
| Plume, observed | 0.26 | 0.09 | 0.00 | 0.00 | 0.00 |
| Lanmodulin capture, free ions (the plan's prediction) | 1.86 | 1.92 | 1.97 | 1.16 | 0.24 |
| Lanmodulin capture × free-ion fraction | 0.74 | 0.68 | 0.41 | 0.28 | 0.01 |
| Free-ion availability alone | 0.40 | 0.35 | 0.21 | 0.24 | 0.06 |
| Iron oxide, carbonate solution (measured, median of 5 runs) | 0.94 | 0.56 | 0.54 | 0.35 | 0.07 |
| Iron oxide at seawater carbonate level (computed) | 1.28 | 1.36 | 1.24 | 0.88 | 0.40 |

- Free-ion fractions computed from published carbonate, sulfate, chloride and hydroxide constants. Relative to
  La across the carbonate range and both published carbonate sets: Pr 0.27–0.58, Nd 0.23–0.54, Sm 0.12–0.38.
- The plume pattern is steeper than every equilibrium binder, with or without the speciation correction.
- **The sink's own preference** (observed ÷ free-ion ratio): La 1, Pr 0.66 [0.36, 1.16], Nd 0.25 [0.08, 0.62],
  Sm below 0.15. Roughly half of the apparent La-over-Nd preference (in log units) is free-ion availability and
  half belongs to the organisms.
- Pr over Nd: 3.0-fold [2.2, 5.3] in R1 and 3.1-fold in R3; free-ion availability supplies 1.13-fold of it.

## Test 2: mussels with methane-eating symbionts — the same uptake pattern, no protein-type reading

Data: Wang et al. 2020, supplementary tables. 14 methanotrophic mussels (two South China Sea seeps) with gill
and mantle; 3 thiotrophic controls from the same sites.

- Gill La: median 1,822 ng/g, 17 times the mantle. Thiotrophic controls: 0.6 times the mantle.
- Gill excess per unit of dissolved element, relative to La (median, quartiles): Pr 0.110 [0.089, 0.135]; Nd
  0.028 [0.024, 0.037]; Sm 0.006 [0.002, 0.017]; Tb 0.002; Dy 0.001. Pr over Nd: 3.9-fold.
- Fingerprint (the 4 gills with Gd reported). Gill over deep water: u +0.30 to +2.07, Y − Ho −0.03 to −0.12; 1
  of 4 in the protein-type range. Gill over mantle: 0 of 4. Thiotrophic controls: 0 of 3 in either
  normalization.
- Y/Ho mass ratio: gill 40, mantle 49, local sediment 30, deep water 52. Gills carry particle-derived material,
  so their heavy-element pattern mixes sources.
- Deviation from the plan: local seawater data were unavailable. Normalizers are the same animal's mantle and
  the deep Gulf of Mexico background from test 1.

## Test 3: abiotic reference — iron oxide

Data: Kawabe et al. 1999 (two papers), distribution coefficients between freshly precipitated iron oxyhydroxide
and solution for every lanthanide and Y. 16 runs in salt solution (pH 5.3–6.4, free ions), 5 runs with
bicarbonate (pH 8.2–8.6).

| | Step u | Y − Ho (log10) | In the protein-type range |
|---|---|---|---|
| Lanmodulins, pure solution (353) | −0.05 [−0.18, +0.06] | −0.21 [−0.23, −0.18] | 339 of 353 |
| Oxygen-donor chelators, pure solution (50) | +0.92 [+0.28, +1.63] | −0.08 [−0.10, −0.03] | 4 of 50 |
| Iron oxide, salt solution (16 runs) | −1.76 [−2.82, −0.81] | −0.41 [−0.54, −0.37] | 0 of 16 |
| Iron oxide, carbonate solution (5 runs) | +0.18 [−0.02, +0.88] | −0.28 [−0.35, −0.24] | 4 of 5 |
| Iron oxide in seawater (computed, 6 settings) | +0.14 [−0.33, +0.42] | −0.26 [−0.31, −0.11] | 5 of 6 |
| Lanmodulin in seawater (computed, 6 settings) | +0.28 [+0.10, +0.41] | −0.04 [−0.09, +0.11] | 0 of 6 |

(Lanmodulins and chelators: median and quartiles. Other rows: median and full range.)

- The plan predicted chelator-type values for oxides. Observed: neither type. The oxide's step is an additive
  +0.05 log10 for Tb–Lu (positive in 16 of 16 runs) on a locally falling curve, the classical tetrad form.
- Carbonate complexation itself prefers Ho to Y: the free fraction of Y exceeds that of Ho by +0.11 to +0.31
  log10. Every binder's Y − Ho reading rises by that amount in carbonate-bearing water.
- The difference between two binders in the same water is fixed: iron oxide sits 0.22 log10 below lanmodulin in
  Y − Ho in every setting. A reading therefore needs a co-located abiotic reference or the water's carbonate
  chemistry to about 0.05 log10; the two published carbonate sets differ by about 0.1 on this term.
- In carbonate-free water the two numbers separate iron oxide from lanmodulins completely (AUC 1.00 for each;
  16 runs against 353 proteins).

## Test 4: seawater (descriptive)

Deep Gulf water has a Y/Ho mass ratio of 52; marine sediment (South China Sea seeps) 30. The 0.25 log10 excess of
dissolved yttrium
matches the sign and size expected from iron-oxide scavenging in seawater (−0.26). Abiotic chemistry accounts
for it.

## The by-product: nature's sink prefers the largest ion

| | La over Nd | Pr over Nd |
|---|---|---|
| Methane plume (removal from water) | 11-fold | 3.0-fold [2.2, 5.3] |
| Mussel gills (accumulation per unit dissolved) | 36-fold | 3.9-fold |
| Free-ion availability alone | 2.9-fold (1.9–4.3) | 1.13-fold |
| Lanmodulin mean (353) | 0.52-fold (Nd preferred) | 0.97-fold (Nd preferred) |
| Most La- or Pr-leaning lanmodulin of 353 | 0.94-fold | 1.11-fold |

- La binds more strongly than Nd in 0 of 353 lanmodulins and in 1 of 72 chelators (malonate, by 0.04 log10); La
  more strongly than Sm in 0 of 74 chelators.
- The proteins responsible for the natural preference are unidentified here. Published candidates are the uptake
  pathway and the dehydrogenase itself (UNVERIFIED for these organisms).

## Limits

- Free-ion fractions use 25 °C infinite-dilution constants, Davies activity coefficients at ionic strength 0.7
  and an assumed carbonate activity range (log a −5.8 to −5.1). Plume water was at 4.7 °C and 1,100 m.
- The lanmodulin profile was measured at pH 3 (SpyCI-LAMBS library; the 16 reserved proteins stay sealed).
- Gill selectivity uses Gulf of Mexico deep water as the stand-in for South China Sea water; Gd is missing for
  10 of 14 gills.
- One laboratory supplies the iron-oxide data. Manganese oxide, clays and organic surfaces were left
  unexamined (tables behind paywalls).
- The novelty of the by-product against the geochemistry literature is unchecked.

## Ideas this leaves (unverified, nothing started)

1. **A selectivity class missing from the binder tables.** Largest-ion preference at 3-fold per element lives
   somewhere in the uptake chain of marine methanotrophs. Plume-bloom and symbiont genomes are the place to look.
2. **A measured template for the ocean idea** (idea 2 of the earlier list). Biology removes La with Sm
   untouched (rate ratio 0.00); iron oxide removes Sm as fast as La (0.7–2.2 computed for seawater). The pair
   (La deficit, Sm deficit) separates the two processes in any dissolved profile.
3. **Where the original marks could still be read:** two pools in the same water (the carbonate shift cancels),
   or acidic carbonate-free waters, where iron oxide and protein sites separate completely on both numbers.

## Sources

- Shiller, Chan, Joung, Redmond & Kessler 2017, Sci. Rep., doi 10.1038/s41598-017-11060-z. Data: GRIIDC doi
  10.7266/N7CC0XQ6 (rare earths) and 10.7266/N7MS3QQ5 (oxygen, nutrients, metals).
- Kawabe, Ohta, Ishii, Tokumura & Miyauchi 1999, Geochem. J. 33, 167–179; Kawabe, Ohta & Miura 1999, Geochem. J.
  33, 181–197.
- Ohta & Kawabe 2000, Geochem. J. 34, 439–454 and 455–473 (carbonate constants; Table 1 of the second paper also
  lists Liu & Byrne 1998, Millero 1992 and Mironov et al. 1982 values).
- Wang, Barrat, Bayon, Chauvaud & Feng 2020, Geochem. Persp. Lett. 14, 26–30, doi 10.7185/geochemlet.2019,
  Supplementary Tables S-2 to S-8.

## Files

- `PLAN.md`, `ADDENDUM_test1_rule.md` (hashes beside them).
- `scripts/`: `fingerprint.py`, `step0_reference.py`, `parse_shiller.py`, `test1_plume.py`, `transcribe_refs.py`,
  `transcribe_mussels.py`, `tests23_references.py`, `make_figure.py`.
- `results/`: `step0.log`, `test1_plume.log`, `tests23.log`, `test1/`, `tests23/`.
- `data/`: `griidc_download/` (two workbooks, sha256 in `SOURCES.md`), `reference_tables/` (transcribed tables),
  `lit_untrusted/` (downloaded papers, treated as data).
