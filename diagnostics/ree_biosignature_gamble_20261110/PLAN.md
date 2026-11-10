# The gamble: protein fingerprints in natural rare-earth patterns (frozen 2026-11-10, before any natural data was read)

Jacob's charge of 10 November ("1 sounds like a crazy ass gamble, why don't you give it a shot… No huge compute
if possible. No h200."). Lead agent: Claude Fable. Mac, CPU, public data only.

## The idea

Two analyses of 10 November give fingerprints that separate binding by protein sites from binding by small
ligands:

| Fingerprint | Protein sites (353 lanmodulins; public protein series) | Small chelators (71 NIST ligands) |
|---|---|---|
| Step between Gd and Tb | none (+10 mA against a smooth curve; no break at that position) | present (+31 mA; heavy branch extra stable) |
| Yttrium | binds like Tm or Er: 0.2–0.3 log10 below a lanthanide of its radius | at its size twin with oxygen donors; lighter with nitrogen donors |
| Shape (lanmodulins) | capture peaks at Nd–Sm; La below Ce–Sm; heavy elements fall away | monotonic for most |

If natural rare-earth patterns carry these marks where biology has moved the elements, geochemistry could show
where to look for lanthanide-binding proteins.

## Step 0: make the fingerprints measurable on a single pattern (data in hand)

For any log10 pattern over the lanthanides (Ce always excluded, because of its redox chemistry; variants without
La and without Eu reported): fit a cubic polynomial in Shannon nine-coordinate radius plus one additive step for
Tb–Lu.

- **Step h** (log10 units): the step coefficient. **Step in element units u** = h divided by the pattern's
  local slope across one Gd–Tb spacing (12 mA); reported only where that slope is at least 0.02 log10 per step.
- **Yttrium offset y** (log10): yttrium minus the fitted curve at radius 1.075 A on the heavy branch.
- Reference distributions: 353 lanmodulins, the public protein and peptide series, 71 chelators. Report how
  well (u, y) separates proteins from chelators one pattern at a time (area under the ROC curve), because a
  natural sample is a single pattern.
- **Stop rule for step 0:** if single-pattern separation is poor (AUC below 0.8 for u and for y), the marks are
  too faint to read in natural data and the gamble ends there.

## Natural tests (each pattern must be a BINDING or REMOVAL pattern, i.e. a ratio between two pools)

1. **Deepwater Horizon methane plume** (Shiller et al. 2017, Sci. Rep.): fraction of each element removed in
   plume water relative to background water.
   - Protein-capture prediction (lanmodulin-like): removal follows the lanmodulin mean profile, so Nd and Sm are
     removed at least as much as La, and Gd is removed too (about 60% as much as Sm).
   - Known-biology alternative (uptake for the dehydrogenase): La ≥ Ce > Pr > Nd and nothing beyond.
   - Reading: "lanmodulin-like" if removal of Nd ≥ removal of La and removal of Sm–Gd is detectable; "uptake-like"
     if removal falls monotonically from La to Nd and vanishes beyond.
2. **Organisms that accumulate rare earths with the help of methanotrophs or other biology** (for example seep
   mussels with lanthanum anomalies): tissue over ambient seawater.
   - Protein-type prediction: |u| below 1 and yttrium 0.1–0.3 log10 below Ho.
   - Chelator-type prediction: u about +2; yttrium within 0.05 log10 of Ho.
3. **Abiotic references from the literature where tables are open** (iron and manganese oxides, bacterial cell
   surfaces, humic substances, marine particles): the same two numbers. Chelator-type is expected for oxides.
4. **Seawater itself** (dissolved patterns): descriptive; the scavenger's fingerprint is inferred only if measured
   carbonate constants for every element can be obtained.

## Verdict rule (stated now)

- **Has legs:** at least one natural biological binding pattern shows protein-type values in both numbers while
  the abiotic references show chelator-type values.
- **Dead:** step 0 fails; or natural biological patterns cannot be told from abiotic ones; or the needed ratios
  cannot be formed from public data.
- Anything in between is reported as "unresolved", with what data would resolve it.

## Outputs

`scripts/`, `data/` (downloaded public tables with source notes), `results/`, `figures/`, `REPORT.md`.
