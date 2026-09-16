# A small La/Ca challenge panel

**2026-09-15 — curation complete; future calculations remain proposed.**
Jacob authorized evidence curation while the main agent investigates accuracy.
This work read the existing 68-row release, preparation/provenance records and
primary literature. It generated no coordinates, preparations, energies,
thresholds, fits or new scientific analyses. The original ledger is unchanged.

## Recommendation

Retain **GGR versus bovine alpha-lactalbumin** as the immediately usable,
already-consumed development challenge. Prioritize **Mex-LanM WT/4P2A** for
finishing a charge-preserving mutant comparison. Keep **PqqT WT/K142A/K142D**
as the strongest numeric same-assay family, with its charge and site-assignment
confounds explicit. These serve different purposes and should not be pooled
into one accuracy denominator.

Among the reviewed records, no fresh, fully mapped, opposite-label comparison
with matched gross charge **and** donor composition was verified. More existing
crystal chains would not fill that evidence gap.

| Comparison | Experimental label/readout | Why useful | Existing coordinates and readiness |
|---|---|---|---|
| GGR / bovine alpha | Ca direction / qualitative La direction | Existing cores both have ligand charge −3; tests more than charge alone | GGR 1GLG/2FW0/2FVY; alpha 1F6S/6IP9. All consumed development cases; donor/water inventories differ. |
| Mex-LanM WT / 4P2A | Both La-favoring folding responses; mutation strengthens Ca response | Neutral substitutions preserve nominal direct-donor identities | WT 8FNS Nd crystal and 6MI5 Y NMR; no exact experimental mutant structure verified. First curation priority, not executable now. |
| PqqT WT / K142A / K142D | Same-assay La and Ca ITC; all favor La | Quantitative mutant response plus WT/A near-null comparison | WT 9B1U and Gd-soaked 9B1V; mutant Ca/La geometry and measured-site assignment unresolved. Charge changes are a confound. |
| GGR Q142 neutral variants | La/Ca mutant directions unresolved | Q→N would preserve charge and amide-O donor class | WT sources only; full original Ca/La table and mutant structures required. Do not promote an abstract trend. |
| Aqualysin weak Ca site | Direct Ca direction | Additional independent Ca-favoring chemistry | 4DZT has two Ca sites; assayed weak site remains unmapped. Neither is selected. |
| Engineered CaMLBT site 1 | Reported La preference; conditions qualified | One active engineered site avoids a cooperative multi-site label | New primary source identified; no exact experimental bound structure verified. Distinct from native CaM and LBT3. |

## 1. Existing development anchor: GGR and alpha

GGR's primary Table I gives Kd(Ca) = **25 ± 11 µM** and Kd(La) =
**729 ± 4 µM**, measured by competition against bound Tb at 25°C, pH 6,
10 mM PIPES/100 mM KCl. The pSF3 construct restores three missing terminal
codons; it is not a metal-site mutant. Preparation strips bound sugar, but
zero residual occupancy was not measured. Crystal sugar states differ.
[Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2899690/).

The bovine alpha strong-site La-over-Ca label comes from 1996 NMR competition;
its complete assay conditions and a numeric same-assay ratio remain missing.
The 2019 paper maps La to the same strong cleft in 6IP9 compared with 1F6S,
but its La-only ITC cannot supply a matched Ca/La ratio.
[1996 primary](https://pubmed.ncbi.nlm.nih.gov/8652630/),
[2019 primary](https://pmc.ncbi.nlm.nih.gov/articles/PMC6370903/).

The existing v3/extended preparations establish equal core charge, **not equal
donor composition**. Keep all three GGR sources as one observation and both
alpha sources as one observation. Their ordered score contrast is a development
check; it is neither a blind test nor a common-condition affinity measurement.

## 2. Best charge-preserving candidate: Mex-LanM WT versus 4P2A

The discovery paper reports near-millimolar Ca response for WT and micromolar
response after replacing the four loop-position-2 prolines with alanines,
while strong picomolar Ln response including La remains. This is a change in
**binding-coupled folding response**, not four microscopic site affinities or
a reversal from La to Ca preference.
[Primary discovery](https://pubmed.ncbi.nlm.nih.gov/30351021/).

The sequence change is neutral and leaves the nominal direct oxygen-donor
residues intact. Actual coordination, solvent and conformational populations
are not thereby fixed. No experimental 4P2A structure was verified; the later
primary spectroscopy study explicitly generated mutant simulation coordinates
from the WT NMR structure. Those are not experimental mutant coordinates.
[Primary structural-spectroscopy study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8963139/).

**Next curation:** retrieve the original WT/mutant constants, errors, construct
tags and free-metal buffering details. The SI retrieval attempted here returned
403. Preserve WT and mutant within one LanM family, including Hans-LanM when
grouping homologues. Any new mutant model or scoring protocol needs agreement.

## 3. Best numeric same-assay family: PqqT

| PQQ-loaded construct | Kd(Ca), µM | Kd(La), µM |
|---|---:|---:|
| WT | 64 ± 5 | 6 ± 1 |
| K142A | 60 ± 10 | 7 ± 5 |
| K142D | 150 ± 30 | 0.6 ± 0.2 |

These tag-free constructs were measured by ITC at 298 K, pH 7,
30 mM HEPES/100 mM NaCl, with three independent preparations. K142D supplies
a clear comparative response; WT/A supplies a useful near-null comparison.
[Primary Table 1 and Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC11331073/).

This is **not a charge-preserving test**: K→A removes a positive side chain;
K→D changes its sign under ordinary ionization. WT Lys occupies the analogous
PQQ cation position. 9B1U contains no Ca/La, and 9B1V's Gd sites may be
nonspecific. The measured site cannot be placed by borrowing the canonical
PQQ coordinate. Mutant experimental coordinates remain missing. Treat the
three constructs as one family and keep affinity separate from the paper's
K142D-specific catalytic observations.

## 4. Narrow holds worth resolving

- **GGR Q142N:** the primary abstract supports neutral gateway substitutions
  having similar Ca affinities, but its La-specific numeric table and conditions
  were not recovered. Q→N is a promising same-charge/same-donor-class contrast;
  Q→S changes donor chemistry and Q→D/E changes charge. No mutant affinity
  label is assigned. [Primary](https://doi.org/10.1021/bi952430l).
- **Aqualysin:** the weak site has Ka(Ca) = 6200 and Ka(La) = 1100 M⁻¹ at
  22°C, pH 6, in 20 mM MES with the strong Ca site occupied. A separate
  La-binding site X is not the comparator. The assignment to 4DZT Ca A302
  versus A303 is unresolved; coordination number must not decide it.
  [Primary](https://www.jstage.jst.go.jp/article/bbb/66/6/66_6_1281/_pdf).
- **CaMLBT, newly found primary:** engineered Arabidopsis CaM N-domain,
  site 2 disabled. Reported Kd(Ca) = 33.2 ± 4.5 µM and Table 1 Kd(La) =
  437 ± 259 pM. La fluorescence uses pH 6 MES/KCl with NTA competition;
  Ca-specific SI conditions remain unverified. Prose reports a slightly
  different La estimate; both are retained. No exact bound structure was
  verified. This is an input-gated candidate, not a replacement label for
  7CCO/7CCN or native CaM. [Primary](https://doi.org/10.1002/ejic.202500468).

## Keep other strata separate

Aequorin remains an ordered site vector against protein-level evidence.
Calbindin's proposed opposing site labels remain unverified; parvalbumin is
cross-study support. CD2/CaM graft measurements have assay/pH and structure
gaps. ConA has conditional exclusion evidence; HRP has unresolved La-specific
direction/site mapping. Lanpepsy's multi-ion state cannot provide five
independent positive labels. LanD/LBT3 within-lanthanide preference is not a
La/Ca label. B2/C5 measure competitive retained metal at protein level in a
preprint, with holo-state and stoichiometry gaps. PQQ functional/cofactor
association stays separate from solution affinity. These exclusions are
documented in the compact JSON; prior labels and records were not rewritten.

## Artifacts and stopping point

- `CHALLENGE_PANEL_CURATION.json`: six proposed comparisons, exact inherited
  labels/conditions, priorities, grouping, structure gates, exposure status,
  primary links and SHA256 pins.
- `PRIMARY_WEB_CAPTURES.json` and `CAMLBT_CONSTRUCT_CAPTURE.json`: new source
  retrieval records, including failed SI access.
- Original 68-row ledger and primary local source pins are recorded in the JSON.

**No new scoring panel was launched.** The next useful data work is the
LanM mutant table/construct completion and the exact GGR neutral-mutant table;
PqqT needs site evidence. A calculation plan must specify eligible constructs,
physical states and grouping before new preparations or scores are authorized.
