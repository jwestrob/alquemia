# Nikasha in the integrated PLM manuscript

## September23 development update

The released method below remains available. The completed fixed-pocket plus
adaptive candidate now provides a real individual-structure improvement:
205correct/0wrong/1inconclusive/19unavailable across the frozen225 sources.
On205 shared sources, released201/2/2 becomes204/0/1. Both previous wrong calls
are corrected; C5AXV8 gains an inconclusive call. All25 canonical references and
three consumed crystals retain their expected classes under the candidate's own
frozen canonical-only bands. Strict aggregate coverage declines, and structural
spread is not uniformly reduced. These are repeats of25 consumed protein groups.

The candidate fixes surrounding fragment membership across the declared ten
source folds, selects four source-mapped angular motions using native MACE
forces, and proposes bounded motions independently for Ca and La. Both metals
then compete over the same origin/Ca-proposal/La-proposal geometry pool using
the existing composite expression below. This is finite local energy selection,
not an equilibrium ensemble or exact composite-energy minimization. Proton and
water inventories stay fixed; no entropy term is added.

Report this as a supported research improvement while its practical input and
numerical profile are qualified. Do not assign its ten-fold membership evidence
to the three-fold scanner automatically. Current PLM exports still contain their
original DFT results. [Completed comparison](../union_adaptive_20260923/TRANSFER225_REPORT.md),
[editable three-method figure](../overnight_discrimination_20260923/FIGURE.md).

The subsequent uniform stopping-policy candidate completes the same225-source
panel with206correct/0wrong/2inconclusive/17unavailable. On207 common sources,
released203/2/2 becomes205/0/2. It retains both error repairs and recovers two
optimizer failures; Q4W6G0 gains an inconclusive call. Strict aggregate coverage
recovers to23 La4,22 Ca5 and21 balanced groups, all correct; all94 available
La triples are correct. The new profile has its own canonical-only reference.
Large nativeGFN2 component changes at nearly identical geometries prevent a
claim of numerical equivalence or uniformly improved robustness. Its new
execution costs165472 allocated core-seconds and2312 requested GPU-seconds,
excluding reused origins/pools and prior qualification. The separate fresh
ten-A0A3-source comparison completes10correct versus released8correct/1wrong/
1inconclusive, without archived energies/forces. Including fresh preparation,
same-host32CPU/H200 allocations take358s versus451s. Batching, model-load count
and scalar solver parallelism differ; this measures the actual execution paths,
not an isolated scientific-component speedup. All160GFN/332MACE evaluations
succeed; eleven tests pass. [Fresh report](../pqq_union_execution_20260923/REPORT.md).

A separately prepared three-La-structure-context pilot retains25 canonical and
three crystal calls and repairs the difficult A0A3Ca1 probe using169 atoms rather
than193. Full transfer subsequently gives92correct/2inconclusive/6unavailable
triples, versus94/0/6 for both released static and ten-fold precision scoring.
Both newly inconclusive triples are the same Ca control, A8R3S4. The pilot does
not qualify the reduced-context policy for routine use. Its general preparation/
execution adapter is preflighted on existing PLM inputs, with molecular launch
held. No default or PLM result has been replaced.

Uniform two-pass native continuation preserves all25 canonical and3crystal
old-band calls and repairs the specified Q4 fold abstention, but five canonical
pools fail numerical settling. Therefore its new reference is unavailable; this
is a numerical diagnostic with a useful repair, not a calibrated release.

## Current supported method

Describe Nikasha as an affordable, structure-sensitive La/Ca **class discriminator**.
The released PQQ route is `fast_PQQ_OMOL_GFN2_ALPB_v1`; preserved native
r2SCAN-3c/CPCM DFT is an explicitly selectable comparator. The new name adds the
thin `scripts/nikasha` entrypoint and does not change historical protocol IDs.
Use [the current guide](../../docs/NIKASHA.md) for runnable source selectors,
single-source scoring and the supported three-source summary.

1. Pin protein sequence, actual source fold/crystal, assembly and metal/PQQ/donor
   selectors. Preserve the released physical source preparation, cap mapping,
   protonation and water state. Retain unsupported preparations in denominators.
2. Evaluate each endpoint with its own metal identity and charge at the same
   prescribed coordinates. The native float64 OMOL checkpoint is
   `MACE-omol-0-extra-large-1024.model`, SHA256
   `9b64b4fd5153ca578c694abc57806d8111050de6ff652e695c9b525bc4d36469`.
   Software/environment identities and exact inputs are pinned in each manifest.
3. Add the matched native GFN2 water-ALPB minus vacuum energy, using ORCA6.1.1:
   `E_M = E_OMOL,vac,M + E_GFN2,ALPB,M − E_GFN2,vac,M`.
   The electronic-state and solver settings belong in the exact technical
   supplement copied from the released input recipe. This is a composite energy
   descriptor; no molecular-dynamics, FEP or binding entropy is implied.
4. Form `R = E_Ca − E_La`, converting each component once. LargerR is more La-like
   on this protocol's scale. Use its frozen canonical calibration; do not inherit
   the older DFT aquo gauge or interpret a universal zero as a preference threshold.
5. Return correct/wrong/inconclusive/unavailable separately for labeled controls.
   For environmental proteins return the computed support category, source spread
   and unresolved reasons. A fold-conditioned hypothesis is not a measured metal
   occupancy state; the output is notKd, a probability or a substrate assignment.

The measured historical release cost is about42seconds/source on oneH200 plus
32CPUs, excluding folding. Do not substitute this for newly measured search costs
or claim a matched-hardware DFT speed ratio from different allocation receipts.

## Reference results currently supportable

Separate canonical calibration and consumed crystal transfer from structural
replicas. The completed challenge includes225 noncanonical source structures
from25 protein groups, with Ca- and La-conditioned folds retained separately.
On the207 sources covered by both methods, the released context composite has
203correct/2wrong/2inconclusive versus DFT198/3/6. This is a developmental structural
robustness result, not207 independent biochemical observations or prospective
validation. All94 available three-fold subsets are correct for both methods;
those subsets are correlated. [Exact report](../nikasha_recovery_20260922/REPORT.md).

Report coverage, class direction and composition limitations explicitly. The
canonical calibration's motif/composition separation limits claims of incremental
electronic information. Claims of general metalloprotein affinity prediction are
outside this evidence. Research accommodation variants remain separate until
their actual reference transfer and cost justify selection.

## Environmental application and existing tables

The [joinable export](../nikasha_plm_export_20260922/REPORT.md) retains176 proteins
and200 original genes, current tree/domain identities, motifs, genome evidence and
the actual transcript exports. **Those delivered PLM scores are the existing DFT
protocol**, not a newly computed MACE cohort. A future candidate result must be
attached with its own protocol and source identities; do not relabel old scores.

Distinguish sample-wideTPM, the existing length/read-depth-normalized expression
profiles, and within-genome/campaign-adjusted correlations. No general
genome-normalized expression matrix was located, and none was invented.
Keep source-reference mapping, missing transcripts, catalog aliases and
gene/genome evidence strengths visible. Do not sum independent assembly mappings
into a community total or infer a replicated depth/month/site causal effect.

## Exact supplement and figure contents

- Pin released protocol/reference files, sources, model/software hashes, parameter
  templates, all preparation failures and per-endpoint execution receipts.
- Provide reference/source-level decisions plus strict protein-group summaries;
  show identical-coverage comparisons and all excluded/unavailable counts.
- Keep electronic contrasts, source variation and accommodation work distinct.
  Any geometry-response illustration must identify its real source and whether
  its metal preference is experimentally known.
- Use editable vector figures generated from the pinned result tables. Present
  PLM phylogenetic placement and transcriptional evidence alongside conditional
  Nikasha support, preserving the broader manuscript's ecological findings.

This outline describes what the existing evidence supports. It does not freeze a
research challenger as the selected manuscript protocol or claim that a new
MACE-based PLM cohort has already been calculated.

## Delivered scientific figures and research conclusion

The editable [structural comparison](STRUCTURAL_COMPARISON_CAPTION.md) shows
the matched207 results and strict100-triple denominator. The independent
[fixed displacement figure](FIXED_DISPLACEMENT_CAPTION.md) shows all nativeDFT,
MACE and composite differential responses for the two tested directions in all
three contexts. It validates those particular displacements, not the later
adaptive/joint proposals. Captions include biological and replication limits.

The later [actual accommodation example](ACCOMMODATION_EXAMPLE_CAPTION.md)
plots the cross-scored common pool and source-mapped metal–heteroatom distances
for the difficult A0A3 Ca-conditioned fold. It makes the additional local response
visible while separating it from the preceding pocket-membership change.

The September22 shared-pool and adaptive experiments remain preserved in their
exact [delivery record](../nikasha_recovery_20260922/DELIVERY.md). The successful
September23 combination described above supersedes the earlier conclusion that
no accommodation variant had demonstrated useful improvement. It has not changed
production defaults. Distinguish the completed accuracy result from the ongoing
qualification of a practical execution profile and smaller input inventory.

The [technical supplement](TECHNICAL_SUPPLEMENT.md) records the released model/
software identities, actual nativeGFN2 inputs, unit conversions, frozen bands
and exact completed example artifacts.
