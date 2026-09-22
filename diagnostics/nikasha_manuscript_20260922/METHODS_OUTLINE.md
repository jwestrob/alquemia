# Nikasha in the integrated PLM manuscript

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

The completed shared-pool and adaptive experiments are research results, not
additional production-method components. None earned routine promotion. Keep
the paper's supported scoring method separate from these failed extensions;
their exact [delivery record](../nikasha_recovery_20260922/DELIVERY.md) preserves
all failures, methods and costs.

The [technical supplement](TECHNICAL_SUPPLEMENT.md) records the released model/
software identities, actual nativeGFN2 inputs, unit conversions, frozen bands
and exact completed example artifacts.
