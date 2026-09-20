# Test whether the folding metal controls the classifier

## Why this test

A single modeled coordination site can encode the metal supplied to the folding
program. The desired classifier should retain biological class across reasonable
structural samples, rather than merely report that modeling choice. Existing
reference predictions allow a controlled test without assigning labels to PLM
predictions or generating new folds. This directly tests a source of error that
the selected-reference calibration cannot expose.

## Fixed inputs and primary comparison

Use all250 already saved AF3 geometries for the25 canonical PQQ proteins:
five La-conditioned and five Ca-conditioned samples per protein, seed1 and the
same actual model build, sequence, PQQ input, MSA and templates. The curation
agent verifies these identities and pins the existing artifacts. Keep the25
canonical source matches as replay controls; the primary structural transfer
set contains100 other La-conditioned and125 Ca-conditioned samples. These are
structural repetitions of25 known calibration proteins, not250 independent
biological observations or prospective blind proteins.

Run the unchanged fixed-core and complete-context preparation policies on every
source, retaining any failures in the full denominator. No alternate source
selection, new geometry optimization, water, protonation policy or method is
introduced to obtain favorable labels. A rawCa source needs an explicit matching
metal selector; replacing the selected ion for paired scoring preserves its actual
source coordinates. Record source-metal identity separately from scored metal.

Evaluate nativeOMOL plus the same nativeGFN2 ALPB-minus-vacuum contribution in
both core and context representations. This is at most1000MACE endpoints and
2000GFN2 endpoints across all250sources and two representations. Reuse exact
compatible canonical outputs when justified; report reused versus fresh work
separately. Preparation runs in independent processes with pinned seeded behavior.
Use existing allocation-aware ORCA and nativeMACE workers, no new workflow engine.
No new DFT calculation is part of this first comparison.

## Frozen candidate summaries before energies

Report every individual source's contrast and the fixed existing representation-
specific calibration-band transfer, including errors, inconclusive outcomes and
unsupported preparations. Retain standalone nativeMACE and the twoGFN2 components
to distinguish a solvent effect from source-geometry effects.

Primary candidate descriptors, with no fitted coefficients:

1. Median contrast across the four unselected La-conditioned samples.
2. Median contrast across the five Ca-conditioned samples.
3. Equal-weight mean of those two medians, a folding-metal-balanced descriptor.

For each candidate pool require every declared member to be supported/scored;
an incomplete pool is unavailable rather than silently reduced to favorable
members. Report within-protein contrast ranges and change between conditioning
arms. An ensemble summary is a robust descriptor, not a thermal population or
binding free energy. Old numerical bands are an explicitly frozen transfer test,
not an automatically validated new calibration. No new threshold fit in this
phase. Native labels remain PQQ functional association.

Within each protein, verify core atom composition, charges, cofactor and water
state. Context membership/charge can change under the frozen geometric selection
rule; retain those changes and do not attribute their full effect to nuclear
response. Primary causal geometry interpretation uses the fixed core. Context
summaries are operational robustness of the existing selector/scorer combination,
not fixed-boundary environmental physics.

## Decisions and later work

The useful outcome is evidence that a cheap, declared treatment of alternative
real structures reduces wrong or inconclusive reference calls and conditioning
sensitivity while preserving coverage. A larger fitted calibration gap alone is
not success. If most sources already classify reliably, that limits the value of
ensemble scoring on this task and redirects effort to chemical states outside it.

Do not change production/default scores or rescore the PLM cohort in this phase.
Any subsequent calibration, DFT adjudication or per-site structural intervention
gets a recorded finite manifest within the authorized overnight research goal.
The separate source-to-score fast-mode release remains independent.
