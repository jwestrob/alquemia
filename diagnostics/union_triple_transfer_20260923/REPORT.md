# Three-source unions lose two supported triple calls

The practical three-source union retains **92/94 supported triple calls**, with
two new inconclusives from the Ca-associated A8R3S4 control. Released static and
the ten-fold precision candidate retain94/94 on these same declared triples.
There are no wrong calls, and all six old preparation failures remain unavailable.
The favorable canonical25/crystal3/stress pilot did not guarantee fold-transfer
fidelity. Do not replace the ten-fold candidate with this policy on these results.

## Frozen comparison

| Method | Correct | Wrong | Inconclusive | Unavailable | Total triples |
|---|---:|---:|---:|---:|---:|
| Released static context | 94 | 0 | 0 | 6 | 100 |
| Ten-fold union/adaptive, old strict stopping | 91 | 0 | 0 | 9 | 100 |
| Ten-fold union/adaptive, qualified numerical candidate | 94 | 0 | 0 | 6 | 100 |
| Three-fold union/adaptive, same numerical candidate | 92 | 0 | 2 | 6 | 100 |

All four three-of-four noncanonical La-fold subsets per protein remain included;
none was chosen after scores. These100 overlapping triples represent25 consumed
reference proteins, not100 independent validation proteins. Across125 distinct
source/context pairs (98 distinct folds), released and ten-fold precision methods
give124correct/1inconclusive; the new policy gives121correct/4inconclusive. Both
mathematical and operational pool policies agree on these counts.

The new three-fold canonical25 reference was frozen before transfer, using only
the original designated canonical geometries. Its Ca_max−405463.71230032056 and
La_min−405456.3869114709 happen to equal the ten-fold numerical reference because
unchanged canonical sources determine the extrema. No boundary was refitted,
widened or padded. Crystals and the separate Ca-conditioned A0A3 stress source
remain excluded from fitting and transfer denominators.

## The observed regression

Both inconclusive triples are A8R3S4:triple_008 (samples1/2/3) andtriple_010
(samples1/3/4). Their median is−405460.145287161 model kcal/mol,3.5670131596 above
the frozen Ca edge. This is not a rounding-sized boundary crossing. They share
the same168-atom union, reduced from174 atoms; the La-state charge remains−2.

Three source/context pairs newly abstain: A8 sample1 with168 atoms, sample3 with
168 atoms and sample1 with149 atoms. The last does not change its triple median.
The fourth individual inconclusive, C5AXV8 La sample3, is inherited unchanged.
The two offending168-atom A8 origins are Ca-supported before accommodation;
selected accommodation shifts R toward La by+9.0227/+8.0513 model kcal/mol. Thus
the failure concerns response in the smaller contexts, with composition, cavity
and possible mode-selection changes together. It does not establish a unique
missing-residue mechanism or a biological affinity change.

All61 triples with exact ten-fold membership replay identically. Among33 changed
triples, structural score ranges shrink for12 and grow for21; median range grows
from9.3163 to10.7915 model kcal/mol. These are correlated descriptive comparisons,
not a statistical validation claim. Full unrounded components, selected work,
each member's value and strict missing-member accounting are in
[the comparison](../../workspaces/union_triple_transfer_20260923/COMPARISON_v1.json).

## What actually ran

Exact receipt/map/state/profile checks accepted70 complete pool reuses. For55
changed contexts,36 paired origins reused actual compatible original-local
energies/forces;19 required38 MACE and76 native GFN2 calls. All110 single-start
searches,110 cross-MACE and440 candidate GFN2 calls completed. Search work used
1669 fresh MACE evaluations: **1817 fresh MACE and516 GFN2 total**, zero scientific
failures or repeated completed chemistry. No DFT, new folds, water/protonation
changes or additional Ca stress probes ran. Seventy-five searches retain boundary
flags; these are admitted finite candidates, not unconstrained minima.

Each metal sees the same origin/Ca-proposal/La-proposal pool. Native OMOL supplies
the vacuum energy and proposal forces; fixed-state native GFN2 ALPB−vacuum
supplies the solvent correction. E(Ca)−E(La) is an electronic model contrast,
not a complete binding free energy or equilibrium population. The frozen
four-mode rule,ftol1e−8 Hartree-equivalent,200 iterations,±0.8rad/0.8Å limits and
native300K/MaxIter500 scalar recipe are unchanged.

Two CPU allocations stopped in preflight before GFN execution because a second
host reproduced four selector diagnostic fields with differences up to1.066e−14.
A read-only replay confirmed identical selected subspaces and every other task
field for all55 pairs. A narrow1e−12 comparison for those four diagnostic floats
recovered the unchanged scalar manifests. Original frozen selectors, coordinates,
forces, parameters and receipts remain intact. See
[the technical recovery](SELECTOR_COPY_RECOVERY.md); this is not an energy or
optimization-tolerance change.

## Cost, tests and use

Whole allocations consumed **31,153 CPU-seconds and644 GPU-seconds**, including
both failed3-second CPU allocations and the one-CPU metadata replay. Nested runner
accounting is not added again. Local preparation/reuse audits, tests and comparison
are additional; the main reuse audit took173.195s wall with four workers. These
figures exclude historical reused pools and folding, so they are not a fully fresh
scanner timing. Maximum search cost was27 fresh MACE evaluations; median was15.
See [costs](COSTS_v1.json) and [scheduler receipts](SACCT_FINAL.txt).

Ten real-fixture tests pass in16.704s, zero skips. They verify125/100 denominators,
all61 exact replays, strict median/missing-member propagation, actual pool algebra,
source/state/protocol identities, the host-copy repair and rejection of corrupted
real records. Molecular execution is separately evidenced above. No scientific
executable was replaced by a fake result. [Commands](COMMANDS.md) replay saved
results; production and explicit historical methods remain unchanged.

Next: audit existing A8 membership/components/candidate selections and native
solver traces across all four triples. No additional chemistry, cutoff, threshold
or model change is authorized by this completed transfer report. The unknown-
protein three-source preparation API is useful plumbing, but this result does not
qualify an automatic scoring rollout.
