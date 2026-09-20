# Does physical donor accommodation improve structural robustness?

2026-09-20. Root's contained continuation under Jacob's explicit overnight goal
and permission to pursue promising directions. The 30-source proposal experiment
has completed and its outputs have been inspected. This is development on consumed
structures, not a blind test. No production changes are authorized by this plan.

## Question and motivation

Native MACE proposals plus actual solvent-aware energy selection reproduce the
large differential accommodation effect from the full composite search at much
lower cost. All reference classes remain ordered, but Q9Z4J7 becomes inconclusive
under the previous method's bands. Calibration gap alone is not an accuracy gain.
The useful next test is whether the same physical rule makes calls more robust
to source geometry on the existing known-class reference fold challenge.

## Fixed population and algorithm

Use every one of the 225 primary, noncanonical cases in the original all250 source
inventory. The 25 canonical coordinates are calibration only. Keep both folding
metals separate, all 25 protein groups and all unsupported inputs. At this point
208 primary preparations are supported and 17 unavailable; one supported context
also has a failed archived solvent endpoint. No replacement, rescue SCF, new fold,
new protonation, water, context membership or sample selection.

Prepare the same Glu chi3 and actual extra-Asp chi2 mappings for all supported
primary cases, reusing exact compatible mappings already prepared. Other residues
and all inactive coordinates stay fixed. Mapping preparation makes no energy calls.
Use the unchanged `native_MACE_proposal_primary_composite_selection_v1` objective,
starting point, bounds, optimizer, overlap checks and two-member energy selection
from PROPOSAL_PLAN.md. No warning- or result-directed subset. Missing q0 or proposal
energy remains unavailable; the original score is not a fallback success.

The finite manifest contains at most 416 metal-specific proposal tasks and at most
832 new primary native-GFN2 singlepoints. Exact archived q0 energies are reused;
unavailable q0 prevents that endpoint's optimization. Actual supported counts and
failures are retained. No new DFT, numerical gradients, Hessians or native minimum
claim. One warm H200/32CPU/200000MiB job proposes geometries; up to four disjoint
64CPU solvent manifests each use eight concurrent eight-rank tasks. Use existing
runners with actual task-path dry-runs and receipts. No time or project-cost cap;
the algorithm and source set are finite. No additional poses or SCF retries.

## Distinct calibration and evaluation

Report the original frozen bands unchanged as a transfer check. Also create an
explicit **new developmental reference** for this new geometry-selection method,
using only the already designated 25 canonical calibration structures from the
completed 30-source experiment. Apply exactly the previous calibration rule:
Ca-supported through the largest canonical Ca score, La-supported from the smallest
canonical La score, provided the same existing minimum-gap condition is satisfied.
Record all 25 source/result hashes, model/selection settings and class extrema.
This is a new calibration, not preserved old thresholds or test-set optimization.
No alternate cutoff, fitted coefficient, selected subset or re-calibration after
fold results. The three consumed crystals and two unknown PLM sources do not
determine bands. If calibration is unsupported, calibrated calls remain unavailable.

Compare the proposal descriptor with unchanged native/core/context/composite and
preserved DFT on identical source IDs. Report correct/wrong/inconclusive/unavailable
denominators for La-conditioned and Ca-conditioned folds separately; transitions,
raw contrasts, within-protein spread and class ordering. Repeat the already
declared strict La4/Ca5 summaries and every three-of-four La subset (100 correlated
triples), with no best-triple selection. Missing members invalidate that summary.
Distinguish increased coverage from accuracy and single-fold gains from grouped
robustness. Homologs and structural repeats are not independent biological labels.

The goal is useful robustness on known classes with affordable throughput. A
shift toward La for an unlabeled PLM source, a fitted calibration success, or faster
optimization alone does not satisfy it. A failed robustness test is retained.
