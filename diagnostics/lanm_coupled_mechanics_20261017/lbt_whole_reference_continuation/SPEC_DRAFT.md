# Whole-LBT fixed-deformation continuation — unbound source draft

17 October 2026 UTC. This prepares the remaining execution path in the existing
four-cell proposal while scout 13051 runs. It does not authorize execution,
establish scout admission, or replace the post-scout cost review. The original
SPEC.md, EXECUTION_SPEC_V1.md, REFERENCE.py, inventory and actual preflight 13049
remain immutable. Root owns integration and all submissions.

## Scientific question and fixed scope

Does whole-peptide PBE0 support the restoring forces and deformation work assigned
to the exact 7CCO La-proposed before-to-after motion by the two existing models?
The original selection is unchanged: all 256 atoms, two fixed coordinate arrays,
and both La and Dy on each array. No source or coordinate is selected from a new
reference outcome. Both endpoints are nonstationary.

Scout `before_La` already runs as 13051. At most three new native EnGrad calls
remain: `after_La`, `before_Dy`, `after_Dy`. Each run executes one declared cell
in a fresh directory. No scout repeat, retry, solver change, optimization,
finite-difference call or additional molecular component is included. Failed
cells remain unavailable in the full four-cell denominator.

Use the exact existing vacuum PBE0-D4/def2-TZVP/def2-J/RIJCOSX, DefGrid3,
VeryTightSCF, PModel and NoAutostart recipe, native executable, metal basis/ECP
bytes and analytic-gradient parser. Preserve charge -2, target index 255,
physical La/Dy multiplicities 1/6 and effective frozen-f singlets. Whole-system
state and atomic order are checked against the original actual preflight.
Twelve MPI ranks and allocation-derived memory are unchanged; no scheduler
memory, time or priority overrides are introduced.

## Required evidence before execution review

The draft must fail closed unless its configuration pins all of the following:

1. The original adapter and bound scout configuration, actual preflight 13049,
   and their original scientific inputs.
2. The actual completed 13051 result and its independently replayed original
   saved collector. Native admission must pass; all residual and physical-state
   caveats remain, without a newly invented all-residual acceptance rule.
3. A root continuation decision, written only after reviewing the actual scout
   forces, native outcome and measured cost. It binds the exact draft source,
   specification, scout/collector/preflight pins and cost receipts, permits
   only the three remaining tasks and retains a three-new-call ceiling.
4. The completed whole-model collection 13142, supplying separately admitted
   after-Dy ALPB and the available same-coordinate model comparisons.
5. For each native run, the actual passed zero-call continuation preflight under
   this exact implementation/configuration binding.

No runnable configuration or continuation decision is created in this source
preparation phase. Static source review can occur now; execution-version review,
real-artifact Slurm checks and scheduling remain after actual scout admission
and cost review. A source draft passing syntax checks is not a tested runner.

## Immutable replay and supplemental solvent

Replay original 13049 with its original three-available/one-missing solvent
inventory and four exact states/decks. Do not edit or reinterpret that parent.
Then admit after-Dy ALPB from the independently completed 13142 collection as
separate supplemental evidence, with full coordinates, physical IDs, metal,
fixed charges, native output and 256-atom gradient checked. This performs zero
new solvent, MACE, GFN or MM calls.

The original before-La native result is replayed and reused. Remaining singleton
manifests are explicitly distinct from the original scout allowlist. Every
stage requires Slurm; source creation and syntax inspection alone run on login.
No stage submits jobs or waits inside an allocation for an unfinished parent.

## Collection and interpretation

Preserve all four reference slots and inspect actual raw native artifacts under
the unchanged admission policy. A missing or failed required cell leaves its
dependent work comparison unavailable. Compute each metal's after-minus-before
vacuum energy in native Hartree before conversion; add the matching fixed-FF
ALPB difference separately. Report Dy-minus-La deformation work and compare
only matched deformation work across methods, never absolute method totals.

Report available full force tables, Phe1/CA vectors, residual distributions and
projections along the original frozen full displacement and initial Phe1 CA
force directions. Preserve missing before-Dy whole-MACE forces and after-Dy
hybrid evidence explicitly. Do not manufacture a full comparator matrix or
launch calls to fill those diagnostic gaps.

The original 1 kcal/mol response resolution is a reporting resolution, not an
uncertainty estimate or a force-accuracy threshold. Preserve all printed SCF
flags, unresolved terminal chemistry, vacuum-dianion/electron-binding concerns,
post-added solvent limitations, frozen-f/SOC omissions and the absence of an
independent numerical DFT gradient qualification. This finite comparison can
inform model choice; it cannot alone establish a LanM affinity or preference
classifier, unique MM parameter defect, or production promotion.

## Required validation after the scout

Check exact original-state/deck/solvent replay, actual scout/collector replay,
supplemental solvent provenance, singleton allowlists and complete accounting.
Meaningful negative controls must reject missing or changed scout evidence,
an altered cost/continuation decision, an attempt to repeat before-La, and a
wrong-coordinate or wrong-metal supplemental receipt before a molecular call.
The actual four-cell collector must independently replay the raw results and
retain every unavailable slot. Record actual scheduler CPU costs separately
from application time and never count reused chemistry as new work.

Status: source preparation only; no new native or numerical job is submitted.
