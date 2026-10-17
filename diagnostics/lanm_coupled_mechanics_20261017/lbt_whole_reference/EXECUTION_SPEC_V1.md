# Whole-LBT before-La scout: execution declaration v1

This implements only the first scout in immutable `SPEC.md` (`02b45065...`) and
`INVENTORY.json` (`48d9f40f...`). Root owns Slurm execution and the decision to
launch. Preparing this source does not submit or start a reference calculation.

## Exact scope and stages

`REFERENCE.py` prepares all four declared whole256 states at the fixed 7CCO La
hybrid before/after coordinates. Its hard-coded execution allowlist contains
**only `before_La`**. There is no command argument for another task, no remaining
panel stage, optimizer, parameter sweep or automatic retry. A later finite
continuation needs a separately reviewed execution version after actual scout
admission and cost review. No numerical derivative calls are added here.

Proposed commands inside root-owned Slurm allocations:

```text
python REFERENCE.py --stage preflight --config CONFIG_PREFLIGHT.json --output FRESH_PREFLIGHT
python REFERENCE.py --stage scout --config BOUND_SCOUT_CONFIG.json --output FRESH_SCOUT
python REFERENCE.py --stage collect --config BOUND_COLLECTION_CONFIG.json --output FRESH_COLLECTION
```

Preflight and scout require exactly 12 actual MPI slots, one node, one CPU per
rank and memory derived from the actual allocation by the existing pinned
runtime helper. The collector performs saved-data analysis only and may use one
CPU; it still requires Slurm. No stage submits a job or changes scheduler policy.
Output directories must be new. `BOUND_SCOUT_CONFIG.json` must contain the actual
passed preflight RESULT pin. The collector additionally binds the actual scout
RESULT; a missing/changed parent fails before any molecular call.

## Zero-call preflight

1. Verify the immutable inventory, source/terminal/search/pool joins, native
   asset/helper/Python hashes, and actual admitted nonconverged 13017 endpoint.
2. Reproduce the exact 256 physical IDs/coordinates, including the three original
   terminal completions H2, H3 and OXT. Materialize the cell's La/Dy element at
   physical target index255; never use hybrid core index0 or the 172-atom XYZ.
3. Derive all-electron count from actual nuclear charges and total charge−2;
   subtract only the declared La46/Dy55 ECP core. Require positive even explicit
   counts, matching La/Dy valence counts and effective singlet populations.
   Preserve physical La1/Dy6 metadata separately. The earlier 48-atom/200-electron
   and 172-atom prospective counts are not imported.
4. Check pinned metal basis/ECP/auxiliary blocks and CHNO light-atom coverage
   under the existing def2-TZVP/def2-J route. Write exact round-tripping full XYZ
   files and the PBE0-D4 vacuum `EnGrad` decks with `* xyzfile -2 1 physical.xyz`.
   No external charges, caps, MM, GFN, CPCM/SMD, extra dispersion, changed SCF
   iteration limits or geometry optimization are present.
5. Replay the three archived full fixed-charge ALPB receipts: raw native scalar
   fields, hash-linked XYZ/charge/physical provenance, binary command and native
   return code, and complete finite256-atom gradient. The fourth after-Dy solvent
   stays explicitly unavailable; no solvent computation is performed.
6. Freeze the full before-to-after displacement chord and saved initial Phe1 CA
   hybrid-force direction before reference outcomes. Reject three meaningful
   negative state controls: charge−1, missing target atom and wrong target metal.
7. Save all four state/deck signatures, complete denominator, limitations and
   allocation evidence. Preflight reports zero molecular calls. It does not
   establish native SCF convergence, basis performance or force accuracy.

The scout reconstructs this preparation in a fresh directory and requires exact
task signatures, state, directions and solvent replay against actual preflight.
Only then can the existing native runner receive the one before-La task.

## Native reference and admission

The fixed recipe is the original vacuum
`PBE0 D4 def2-TZVP def2/J RIJCOSX NoAutostart DefGrid3 VeryTightSCF EnGrad`,
PModel guess, exact La lcecp1-TZVP/AutoAux-J and ORCA6.1.1. The default native
iteration/convergence settings are unchanged. The narrower whole-state/deck
wrapper replaces the old helper's hard-coded charge−1; the existing analytic
parser is reused with the newly derived whole state.

Admission retains native normal termination and convergence, actual charge and
derived whole electron count, effective multiplicity, RHF/ECP identity,
ECP/dispersion analytic-gradient evidence, finite
gradient and exact atomic order. The original parser's energy correspondence
tolerance is `1e-8 Eh` and coordinate correspondence is `1e-6 Å`; neither is
loosened. No count from a different electronic region is imported.
The input recipe, native PModel and vacuum scope, `ConvCheckMode=Total+1el-Energy`,
`ConvForced=0`, actual 12-rank memory policy and ORCA binary must match.

Keep every printed residual and gap/state warning. `all_printed_residuals_pass`
remains an observation under the original native admission semantics. A native
failure or parser failure is retained as unavailable with raw artifacts and no
retry. Runner invocations and actual created native outputs are counted
separately; reused ALPB and every other molecular component count remain zero.

## Saved collector and interpretation

The scout and independent saved collector produce the same before-La force
comparison after replaying native receipt/gradient admission. They use the exact
same raw ALPB as the existing whole-MACE and hybrid point. The reference contains
whole vacuum PBE0 plus that solvent term, with no retained MM/GFN internal term.

Save the full256 force table, Phe1 table/CA vector, global residual distributions,
reference-minus-model force errors and projections along the two frozen
directions. Report raw vacuum reference forces as well as composite forces.
No absolute energy subtraction between methods, force-accuracy pass threshold,
stationary-minimum claim or metal preference claim is made. Three reference
cells remain `not_requested`, and before/after work remains unavailable.

The original terminal uncertainty, vacuum dianion/electron-binding concern,
post-added solvent limitation, frozen-f/SOC omission, nonconvergence and
non-uniqueness of MM-parameter attribution are present in every task/result as
machine-readable caveats. The independent PBE0 calculation is not represented
as the exact MACE training functional. All numerical analysis runs only on Slurm.

## Current state

Source and AST review only; actual preflight, whole-state electron counts,
native calls and new force comparisons have not run. Root binds actual execution
pins after static review. The original proposal and Mex inventory remain intact.
