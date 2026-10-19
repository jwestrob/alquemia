# Saved development comparator

Source-only preparation under the root-owned
[frozen design](../whole_base_delta_design_v1/SPEC.md). No calculation has run.
`CONFIG.template.json` deliberately leaves `candidate_implementation` and
`candidate_result` unbound. Root binds them only after the actual native/assembly
parent exists. This comparator executes only inside Slurm; there is no submitter.

```text
PINNED_PYTHON COMPARE.py --config BOUND_CONFIG.json --output FRESH_DIRECTORY
```

The pinned interpreter is the existing `lanm_qmmm` Python from the admitted
scaffold audit. Use one allocated CPU and its original single-thread environment;
root owns the launcher and actual allocation/cost receipt.

## Producer contract

Agreed with the whole-base-delta source owner before comparison implementation:

- Parent protocol `lbt_whole_GFN_regional_vacuum_delta_v1`, stage `collect_7CCO`,
  status `complete_7CCO` only with all four development rows complete.
- Exact `implementation`, `config`, `inventory`, `preflight_result` pins;
  `new_molecular_calls=0`; preserve parent results and cost/call evidence.
  Producer config's `design` binds the root SPEC; its `specification` is the
  separate native adapter spec.
- `rows` is the ordered eight-cell list from `FEASIBILITY_INPUTS.json`.
  First four statuses are `complete`; the four 7CCN statuses remain
  `not_requested`. An unavailable first-source parent yields no work triage.
- Completed row fields: `cell_id`, `source_id`, `geometry`, `metal`, `physical`,
  `positions`, `mapping`, `physical_ids`, `state` (whole/region, exactly the
  admitted preflight manifest), and `solver_exception_certificate`.
- `component_receipts`: `GFN_full`, `MACE_region`, `GFN_region`, `ALPB` pins.
  `native_energies`: `GFN_full_Eh`, `MACE_region_eV`, `GFN_region_Eh`, `ALPB_Eh`.
- `component_gradients_kcal_mol_A`: `GFN_full`, `MACE_region`,
  `negative_GFN_region`, `fixed_FF_ALPB`. Local terms are independently scattered
  through the exact archived cap map. Candidate gradient is the left-to-right
  sum in that order; base gradient is GFN_full plus fixed_FF_ALPB. Stored fields
  are `candidate_gradient_kcal_mol_A`, `base_gradient_kcal_mol_A`,
  `candidate_energy_kcal_mol`, `base_energy_kcal_mol`.

The continuation producer is still being prepared; the template does not invent
its source hash or a completed result. Source-level schema review must include
its eventual exact emitted fields before binding or running this comparator.

## Reuse and interpretation

The only imported project module is the frozen scaffold audit source. Only its
saved `saved_inputs()` and `describe()` functions run. Those preserve the exact
whole-PBE0/model raw receipt joins, producer-specific ALPB conversion order,
four-state physical identity, named terminal atoms, full256 force definitions
and original two directions. No MM execution, model load, evaluator, native
process or optimizer is called.

Native component energy differences are taken before conversion. The candidate
and base-only work are reported for both metals. The sole challenge triage is
absolute candidate-minus-PBE0 deformation work <=1 kcal/mol for **both** metals.
A scientific triage failure is a completed saved analysis with `eligible=false`;
it does not authorize a 7CCN challenge. An admission/replay failure leaves the
comparison unavailable and exits1. No partial-survivor work is reported.

`FORCE_TRANSITIONS.json` retains complete256 force, model-minus-reference error,
after-minus-before force change and error-change arrays. Original whole-MACE
before-Dy force and old-hybrid after-Dy force/work stay unavailable. Regional
component arrays are admitted saved producer outputs with exact internal
closure, not a new cross-host raw-cap bitwise-reconstruction assertion. Raw
native GFN/solvent identity and full-term conversions are checked separately.

SCF/electronic-state evidence and all reference caveats stay visible. There is
no candidate derivative qualification, force-accuracy pass, Hessian, optimizer,
affinity or protein-transfer conclusion. Pending7CCN reference outcomes are not
read or used by this first-source development comparator.

Validation performed so far: AST syntax parsing and actual saved schema reads
only, without importing this source or executing numerical analysis. Future
real Slurm execution remains required.
