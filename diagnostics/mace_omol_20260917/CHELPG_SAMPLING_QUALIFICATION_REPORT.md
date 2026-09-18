# Native fixed-density CHELPG replay qualifies

Both real GGR2FVY endpoint replays pass the gates frozen in
[the sampling plan](CHELPG_SAMPLING_STABILITY_PLAN.md). No SCF optimization,
gradient, MACE inference, geometric change or new biological score occurred.

| Check | Ca | La | Allowed |
|---|---:|---:|---:|
| Energy change, Hartree | +2.356682671e-8 | -9.358018360e-9 | absolute1e-7 |
| Largest printed charge change,e | 0 | 0 | 2e-6 |
| Largest density potential change,au | 7.818745651e-14 | 9.230116671e-14 | 1e-8 |

Native ORCA6.1.1 confirms matching input geometry/basis, imported occupations,
NoIter and one fixed-orbital energy evaluation. Both terminate normally with
return code0. The utility queries used the newly written densities at the
original exterior validation positions; saved-orbital reuse alone was not
accepted as proof of identity. La Def2-ECP46-core-electron accounting passes.

The generic runner deliberately retains `scf_converged=false` and its failure
status: no iterative SCF was requested. Thus Slurm1201173 reportsFAILED even
though both property calculations completed. The separate NoIter contract
qualifies those actual results; nothing was relabeled as a converged SCF or
rerun to hide the status. Density query job1201174 reportsCOMPLETED.

Actual costs:1201173 41wall seconds,2624allocated core-seconds,
1013.593CPU seconds,3951892KiB peak RSS.1201174 5wall seconds,
320allocated core-seconds,6.196CPU seconds; Slurm reports0KiB batch RSS,
which is not a meaningful measured memory peak. Both allocations were64CPU
despite requests32/2. No GPU. Per-command time-v receipts are retained.

All five real-fixture qualification tests passed in1.348seconds, no skip.
The declared38remaining property replays and38density identity queries can
now proceed; the initial two are reused exactly. The production baseline,
original parser, charge convention and scientific thresholds are unchanged.

Primary artifacts under `workspaces/mace_omol_20260917/`:
`chelpg_sampling_qualification_v1/manifest.json`,
`chelpg_sampling_population_recollection_v2/result.json`,
`chelpg_sampling_density_identity_v1/collection_job_1201174.json`.
The v1population collection is retained; v2 adds receipt-artifact verification
without changing scientific values or executing another calculation.
