# Independent water-orientation validation assessment

**Recommend a separately declared four-cell water substage.** The frozen17605 water candidate can be tested before donor17621 finishes. No physical donor-fit dependency appears in these configurations: each contains one ion and eight waters, with types500/501,36,37. Neither carboxylate403 nor amide105 is present. The donor evaluator leaves these water inputs unchanged, and their ECP-PBE0 reference Hamiltonian is independent of the donor fit.

## Existing evidence to reuse

- Passing frozen water result `6027e319`, selected raw review `a3147adc`, and its exact six parameters.
- Actual17612 water Cartesian checks: all four dimer/aquo centers reproduce saved energies exactly; largest component discrepancy is9.0e−8 kcal/mol/Å against the unchanged1e−3 limit. The actual review is `3f7b3a9b`. No repeat center checks are needed solely to remove the scheduling dependency.
- Exactly four `validation_aquo8_orientation_{minus,plus}_{La,Dy}` geometries from reserved panel `eed99b0f` and native panel `75d3fc69`.
- Existing17605 central-aquo model anchors and original ECP central-aquo references `f1be7b98`/`c1ff6777`. No new anchors or fragments.

## Required explicit amendment

The current Stage A declaration and both adapters require a passing17621 donor candidate and all12 consumed donor numerical checks before any reserved execution. Leave those gates intact. A new root scheduling declaration must exempt only the four water IDs from those donor dependencies, substituting the actual passing17605 selection and17612 water numerical evidence. Use narrow versioned bindings with exact water input equality; do not supply a fictitious donor candidate or weaken the original16-cell gates.

Declare four unchanged reference E/gradient calls, four model E/force cells (eight analyze/analytic-testgrad programs), four MPI ranks per reference/default resources, and no new Cartesian energy panel. Freeze this before reserved readout.

## Preserve validation meaning

Use the existing full16 scorer with the other12 numerical records unavailable. Report only the water substage's two aquo-referenced energies, three La/Dy/paired orientation works, and six rigid-force fixtures. Keep the original1 kcal/mol energy and2/5 kcal/mol/Å force limits. Full Stage A remains incomplete/false with33 other comparisons unavailable; Stage B remains sealed.

Carry every water result and failure forward without repeating or dropping it. Later full16 integration must use the identical water Hamiltonian and satisfy the original donor gates. These four observations do not expose donor validation outcomes or justify changing the ongoing donor fit. If they drive water redevelopment, they become consumed development evidence; replay cannot restore independent-validation status.

Assessment receipt: `ASSESSMENT.json`, SHA `8979d95569ac54bde5b60b8d276ed30316778e9adc65948171e1866d32f68d87`. No gates edited, candidate bound, reserved E/F opened, scheduler query or native calculation performed.
