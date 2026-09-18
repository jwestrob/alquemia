# GGR discrepancy: distributed field and source self-solvation

The saved-output decomposition reproduces all three GGR structures and their
executed rigid variants. All1072 component/aggregation/rotation checks pass;
maximum residual1.50e-13 kcal/mol. Two real-output tests pass in2.463seconds.
No new quantum, MACE, density query or solver call was needed for this audit.

Of the +13.421788 kcal/mol direct-field shift (2FW0 minus2FVY), only+1.399594
comes from residues containing quantum-source atoms or cap projection anchors.
The remainder is distributed among unmodified protein residues. The largest
individual contributions include opposing ~3.11kcal terms from residues139
and137, and -2.88 from203. Several more remote charged residues contribute
around1kcal each. These terms are an additive accounting decomposition of this
functional, not independent residue binding energies or a unique causal test.
The induced contribution is-3.206576kcal. Full per-atom/source-residue/radial
records, unmatched IDs and every boundary charge redistribution remain pinned.

An additional analytical decomposition of the actual native GK expression
splits its+12.600209kcal shift into **+18.129021 source self-solvation** and
**-5.528812 source/environment cross-solvation**. This uses actual projected
charges, native effective radii and whole physical geometry. GK source moments
are monopoles, so the source self term is the native screened monopole kernel.
The residual cross term includes the environment multipoles. All five GGR/alpha
analytical rigid checks pass. Independent source-only native evaluation is
prepared separately and is not yet claimed as executed or validated.

All three GGR metal effective Born radii are near30A. The pinned native
`born.f:tanhrsc` explicitly approaches30A; OpenMM's reference implementation
uses the same bound and formula. Thus this is a property of the selected model,
not demonstrated evidence of a porting bug. Small input-radius sensitivity does
not by itself validate a saturated effective-radius approximation. Do not tune
or disable rescaling merely because a result is inconvenient.

Sources: native Tinker commit87050685eff8840d312e2a332cc82c33f63c7c3d,
`source/esolv.f:egk0a` and `source/born.f:tanhrsc` (exact file pins in the GK
accounting config); [OpenMM reference implementation](https://github.com/openmm/openmm/blob/8.5.1/plugins/amoeba/platforms/reference/src/SimTKReference/AmoebaReferenceGeneralizedKirkwoodForce.cpp).

Localization used13.066341wall/12.025296CPU seconds. The additional GK
accounting used8.335178wall seconds; its full CPU/RSS receipt is in the output.
Neither is a new predictor. Baseline and the failed frozen hybrid are unchanged.

## Next experiment

Verify source-only and empty-cavity energies directly against the unchanged
native solver, as declared in GK_SOURCE_NATIVE_VALIDATION_PLAN.md (15static
calls, no new quantum). If accounting holds, test a physically justified local
solvation/reference improvement. The present data do not justify changing
labels, omitting2FW0, retuning radii or attributing all error to cut bonds.

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/mace_density_localization.py --config workspaces/mace_omol_20260917/ggr_density_localization_config_v1/config.json --output workspaces/mace_omol_20260917/ggr_density_localization_review_v1
```

This reconstructs archived outputs only. Use a new output directory.
