# Root-ready donor CP scout and finite panel

**No native calculation or submission occurred during preparation.**
The exact22-cell inventory is frozen in the matching prepared_v1/PANEL.json.
The scout is `cc__acetamide__La__complex`:10centers/201AOs/40electrons, RHF,
NoFrozenCore, conventional CCSD(T), exact matched-dry coordinates, effective
singlet. Use4MPI/one thread each, standard_expand/defaultmemory/time.

## Immutable bindings

- Workspace panel: `workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/donor_cp_v1/prepared_v1/PANEL.json`
- Panel SHA256: `1a4bca36f2e5f0d0e223b235b0915708a455627444bcb355fc7cbcba2930a701`
- REFERENCE.py SHA256: `f4baaeddce09cb26551e845f05e7f3c1bfff604a84fc02c666908ada2420dea2`
- SPEC.md SHA256: `697a9d9037901a40be81acb0b4a8f6158ebccbac2a9c8bdf9bc15a00428e40eb`
- validation_v2/VALIDATION.json SHA256: `b1c137e222e04eb6e4ee59197a571aa1f4643e49e1ff8830412d6505db9d54f8`

Use the owned RUN.sbatch wrapper with absolute arguments:

    REFERENCE.py run --panel PANEL.json --cell cc__acetamide__La__complex --out NEW_EXECUTION_DIRECTORY

Actual grant determines MaxCore; no memory/time override, duplicate input,
automatic retry or solver fallback. Review measured scout cost and native
admission before releasing the remaining11CC cells. The ten declared DFT
fragment cells are scientifically independent of the CC scout outcome, though
root decides their release. All four existing complex results are reused.

## Preparation and parser evidence

Static preflight checked every center, charge, electron count, ghost mask and
basis count from actual deck text. Acetamide complex/ion/ligand counts40/8/32;
acetate counts40/8/32 with ligand charge−1. Basis functions remain201/189 in
all CP roles,154/142 for the two real ligands. Metal ghost ECP is removed while
its47 orbital functions remain. D4 is confined to real-ligand cells and reused
real complexes; no ghost dispersion is admitted.

Parser validation passed45records: four real donor PBE0 outputs (including
all8/10basis centers), six real canonical water outputs,22positive controls
and13failure controls. The final center/basis beyond the water-era four-center
limit, anionic charge, correlated-electron count, frozen core, CC residual,
ghost ECP and ghost D4 are explicitly covered. No native preflight was repeated.
Validation used one local thread and no molecular calls. Validation_v1's
failed negative-control assumption about the final atom identity is retained
with its source snapshot; the native adapter and decks were unchanged.

## Collection and interpretation

The runner includes a complete22-result collection path:

    REFERENCE.py collect --panel PANEL.json --out NEW_RESULT.json RESULT1 ... RESULT22

The collector uses the pinned completed water CP baseline and compares each
metal's donor-minus-water interaction, then Dy-minus-La. Electronic PBE0 and
PBE0-D4 remain separate. The predeclared0.25kcal/mol diagnostic resolution
applies to each donor's paired PBE0-D4-versus-CC exchange difference. This
does not establish diffuse/metal-basis, frozen-f, spin-orbit, force or ensemble
accuracy. Keep the original native SCF stopping logic and raw residuals.
