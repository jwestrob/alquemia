# Actual adaptive response of the difficult A0A3F2YLY8 Ca-conditioned fold

The figure uses the completed, uniform-precision common pool for
`a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1`. It is a deliberately
identified development example of a repaired reference error, not an independent
validation case or a claim that every protein responds this way.

**A.** Both metals are evaluated at the same three physical geometries: origin,
Ca proposal and La proposal. Bars show the native MACE plus native GFN2
ALPB-minus-vacuum composite, relative to each metal's own origin. The black
outlines identify the operationally selected row minima. Ca gains5.40234 and
La gains10.59938 model kcal/mol, so accommodation increases the paired Ca-minus-La
contrast by5.19704 model kcal/mol. The solvent correction partly offsets the
native MACE energy reductions. All raw components are retained in the exported
CSV. No cross-element absolute electronic energies are interpreted as affinity.

**B.** Actual mapped source-atom distances show what moved. The selector admits
Glu197 chi1, Asn285 chi1 and Asp329 chi1/chi2. The Glu197 and Asp329 carboxylate
orientations change, with different Ca/La responses; the shown PQQ atoms remain
fixed. All mapped N/O atoms initially within3.5Å of the metal are displayed.
This cutoff selects labels for the illustration only; it assigns neither bonds
nor coordination numbers. The source identity, atom indices and unrounded
distances for all three geometries are exported.

The context membership, chemical state and water inventory are identical within
this figure. The preceding context expansion is a separate intervention; on this
case it already repairs the original static call. This panel illustrates the
additional accommodation response, not proof that accommodation alone repaired
the error. Across the full panel, the combined method preserves both error
repairs with fewer new inconclusives than context expansion alone.

The Ca proposal reaches the prescribed0.8Å heavy-atom displacement boundary.
Neither proposal is claimed to be an unconstrained stationary minimum of the
composite Hamiltonian. MACE-vacuum forces generate proposals; exact composite
energies select them. These finite candidates are not an equilibrium ensemble,
and their energy differences are not binding free energies or entropy estimates.

## Artifacts and exact reproduction

Actual output: `workspaces/discrimination_transfer_figures_20260923/accommodation_example_v1/`.
Editable SVG/PDF, PNG, distances CSV, component CSV and hash-pinned receipt are
present. The PNG was visually inspected. Rendering also checks every relative
energy against raw matrix components and the archived selection work. No new
molecular evaluation or threshold fit was performed.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/plot_accommodation_example.py \
  --collection workspaces/slsqp_precision_20260923/pool_v1/collection_final.json \
  --case a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1 \
  --output workspaces/discrimination_transfer_figures_20260923/accommodation_example_v2
```
