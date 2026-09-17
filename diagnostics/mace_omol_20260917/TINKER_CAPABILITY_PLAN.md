# Next engineering step: supported pair damping and frozen sources

Declared after the three AMOEBA framework preparations, before any Tinker
calculation. The active goal's standing authorization applies. Preserve the
baseline, existing environments and immutable results. This step does not
change the failed prior hybrid's definition or reclassify its results.

## Question and finite scope

Can a maintained Tinker implementation supply the published pair-damping
capability and freeze QM-source induced variables while preserving justified
damping, full-system boundary atoms and compatible GK accounting?

1. Inspect/pin maintained Tinker/Tinker9 sources and release/build requirements;
   inventory existing local installations without changing them. Verify actual
   POLPAIR, POLARIZABLE, multipole axes, covalent scaling, GK and no-cutoff
   behavior, including which combinations are supported by a GPU backend.
2. Retrieve the primary La parameter supplement and exact parameter files,
   if accessible. Distinguish AMOEBA09, AMOEBA2018 and AMOEBA22. Verify compatible
   Ca coverage and applicability to source damping. Do not infer unpublished
   parameters from biological classifications or inaccessible prose.
3. If needed and supported, build one pinned CPU/reference implementation in
   a new isolated workspace, recording compiler/dependencies/build costs.
   Do not install over an existing environment or change shared permissions.
   A GPU build is deferred until the backend combinations are understood.
4. At most three topology/parameter exports, using the exact already frozen
   GGR1GLG/alpha1F6S/alpha6IP9 physical inventories and successful AMOEBA
   mappings. Check identities, existing hydrogens, disulfides, waters, axes,
   all force-field terms and units. Parameterization/format verification only;
   no source geometry, protonation, cofactor or charge-state changes.

No new scientific energies, forces, minimizations, trajectories, DFT, MACE,
charge fits or training are part of this capability step. Source-level tests,
compilation and format parsing do not establish predictive success. Do not
invoke an energy-analysis program as an alleged format-only check.

## Required decision

Resolve the concrete gaps in AMOEBA_ACCOUNTING.md, or identify their exact
unsupported status. Merely finding the keyword is insufficient. The current
Tinker `kpolar.f` inspected through the primary repository does parse POLPAIR
and a list of polarizable atoms; its actual downstream damping/response and GK
semantics still need verification. No Tinker executable was found on PATH in
the preceding check; this is not an exhaustive installation search.

Before any scoring experiment, declare one complete model, exact physical
preparations, source and cavity parameters, core subtraction, expected call
inventory, numerical/field/partition gates and actual resource execution
settings. Reuse the existing runner and allocation policy. No automatic
benchmark expansion, fitted decision threshold or production promotion.

If the published direct AMOEBA ion model is coherent but the frozen-QM hybrid
is not, document that distinction and declare a separate experiment before
evaluating it. A direct AMOEBA score alone does not complete the MACE goal.
