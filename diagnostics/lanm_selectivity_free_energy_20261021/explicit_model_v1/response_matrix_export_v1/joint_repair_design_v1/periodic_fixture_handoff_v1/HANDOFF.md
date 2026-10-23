# Existing aqueous boxes: static implementation handoff

The two existing boxes supply concrete, unchanged geometries and topology. **They do not yet supply a qualified periodic Hamiltonian for the custom response.** No key, coordinate, parameter or engine source was changed in this handoff.

## Exact fixtures and prior admission

Both files under `workspaces/lanm_selectivity_free_energy_20261021/aqueous_water_mapping_v1/prepared_v1/{La,Dy}_aqueous_BIO/` contain 11,958 atoms: one M3+, 3,980 intact waters, seven Na+ and ten Cl−. The exact monopole sum is zero. Both use the 49.324 Å orthogonal cube and identical coordinates, bonds, water/ion order and nonmetal types; only the metal name/type differs. Every water has reciprocal O–H bonds and no metal bond.

| Artifact | SHA256 prefix |
|---|---|
| La `system.xyz` | `b59954d7` |
| Dy `system.xyz` | `1b7c52d2` |
| Both mapped `system.key` | `c2f77a0a` |
| Unchanged `amoebabio09.prm` | `eda92c36` |
| Mapping manifest | `0018cc60` |

The XYZ files are byte-identical to `paired_legs_v1/prepared_v1/{La,Dy}_aqueous.xyz`. Full paths and hashes, original preparation and all relevant admissions are in `INPUTS.json` and `INVENTORY.json`.

The only actual BIO numerical admission is job17668: two nonperiodic 25-atom aquo8 controls, four native programs. They use the older frozen17605 water-pair candidate. The known Urey energy/gradient difference was independently reproduced; remaining energy and gradient differences were about 1.81e−9 kcal/mol and 1.48e−9 kcal/mol/Å. Both full boxes remain prepared-only, unequilibrated and periodically unqualified. No new physical/reference outcome was read.

## Exact mapping required; no mapping applied

| Role/operator | Qualified custom vacuum namespace | Existing BIO namespace |
|---|---|---|
| La type/class | 500/500 | 900/900 |
| Dy type/class | 501/501 | 901/901 |
| Water O type/class | 36/34 | 402/73 |
| Water H type/class | 37/35 | 403/74 |
| Metal–O VDWPR | class 34 with 500/501 | class 73 with 900/901 |
| Metal–O POLPAIR | type 36 with 500/501 | type 402 with 900/901 |
| Response record and graph | metal 500/501, O36, H37 | would need explicit metal900/901, O402, H403 support |
| Sodium / chloride | absent from qualified vacuum fixture scope | types405/414, classes76/85 |

The full pinned BIO library assigns types1–416 and classes1–87. **500/501 themselves do not collide**, even in this full library; neither do reserved900/901. Retaining500/501 could therefore be an explicit alternative namespace decision, but it would require a separately declared rewrite of the metal atom and every corresponding reference. Nothing was rewritten here. This namespace observation is not an extension of the current supported topology.

Copying the custom SI key wholesale is unsafe: BIO36/37 are Ser OG/HG; BIO401 is a phosphate oxygen; BIO402/403 are water O/H; BIO404 is Li+. Custom donor401–404 and amide104/105 would overwrite unrelated BIO chemistry. Small-system class numbers also cannot be substituted for BIO class numbers. The exact role records and full-library namespace audit are retained in `NAMESPACE_AND_REPRESENTATIVES.json`.

The existing mapped keys are not a current-candidate binding. Their old17605 metal–O values differ from the admitted17676 baseline used by the custom-response qualification:

| Metal | Existing BIO R / epsilon / POLPAIR | Custom qualification baseline R / epsilon / POLPAIR |
|---|---|---|
| La | 3.9061367164924641 / 0.20765190997435956 / 0.28782327521600637 | 3.9193948492218289 / 0.15863887317162664 / 0.17662704152555941 |
| Dy | 3.5865909437359376 / 0.31732480566797089 / 0.28375642691227121 | 3.6401663222400971 / 0.20247982331452399 / 0.15457541153280463 |

No live fit outcome or later coefficient selection is substituted. A future composition must bind the intended exact candidate explicitly. Water alpha remains O0.837/H0.496 Å³, Na0.120 and Cl4.000; intrinsic metal alpha is La1.6597/Dy0.728. BIO keys use common metal mass138.91 Da, while the custom Dy source mass is162.5 Da. That mass difference does not affect a static potential but must remain explicit before dynamics.

## Urey difference and cancellation scope

After role renumbering, the compared water multipoles, local-frame definitions, alpha/Thole, vdW, bond and angle records match. The existing water-model difference is the internal H–H Urey term: k=−7.60, rSI=1.5337 Å and rBIO=1.5537 Å, with native unit1 and no cubic/quartic terms. Its difference per water is

`delta U = k [(rHH − rBIO)^2 − (rHH − rSI)^2]`.

At identical geometry it changes internal hydrogen gradients by equal-and-opposite central forces; it does not change water intermolecular multipoles or polarization parameters. It cancels in La-minus-Dy energies when the same water internal geometries occur, and in balanced interaction cycles that subtract the same water geometries. It also contributes zero work to a rigid-water translation or rotation because rHH is unchanged. These cancellation statements do **not** establish full Cartesian-force parity: the individual H gradients differ. Flexible relaxation/dynamics can alter water internal distributions and remove the cancellation; the BIO value must not be overwritten to force equality.

## Periodic settings: actual versus previously proposed

The actual mapped keys contain `polar-eps 1e-8`, reporting precision and one-thread settings, plus metal definitions/pairs. They contain **no Ewald enablement, grid, cutoff, surface or periodic iteration policy**. The cell is present in XYZ; that alone does not bind a PME method.

`paired_legs_v1/REPORT_v2.md` proposed, but did not execute or bind: conducting/tinfoil, Ewald real cutoff9 Å, permanent/polarization alpha0.4 Å−1, 64³ grid, both orders5; vdW cutoff9 Å/taper8.1 Å and no long-range vdW correction; mutual tolerance1e−7 and100 iterations. The latter tolerance differs from the actual mapped/current-qualified1e−8 setting. These remain prospective choices requiring an explicit later declaration and native readback.

The pinned current `kewald.f` starts with TINFOIL but interprets any `EWALD-BOUNDARY` keyword as VACUUM. Do not encode `EWALD-BOUNDARY TINFOIL`. All PME, cutoff and surface settings must be bound and read back explicitly in a future static test.

## Exact cell and literal-representative bounds for the tail definition

For both boxes, H=49.324 I Å, so the lattice singular-value lower and upper bounds are exactly **49.324 Å**. Volume is119998.237484224 Å³. For every literal water-O position relative to the literal metal:

- Exact maximum squared displacement: **1891.601855145683 Å²**.
- Upward-rounded bound: **D=43.492549421087 Å**, verified by squaring the decimal upper value.
- This is a fixed-fixture representative bound, not a bound on future trajectories.

Do not substitute the nearest-image half-diagonal42.715837016264 Å: it is too small for these literal representatives. O coordinates span49.335848 Å in z, slightly greater than L, and do not all share one axis-aligned fundamental cube. Some atoms extend further because waters were retained whole. No wrapping or pair-specific lattice reindexing was performed. The generic canonical-cell diameter85.431674032528 Å is a separately labeled geometrical value; the actual D above is justified directly from all3980 literal displacements, without assuming a canonical-cell representative convention. No tail number was evaluated.

## Concrete next static implementation inputs and missing support

Use these two exact geometry/topology artifacts as prospective static fixtures, with a future explicit whole-model composition decision for namespace, candidate pair/response coefficients and BIO intramolecular terms. Then bind periodic numerical settings and salt treatment. The current admitted refactor still rejects periodic/Ewald/list routes, n>28 and all non-water/non-metal types; it accepts response records only for500/501–36 and expects H37. It therefore cannot execute these boxes unchanged. Source support for mapped water/metal selectors, unchanged Na/Cl background interactions, a declared periodic image convention and an adequate stability check must precede any static evaluation. Relaxing guards alone is not implementation qualification.

Preparation here used one metadata thread, 0.746 CPU-seconds for the main inventory, zero native calls, zero geometry/parameter changes and zero tail evaluations. The supplementary metadata record preserves the rejected centered-cell assumption and a prior prose-tokenization failure; neither changed scientific inputs.
