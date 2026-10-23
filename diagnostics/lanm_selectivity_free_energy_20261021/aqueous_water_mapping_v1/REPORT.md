# Aqueous frozen-water BIO mapping prepared

**Prepared only; no molecular calls or production qualification.** The exact passing17605 six-water candidate is mapped into the unchanged BIO namespace. The previous17612 record supplies the exact selected aquo8 input/native authorities; no donor-fit17621 result was read.

## Concrete outputs

- Two25atom BIO aquo8 controls at the identical SI coordinates, with unchanged bonds and full index/type map. Each pins the original selected SI analyze/analytic outputs; no SI rerun is needed.
- Two11958atom neutral aqueous boxes retained byte for byte:3980waters,Na7/Cl10,49.324Å cube. Periodic numerical policy and dynamics remain unbound.
- Every copied `amoebabio09.prm` is byte-identical to the original `eda92c36` file. Keys add only intrinsic metal definitions900/901 and water-O VDWPR(class73)/POLPAIR(type402). No nonmetal or donor parameter is overridden.

## Exact pair coverage

| Metal | Recipient type/class | Rmin Å | epsilon kcal/mol | Thole | Policy |
|---|---|---:|---:|---:|---|
| La | 402/73 | 3.90613671649 | 0.207651909974 | 0.287823275216 | Frozen17605 override |
| La | 403/74 | 3.52219283708 | 0.0430617273533 | 0.25 | Native mixing; matches SI water H |
| La | 405/76 | 3.57033382368 | 0.468767032351 | 0.25 | Native mixing; unqualified salt/ion transfer |
| La | 414/85 | 4.02497204769 | 0.53035279646 | 0.25 | Native mixing; unqualified salt/ion transfer |
| Dy | 402/73 | 3.58659094374 | 0.317324805668 | 0.283756426912 | Frozen17605 override |
| Dy | 403/74 | 2.99814971034 | 0.0509001790932 | 0.139 | Native mixing; matches SI water H |
| Dy | 405/76 | 3.10470303823 | 0.866927625701 | 0.139 | Native mixing; unqualified salt/ion transfer |
| Dy | 414/85 | 3.78121139259 | 1.02736435749 | 0.139 | Native mixing; unqualified salt/ion transfer |

Intrinsic metal charge+3, zero higher multipoles, polarizabilities, intrinsic Thole and vdW values are preserved from the selected source keys. Both mapped masses are138.91Da under the previously declared constant kinetic-mass convention; source Dy162.5Da is retained in metadata. This mass convention has no static E/F effect.

## Physical difference retained

The sole normalized water-record difference is the existing Urey–Bradley equilibrium H–H length: SI1.5337Å versus BIO1.5537Å, both with−7.60 force constant. Global bonded/nonbonded settings and the other corresponding water records match. BIO1.5537 remains unchanged. Root will separately declare and execute the small-cell numerical comparison with analytic Urey energy/gradient accounting; this preparation does not compute that correction or claim total-E/F equivalence.

Water-H interactions retain source native mixing; Na/Cl pairs use explicit recorded CUBIC-MEAN/HHG/minimum-Thole mixing. Their existence is not evidence of accurate ion association. The boxes are unequilibrated starting configurations, not hydration ensembles. Bulk hydration, periodic electrostatics, polarization under periodic conditions, salt association, an alchemical path, protein donor mapping and sampling remain unqualified.

## Validation and cost

Metadata preparation: 1.709s wall/1.620s CPU, one thread. 83 metadata/negative checks pass (0.367s wall/0.335s CPU), including namespace collision, donor/salt overrides, water-H/O confusion and altered BIO polarizability/Urey records. Zero native, energy, force or dynamics calls.

[Manifest](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/aqueous_water_mapping_v1/prepared_v1/MANIFEST.json) (`0018cc60`); [controls](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/aqueous_water_mapping_v1/metadata_controls_v1/RESULT.json) (`1798d0d1`); [ready binding](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/aqueous_water_mapping_v1/READY.json). Sources and actual inputs are pinned there.
