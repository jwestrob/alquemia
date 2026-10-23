# Fixed-geometry electrical response: CPU feasibility

## Decision

**Supported by the installed CPU engine, without a build or Hamiltonian change.** A copied, exact selected input can use `EXTERNAL-FIELD Fx Fy Fz` and the installed `analyze` executable. The field enters the mutual induced-dipole solve, and native total energy includes both permanent and induced field coupling. This supports a prospective fixed-geometry static polarizability comparison against the same PBE0 electronic state. It does not establish that the model response is too weak or too strong.

**Cheaper preceding observation:** `analyze system.xyz EM` can report the full zero-field permanent-plus-induced molecular dipole without changing the Hamiltonian. Existing QM dipoles can be collected first. Neither a zero-field dipole nor a population partition is a polarizability measurement. Symmetry can make a total dipole small despite substantial redistribution.

This is a source/metadata feasibility record: zero molecular calls, builds, submissions or engine edits. No field magnitude, new fit, parameter sweep, or execution matrix is released. Source/input/binary pins are in `SOURCE_PINS.json`. The selected source remains failed-science 17676, not either positive EXCHPOL perturbation.

## Exact installed route

Source paths below are relative to `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/tinker_source/source/`.

| Requirement | Source evidence |
|---|---|
| Static uniform field and units | `prmkey.f:402` reads three `EXTERNAL-FIELD` components and divides by `elefield`; `EXFLD-FREQ` is a distinct, excluded option. Input units are MV/cm. |
| Self-consistent response | `induce.f:204–248` adds the field to **both** permanent driving fields before constructing direct dipoles and performing the ordinary mutual solve. |
| Permanent field energy | `exfield.f:35–99` implements `−K Σ(q r + μperm)·f`; `empole.f`, `empole1.f`, and `empole3.f` call the energy, gradient, and analysis versions respectively. A uniform field has no direct quadrupole coupling. |
| Induced field energy in the proposed executable | `analyze` calls `analysis` and `epolar3`; the nonperiodic `epolar3a` calls `induce` and explicitly adds `−K/2 Σ μind·f` at `epolar3.f:553–565`, alongside its usual pair polarization energy. Do not append another external-induced energy correction. |
| Other native energy/gradient routes | `epolar.f:1797` and `epolar1.f:8445` use `−K/2 Σ μind·udirp/α`; `udirp` already contains the external field. The ordinary gradient route also calls `exfield1`. This establishes source support, not a newly executed finite-field force qualification. |
| Zero-field total moment | `moments.f:101–218` rotates permanent multipoles, calls `induce`, adds induced dipoles, and sums charges plus dipoles. `analyze.f:766–774` prints charge and the three total dipole components. |
| Higher-precision induced moments | `SAVE-UINDUCE`, parsed by `optinit.f:56`, makes `analyze.f:995–1024` save each polarizable site's induced components in Debye. This is a **key option**, not an `analyze U` mode. |

The installed `analyze` symbol table contains `prmkey_`, `induce_`, `exfield3_`, `epolar3a_`, `moments_`, and `saveyze_`. Its hash matches the pinned, previously used CPU executable. No executable was invoked for this review. The public [Tinker User's Guide](https://dasher.wustl.edu/tinker/downloads/tinker-guide.pdf) also documents the field in MV/cm; the pinned local implementation governs these conclusions.

The standalone `polarize.f` source is not the minimal route: no corresponding binary is installed in this build; it uses a hard-coded internal field and a duplicated `moluind` solver rather than the production `induce` path. In particular, that duplicate does not apply the EXCHPOL route. No build is proposed.

## Units, charge and precision

Let `F` be the field in MV/cm, `L = elefield = 1439.96455`, `f = F/L` the internal field, and `K = electric/dielec`. The admitted water inputs use dielectric 1 and `electric = 332.0637133`. Native induced dipoles before output are in e Å; `D = 4.803206802` converts them to Debye. Native energies are kcal/mol and site polarizabilities are Å³.

For a single axis, with equal positive and negative field magnitude `F`:

```
alpha_jj [Å³] = -(L²/K) * (E(+F ej) + E(-F ej) - 2 E0) / F²
alpha_ij [Å³] = (L/D) * (sum u_i(+F ej) - sum u_i(-F ej)) / (2 F)
```

The second expression uses saved induced components in Debye. At fixed geometry, unchanged permanent moments cancel from this difference. Six signed axial fields provide all tensor columns through induced moments, plus the three energy-curvature diagonals. Six energy values alone provide only the diagonal/trace, not off-diagonal elements.

For the quantum comparison, use a single explicit conversion convention: one polarizability atomic unit is `bohr³ = 0.1481847111705721 Å³` with the pinned `bohr = 0.529177210544 Å`. An exactly native-energy-consistent conversion from a quantum field `h` in atomic units is `F = h * hartree * L / (K * bohr)` at dielectric 1, approximately `5142.2067588 h MV/cm`, with `hartree = 627.509474063`. The conventional `L/bohr²` conversion differs here by only about 7e−10 relatively because native constants are rounded; record the convention instead of mixing constants silently.

For these water/ion graphs, the direct and energy permanent-field scaling agree: the complete water polarization group has zero intramolecular direct and 1–2/1–3 permanent polarization scales; distinct waters and the unbonded metal use unity. Thus the ordinary induced response is consistent with the field-energy derivative for this scoped model. Do not generalize this check to arbitrary protein scaling conventions.

All clusters have charge +3. At a fixed laboratory origin, a translation by `a` changes the dipole by `Q a` and changes field energy by `−K Q a·f`. It changes the odd linear term, **not the polarizability curvature**. Keep exact coordinates and one origin across signs, methods and fragments; use no periodic boundary or independent recentering. A uniform field on a free charged molecule is a fixed-nuclei electronic probe, not an equilibrium sampling setup.

There is a separate **printed-moment origin trap**: `moments.f` subtracts each structure's center of mass. Convert its printed moment to the shared laboratory origin as

```
mu_lab [Debye] = mu_printed_COM + D * Q * COM[Å]
```

using the actual native masses, active-atom set and coordinate frame. Keeping the same parent coordinates alone does not make different fragment COM moments comparable. A common-origin charge-neutral cycle is then origin invariant.

Precision is fixed by the implementation: energy is printed to 8 decimals in kcal/mol even with `digits 10`; total dipole components to 3 decimals in Debye; saved induced components to 6 decimals. The field magnitude must be declared prospectively with a rounding and induction-error budget. For `n` saved sites, each summed induced component has a worst-case print-rounding error `n × 0.5e−6 Debye`; the signed difference has twice that bound. The energy numerator has at most `2e−8 kcal/mol` print-rounding error before solver error. The existing induction residual limit is 1e−8, but it is not by itself a guaranteed molecular polarizability error bound. No field magnitude is chosen here.

## Minimal faithful adapter, if root declares it

1. Bind each exact selected XYZ/key/prm and the existing raw zero-field energy. Copy them to a fresh private cell directory. Reject existing external/time-varying fields, periodic/solvent settings, partial active-atom selections, or unrequested EXCHPOL terms. Preserve all molecular parameters and coordinates. The selected compact DATA's `baseline_model[gid].provenance.effective_input` is the authority, not the older training input.
2. Append exactly one static field record plus reporting-only `polar-print` and `save-uinduce`. Run the pinned single-thread `analyze system.xyz EM`; no testgrad, dynamics or geometry relaxation is needed to obtain the electrical response. Retain energy, total moment, all site dipoles, raw output, field/units, input pins, residuals and process receipt.
3. Use the existing strict native energy parser and all-summary finite induction admission. Require every expected site/component and a complete declared denominator. Retain the field-energy and dipole derivative estimates separately; require finite tensors, symmetry within a predeclared numerical budget, zero-field energy replay, and energy/moment sign-and-unit consistency. A source-only isolated-site check predicts `muind = alpha f` and `Delta E = −K alpha |f|²/2`; it is an appropriate prospective sign/factor control, not an executed result.
4. The smallest balanced **cluster** comparison would be the four already consumed W8/W9 × La/Dy structures at six signed axial fields: 24 new model analyze calls and 24 matching quantum field single points, with exact zero-field energies reused. Root must freeze the actual matrix, field size and reference field convention first. This count excludes any separately declared zero-field moment extraction, fragment response, or numerical control calls; none is automatically added here.

The fixed-geometry baseline induced-dipole equations are linear in the applied field: their response matrix depends on fixed geometry, polarizabilities, damping and mutual interaction, not on the permanent driving field. A field response can therefore constrain that matrix without another parameter fit. Agreement would not establish correct permanent local fields or the full coordination energy; disagreement would identify a directional electrical-response mismatch, not a unique missing term.

## Saved zero-field redistribution can precede finite fields

The completed common-basis result reported by root has reference `B = +9.036179` versus model `+2.797568 kcal/mol`; basis equalization strengthens the discrepancy. Positive EXCHPOL and alpha probes did not rescue it. These observations motivate an independent electrical observable, but do not identify the remainder uniquely as polarization.

For each parent metal/shell geometry, the matched dipole counterpart of the energy remainder is

```
mu_B = mu_full - sum_i mu_dimer_i + (N-1) mu_ion
       - mu_shell + sum_i mu_water_i .
```

Its coefficients cancel charge, provided every value is expressed at the same origin. Use the same full-cluster basis and exact fragment masks for the existing QM quantities; retain raw electronic/nuclear/total dipoles and native units, and verify the native origin convention before combination. Mulliken/Loewdin populations are supplementary, basis/partition-dependent evidence. Do not invent Hirshfeld data or interpret populations as uniquely observable charge transfer.

For the classical counterpart, unchanged permanent atom contributions cancel from this exact matched cycle. After verifying atom mapping, types, geometry and permanent multipoles, `mu_B` can be computed from **induced-dipole sums alone**, avoiding fragment COM correction. The isolated model ion has zero induced dipole at zero field. Existing analyze-E outputs do not automatically contain the necessary full dipoles: collecting missing model moments would require new, separately declared model calls, even though the Hamiltonian is unchanged. No such calls were made here.

Zero-field `mu_B` probes redistribution in the current internal field; polarizability probes its derivative under an external field. They are complementary, and symmetry can suppress either total-vector signal. Neither establishes a global family limitation, donor accuracy, a production model, or a bulk-water/periodic qualification.
