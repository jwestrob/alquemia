# Hydrated donor transfer: current candidate needs physical qualification

This review reads existing results and source only. Job 17335 completed all 20 AMOEBA rows: 18 new hydrated donor geometries and two reused aquo controls. The root comparator has all eight declared PBE0 references despite its `partial_v1` directory name. No calculation was repeated for this review.

The subsequent root-requested [saved force attribution](WATER_FORCE_ATTRIBUTION_REPORT.md) distinguishes shared internal-force mismatch from metal-dependent response: ligand internal distortion dominates the individual errors, while metal plus ligand dominate the matched Dy-minus-La discrepancy. This qualification prevents attributing the large individual RMSE solely to ion transfer.

## Reported matched observations

Central Dy-minus-La donor-replacement contrasts relative to the common eight-water control, kcal/mol:

| Donor | Frozen-f PBE0 | MACE | Transferred AMOEBA |
| --- | ---: | ---: | ---: |
| Monodentate acetate | +0.575053 | −1.737488 | −6.235815 |
| Bidentate acetate | −1.974185 | −4.566071 | −9.811987 |
| Acetamide | +1.774752 | +1.206225 | +9.162894 |

AMOEBA's individual Cartesian-component force RMSE against these references is 8.992–14.095 kcal/mol/Å; its matched Dy-minus-La component RMSE is 1.776–3.920. The corresponding MACE ranges are 1.529–1.865 and 0.522–0.693. These are **component RMSE**, not atom-vector RMS. MACE is closer on those aggregate metrics, but not every observable: the bidentate paired ligand-direction derivative error is +0.827660 for AMOEBA versus +4.899984 for MACE. Acetamide's AMOEBA error on that derivative is −28.775020.

The balanced contrasts cancel the common free fragments; absolute per-metal replacement energies remain unavailable because no new free-water evaluation was made. Fixed clusters and radial displacements do not establish binding free energies or coordination populations. The comparator explicitly retains `reference_accuracy_qualified=false` and `binding_free_energy_claimed=false`.

## What the implementation identifies

At matched coordinates both ions carry the same spherical permanent charge +3. Saved bonded and permanent-multipole energies are identical between metals. Their energy differences therefore arise from the implemented **van der Waals and induced-polarization terms**. For example, central hydrated acetamide has La/Dy polarization energies −148.468578/−139.858309 and vdW energies +88.686024/+40.781398 kcal/mol. The aquo reference also has a substantial metal-dependent vdW and polarization difference. Aquo subtraction removes that reference contribution, but cannot cancel the nonlinear response of the mixed ligand/water cluster.

This behavior is consistent with a material transfer assumption, rather than a demonstrated atom-number error. The La model uses donor-specific `POLPAIR` and `VDWPR` terms. The Dy asset has no corresponding donor fits and originated with different GEM-water multipoles; this pilot transferred it into the La SI water/ligand environment. La and Dy also differ in polarizability, damping and buffered-14-7 parameters. The native source combines ordinary Thole parameters by their minimum (except a zero-value special case), then applies explicit `POLPAIR` overrides (`kpolar.f`, lines 493–518). Its damping depends jointly on the pair parameter, distance and polarizabilities (`damping.f`, lines 98–117). A changed donor environment can therefore alter the paired response even when all mappings are exact.

The large individual force errors can include shared water/ligand internal-force differences. The smaller paired errors remove much of that common contribution but do not uniquely assign the remaining error to a particular metal radius, damping constant, donor multipole or missing physical effect. The successful narrow derivative controls and reproduction of published La dimer energies establish implementation consistency, not transferable physical accuracy. The known La monodentate-acetate underbinding remains relevant.

## Concrete consequence

Do not promote this transferred pair to protein production sampling or copy La corrections onto Dy to fit these inspected contrasts. The prepared Asp/Glu/Asn/backbone fragments now provide actual protein donor parameters and exact coordinates for a finite, separately declared paired qualification. That test should retain both metals and compare energy differences and forces under the same geometry/state convention; small-molecule transfer and native protein donor behavior must remain distinct. A reference with qualified electronic-state and derivative behavior is still needed before treating its residuals as force-field fitting targets. Any subsequent parameter development needs a declared training set and separate validation fixtures.

There is no source-evidence basis here for broad MD or a unique parameter repair. The immediate question is whether one coherent paired ion/water/protein-donor Hamiltonian reproduces the relevant local response, including mixed-shell polarization.

## Evidence pins

- Root [comparison result](../../../workspaces/lanm_selectivity_free_energy_20261021/hydrated_comparison_v1/partial_v1/RESULT.json): `e096e732d70a7e352466ab8a1dc46ff19de8427416a0f21e6dcce3a960352c7c`.
- Native [hydrated AMOEBA result](../../../workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/hydrated_mixed_run_v1/RESULT.json): `2d55465328abcd60dcf0ca5b14761d3f4c96314310c7ec0381a2696d8c3f8674`.
- [Fixed-panel report](PANEL_REPORT.md) preserves exact La/Dy asset provenance, mapped parameters and numerical controls.
- [Protein fragment report](PEPTIDE_FRAGMENT_REPORT.md) records the eight prepared configurations and zero-energy boundary.
- Native Tinker source commit: `87050685eff8840d312e2a332cc82c33f63c7c3d`.
