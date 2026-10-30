# Independent method review: bounded La/Dy ALMO quartet

2026-10-23. Scope: supplied marching orders, official method documentation, and read-only inspection/hash checks of archived inputs and results. No electronic calculation, numerical replay, scientific-runtime import, job, or alternative method was performed or constructed.

**Conclusion:** the requested four-state comparison is scientifically interpretable if the exact Hamiltonian, fragment states, ALMO polarization constraint and energy bookkeeping are preserved. **Installed executable/version and licensed access are unverified in this review, pending Singer's software-access record.** Documentation establishes a possible implementation route, not local availability or an admitted calculation.

## Exact archived inputs and states

Let `W` denote `workspaces/lanm_selectivity_free_energy_20261021/`. The original executed aquo directories are `W/reference_chemistry_v1/ecp_remaining_v1/tasks/ecp_pbe0__aquo8_central_{La,Dy}/execution/`; mixed-shell directories are `W/hydrated_chemistry_run_v1/reference/tasks/ecp_pbe0__hydrated_acetate_mono_central_{La,Dy}/execution/`.

| Case | Atoms | Total charge / effective multiplicity | Effective electrons | SHA256 of original `geometry.xyz` |
|---|---:|---|---:|---|
| aquo8_central_La | 25 | +3 / 1 | 88 | `371533e94b5e7b2bbe61f1ecec77d1ddbfdbc4edd0beea98af77dad23bf39c4a` |
| aquo8_central_Dy | 25 | +3 / 1 | 88 | `8448c74ce91cb7bfdd814566a17f1bc0e219a8f8b568747cc488fb5756781c62` |
| hydrated_acetate_mono_central_La | 29 | +2 / 1 | 110 | `2501a764ae1b252eaeba4df5c26235335f804de113d0d526ec594995957a7228` |
| hydrated_acetate_mono_central_Dy | 29 | +2 / 1 | 110 | `89f08e66242fcd1798ee89d3ed9bfba01f874066afc25f16083b513a47c58fdf` |

Independent literal-token comparison confirms all nonmetal element/coordinate rows match between La and Dy within each shell family. Preserve original coordinates and atom order; no optimization, regeneration, or replacement by common-basis fragment geometries.

Original `endpoint.inp` hashes, in table order:

```text
da6264c6c74763d3ac0bd22d8cb1ec0f0744c56c4d1051a0534ef3d2e3e66b26
92a6f3264039af1eb4fc18269dbaa919a8f05625d1d9457b5ada8d328658963c
5675a433485b41b78f258e8fb85dcf6d5fe133e7ac0cd6c0b7cdd45646b20af9
570cba44e93ea26bfcd6610fb2568f621dd9125d80633020e78e2f77e5c77163
```

All four decks specify `PBE0 D4 def2-TZVP def2/J RIJCOSX NoAutostart DefGrid3 VeryTightSCF EnGrad`, `Guess PModel`, `MaxIter 300`, and explicit metal orbital/auxiliary basis and ECP blocks. La removes 46 electrons; Dy removes 55, retaining the archived frozen-f effective singlet. The exact metal assets also occur at `W/reference_chemistry_v1/common_basis_decomposition_v1/prepared_v1/assets/La_basis.inc` (`c5298b309b2d897e5de59ddb6b4c7992913dd8cba18b7cd21f9757ce09b167d3`) and `Dy_basis.inc` (`87bb209a42eca365c93d8f69088a57371285623a1d4afa1c1336390279fafb20`). A named default metal basis/ECP is not an interchangeable substitute.

Fragment 1 is M3+, eight explicit electrons, effective singlet. Fragment 2 is the **entire** shell: aquo shell neutral/80 electrons/singlet, or acetate-minus plus seven waters, charge -1/102 electrons/singlet. Physical Dy multiplicity six is not substituted for the archived effective singlet. These calculations cannot resolve omitted explicit-f or spin-orbit effects.

## Documented implementation requirements

Q-Chem documents general orbital bases and user ECPs, including local/nonlocal channel order, removed-core count and radial-power convention. Translation must preserve the actual coefficients and angular conventions; accepting a similarly named ECP is insufficient. This is documented capability, not proof of compatibility of every installed ALMO/version/custom-ECP combination. [Q-Chem 6.4 user-defined ECPs](https://manual.q-chem.com/6.4/subsec_ECP-userdef.html)

The ALMO definition must be explicit. Q-Chem 6.4 distinguishes AO-block polarization (`EDA2=2`) from response-function polarization (`EDA2=1`, documented there as nMDQ-uFERF); merely requesting EDA does not fix an implementation-independent definition. Preserve a single declared choice across the quartet; do not silently switch definitions to obtain convergence. Both fragments must be allowed intrafragment relaxation, while metal-shell delocalization remains excluded in the polarized state. [Q-Chem 6.4 EDA2 controls](https://manual.q-chem.com/6.4/sect_eda2_job.html)

The AO-block polarization/CT partition has a known large/diffuse-basis limitation; response-function constraints define a different polarization subspace. Neither label supplies an experimental partition of the energy. [Q-Chem polarization definition](https://manual.q-chem.com/6.4/sect_ferf.html)

The pilot must establish the requested full-energy agreement with the archived protocol, with electronic and D4 contributions identified. RIJCOSX, grids and auxiliary bases are numerical implementation details that cannot be assumed identical across programs. The acceptance tolerance must be declared before the pilot; no tolerance or substitute protocol is selected by this review. Preserve fragment, polarized and full SCF residuals and any failure. Small printed SCF changes alone are not a rigorous error bar.

## Closed comparisons

At each fixed geometry, retain absolute energies `E_F`, `E_P`, `E_U` (frozen, polarized, full) and `E_R = E_M3+ + E_shell` using the declared fragment basis/reference convention. Define:

```text
FRZ = E_F - E_R
POL = E_P - E_F
CT  = E_U - E_P
E_U - E_R = FRZ + POL + CT
DD(X) = X_Dy,acetate - X_La,acetate - X_Dy,aquo + X_La,aquo
DD(E_U) - DD(E_P) = DD(CT)
DD(E_U) = DD(E_R) + DD(FRZ) + DD(POL) + DD(CT)
```

These identities distinguish absolute CT-suppressed energies from interaction-energy increments. The standard ALMO stages and fragment-reference construction are documented in the [ALMO-EDA manual](https://manual.q-chem.com/latest/sect_eda1.html). For this matched quartet, isolated metal terms cancel in DD, and shell terms should cancel when their geometry, electronic state and fragment basis truly match. Supersystem/ghost-basis conventions can spoil the latter cancellation; report the actual reference DD rather than assume zero. Fragment energies at complex geometries exclude deformation from separately relaxed shells. No additional relaxed reference calculation is requested.

The whole-shell partition retains internal water-water and acetate-water delocalization. Thus “CT suppressed” means **metal-shell** CT suppressed, not a CT-free shell or a ligand-by-ligand CT decomposition.

## D4, BSSE and numerical closure

- Archive raw electronic energies and empirical D4 separately. An EDA frozen-energy dispersion label is not automatically the archived D4 correction. Identical per-geometry D4 added to all three stages cancels in `E_U-E_P`; differing fragment/full or state-dependent corrections require explicit accounting. The full comparison must retain the archived PBE0-D4 total once, not omit or double-count dispersion.
- Keep raw and counterpoise-corrected quantities distinct. Q-Chem exposes fragment super-basis and BSSE options separately. Do not sum corrected CT against an uncorrected total or introduce a new BSSE convention merely because a default enables it. Legacy ALMO documentation also records a general-basis restriction for its BSSE route; compatibility must be checked for the exact chosen driver/version rather than inferred globally. [EDA2 controls](https://manual.q-chem.com/6.4/sect_eda2_job.html), [legacy ALMO BSSE notes](https://manual.q-chem.com/latest/sect_eda1.html)
- Report both per-cell component closure and quartet DD closure, full-energy replay residuals, and all available numerical-convergence diagnostics. Compare CT DD with numerical uncertainty and the pre-existing model discrepancy, not solely a sign reversal of a tiny total contrast.

## Inference and stop condition

A substantial CT DD would establish sensitivity of this fixed-geometry electronic exchange comparison to releasing this particular metal-shell constraint. Large absolute CT energies can cancel in DD; energetic stabilization is not the same quantity as net charge transfer. [Official CT-analysis description](https://q-chem.com/explore/molecular-interactions/ct-almo/)

It would not uniquely identify the force-field error as missing CT, establish protein selectivity or an aqueous free energy, or validate a new response law. Small CT would likewise leave polarization, frozen interactions, geometry and sampling unresolved.

**Recommendation:** proceed only with the requested documented, licensed, compatible ALMO implementation and unchanged quartet, using the authorized La-aquo pilot before the other three. If unavailable or incompatible, report the precise blocker and stop. SAPT/cDFT, population analysis, ETS-NOCV, frozen-density embedding, altered ECPs/bases/spins, or a hand-built localization approximation answer different questions and must not substitute. This review authorizes no fitting, additional geometries, retries, or peptide campaign.

## Free implementation candidate: CP2K, not yet qualified

Following the user's request for a viable free alternative, CP2K's official documentation and current public source were inspected without installation or execution. **CP2K is a genuine native ALMO candidate, but the exact requested combination is not yet qualified.** These observations concern the inspected current documentation/source; they do not establish any locally installed version or compiled-library capability.

- **Constraint and full states:** native ALMO provides block-diagonal optimization, molecular-fragment initialization, and `DELOCALIZE_METHOD FULL_SCF`. Its default `XALMO_X` is a spatially restricted single-excitation correction and is not the required fully relaxed state. The documented closed-shell restriction is consistent with the archived effective singlets. [CP2K ALMO_SCF](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/ALMO_SCF.html)
- **Charged fragments:** the current `FRAGMENT` interface explicitly accepts charge and multiplicity. Thus +3 metal and 0/−1 shell are not excluded merely by the charge input interface. Effective electron counts and the sum of fragment charges must still match the original full problem. [CP2K FRAGMENT](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/ALMO_SCF/FRAGMENT.html)
- **Disconnected whole shell:** `CONNECTIVITY USER` defines molecular groups from coordinate labels, so all shell atoms can share one group without adding fictitious covalent bonds. A native charged-singlet ALMO regression input demonstrates this grouping route. Automatic distance-based molecule detection would instead risk splitting the water shell and changing the experiment. [Topology documentation](https://manual.cp2k.org/cp2k-8_2-branch/CP2K_INPUT/FORCE_EVAL/SUBSYS/TOPOLOGY/MOL_SET/MOLECULE.html), [official charged-fragment regression input](https://github.com/cp2k/cp2k/blob/master/tests/QS/regtest-almo-1/almo-singlet-charge.inp)
- **Original ECPs:** CP2K is not restricted to GTH potentials. `KIND/ECP_SEMI_LOCAL` explicitly supports original semilocal ECPs when the required integral library is present; disabling it instead constructs a nonlocal one-center expansion. Exact archived coefficients, local/nonlocal channel conventions, basis functions and library support remain compatibility requirements. [CP2K KIND/ECP_SEMI_LOCAL](https://manual.cp2k.org/cp2k-2024_1-branch/CP2K_INPUT/FORCE_EVAL/SUBSYS/KIND.html)
- **PBE0 combination remains unverified:** CP2K supports PBE0/Hartree–Fock exchange generally. The inspected ALMO entry and KS-update source showed no explicit hybrid prohibition, and the ALMO density is passed to the common KS updater. This is supporting source evidence, not proof of correct combined hybrid/ECP/charged-fragment EDA operation. No matched PBE0 quartet or corresponding compatibility test was executed. [CP2K features](https://www.cp2k.org/features), [ALMO KS interface](https://github.com/cp2k/cp2k/blob/master/src/almo_scf_qs.F)

**Concrete unresolved frozen-reference issue:** the native analysis documentation offers `FROZEN_MO_ENERGY_TERM SUBLATTICE` and explicitly warns that its reference term lacks an interpretation for charged fragments. M3+ is charged in every requested cell. The molecular nonperiodic case therefore needs a documented/source-supported resolution of this reference convention before CP2K can be represented as delivering the complete requested frozen/polarized/full EDA. A successful charged-fragment SCF, or a periodic regression input using the option, does not settle that scientific interpretation. This is a blocker to claiming compatibility now, not proof that every nonperiodic charged-fragment implementation is impossible. [CP2K ALMO analysis](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/ALMO_SCF/ANALYSIS.html)

No GTH substitution, PBE-for-PBE0 change, approximate XALMO correction, periodic charged-cell reference, or custom replacement EDA was constructed or recommended as an automatic workaround. Original D4 and full-energy agreement requirements still apply. Installation, builds and calculations remain outside this review.

Input-bundle handoff: `inputs/PROVENANCE.json`, SHA256 `677c759c1045e2d7d5ab2dd05977e4a516ca04ff8d520c40111669eabb452d66`, now records the copied quartet and original evidence. Stored `printed_scf_criteria` arrays are empty; their emptiness is not a claim that all residuals passed and provides no new precision bound.
