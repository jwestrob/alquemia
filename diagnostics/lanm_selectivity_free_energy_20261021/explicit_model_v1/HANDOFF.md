# Explicit polarizable model: implementation checkpoint

21 October 2026. Paused after the parent's documentation-only instruction. No additional calculation or submission followed that instruction.

## What is implemented

An isolated native Tinker 26.2 CPU build can evaluate the exact published La AMOEBA22 donor corrections, including `POLPAIR`. Source commit is `87050685eff8840d312e2a332cc82c33f63c7c3d`. Numerical sources are unchanged; the isolated CMake source list was repaired to include omitted `uatom` and `uniquetyp`. The failed first build remains in `build.log`; successful recovery is in `build_v2.log`. Compilation used portable `-O2`, at most four build processes. OpenMM 8.5 is installed, but its standard per-atom Thole interface does not directly preserve these pair-specific corrections; no OpenMM replacement or GPU model was run.

`PREPARE.py` generated 28 matched physical geometries and four free-ligand fixtures from the chemistry team's exact panel. The common parameter set preserves the complete La SI parameter file and key, including its historical water atom types. It appends the author-contributed Dy ion under a new type/class. No ALPB, implicit solvent, MM/MACE combination, parameter fitting, protein simulation or dynamics has run.

## Completed numerical control

Only the four-atom central La–water fixture was evaluated: one native `analyze` invocation and one native `testgrad` invocation. The latter internally evaluates displaced energies and must not be counted as a single energy evaluation. Both used one numerical thread, mutual polarization with `polar-eps 1e-8`, and the prepared input unchanged.

The native output reports intermolecular energy **−91.27177843 kcal/mol**, reproducing the published SI value **−91.27**. Total energy is −89.03447220 and the polarization term is −50.35329054 kcal/mol. The saved full Cartesian analytic/numeric gradient table agrees to approximately 1e-8 in its printed components at the predeclared 1e-5 Å finite-difference step. This is a narrow implementation control; no independently computed maximum-error statistic or paired-model acceptance has been produced. The declared acceptance threshold remains 1e-3 kcal/mol/Å.

Commands, run from `prepared/water_central_La`, with `OMP_NUM_THREADS=1`, `OMP_THREAD_LIMIT=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and `NUMEXPR_NUM_THREADS=1`:

```text
timeout 45 ../../tinker_build/analyze system.xyz EP > smoke_analyze.out 2>&1
timeout 45 ../../tinker_build/testgrad system.xyz Y Y 0.00001 > smoke_testgrad.out 2>&1
```

The first tool invocation reported about 0.52 seconds wall time. A separate native-process accounting receipt and exact build CPU cost were not recorded; do not infer them from successful output. Raw outputs remain available. The full Slurm runner and automatic result collector are not implemented.

## Constraints for the full plan

- **A mathematically executable paired model is not yet a validated paired parameterization.** The repository Dy ion accompanies different GEM-derived water multipoles. Its transfer to the common La-SI water and donors is a declared hypothesis. Do not copy La-specific donor damping to Dy or describe this mixture as an established published La/Dy protein model.
- Preserve the published monodentate-acetate weakness: the La SI reports −466.67 kcal/mol QM versus −417.58 AMOEBA22 for its listed geometry. Water agreement alone does not qualify donor transfer.
- Complete the declared paired water/donor numerical controls, then compare the finite shared geometry panel with the independent reference-chemistry results. Only after those gates should a protein/water model and periodic sampling implementation be considered.
- The experimental target agent's primary construct is NH3+–FIDTNNDGWIEGDELLA–CONH2, peptide charge −4 and La(III)/Dy(III) complex charge −1. The old free-carboxylate model is a separate alternative. Small aquo/donor fixtures do not resolve protein terminal typing or experimental-source ambiguity.
- Classical spherical +3 ions have no explicit electronic spin state. Paired force-field energies do not establish the QM state convention, hydration free energies, coordination equilibrium, selectivity or experimental affinity.

## Exact checkpoint pins

All workspace paths below are under `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/` unless noted.

| Artifact | SHA256 |
| --- | --- |
| Diagnostic `SPEC.md` | `e0f6137be480d95d248040abd92655884bc4b6d832531f7a22767999b1de2fe0` |
| Diagnostic `PREPARE.py` | `5e23fb26585cace8558b181c2ef6efccd17d2f00edd6c7514e8802d6828c0071` |
| `prepared/MANIFEST.json` | `0ff77d11c4c1c6773802ea1767260d0d2c177f2adeec1d7c491ed31f6b414813` |
| `prepared/water_central_La/smoke_analyze.out` | `12bb12056d4110bc669782eb31574070444cdab0ebfecc644fe4da45631428a9` |
| `prepared/water_central_La/smoke_testgrad.out` | `854e4c4b99eb142b7de8d5d726231637ab8cdfe511ad60ef6852e8c7602be95e` |
| `tinker_build/analyze` | `f67420e79f8c1793bf491569f77ab4354219926fe5b4c49a2311e64f217f7a28` |
| `tinker_build/testgrad` | `d82b2872136785cadf0c0763aafef60ba1bdd89283f685850485a07a875f6dc2` |
| `tinker_build/dynamic` | `10070a30335cf99e52f73122ef7eb719d0dae1d6f540567ccbb458026bbc2ab7` |
| `tinker_source/params/amoebail22.prm` | `57d10d07b1067d632092d2637ae2a263865a9c8918148d8bd97efb1988c33209` |
| Shared chemistry `reference_chemistry_v1/prepared/PANEL.json` | `50b4e5e03c555f3d729c9df069beea6578b7382e48c6d6d1ff1bf8a1f8882ca9` |

Primary provenance: [La AMOEBA publication](https://doi.org/10.1021/acs.jpcb.2c07237), its exact ACS SI retained by the chemistry agent, and [author-contributed Dy parameter asset](https://github.com/TinkerTools/tinker/blob/87050685eff8840d312e2a332cc82c33f63c7c3d/params/amoebail22.prm). No new source retrieval or numerical calculation was needed to write this checkpoint.
