# Training Hamiltonian and reference-derivative context

20 October 2026. Source, metadata and primary-paper review only; no model loading, numerical replay or new calculations.

## The force comparison uses different electronic models

The pinned MACE-OMOL-0 checkpoint is the official100M OMol model with charge and multiplicity conditioning. Its release says “MACE model trained on omol 100M with spin and charge embedding.” Installed MACE metadata lists its target as `DFT (wB97M-VV10)`; the original local plan explicitly records vacuum **ωB97M-V/def2-TZVPD**.

The [OMol25 primary paper, Section2.6 and AppendixA](https://arxiv.org/html/2505.08762v1#S2.SS6) confirms ωB97M-V/def2-TZVPD, ORCA6.0.0, RI-J/COSX and DEFGRID3. Diffuse functions were selected because the dataset includes anions. Lanthanide-containing singlets and non-singlets used UKS; singlet metal calculations used a symmetry-breaking initial orbital rotation. The paper also discusses grid-related energy/force inconsistencies and its chosen numerical settings. This establishes the dataset recipe, not numerical qualification of our particular peptide or proof that the checkpoint reproduces its training targets exactly.

Our whole256 reference uses **PBE0-D4, CHNO def2-TZVP/def2-J, RIJCOSX, DefGrid3 and VeryTightSCF in ORCA6.1.1**, with separately pinned La/Dy basis/ECP46/55 assets and effective singlets. Physical multiplicities remain La1/Dy6; the Dy reference freezes4f9 and omits SOC. Exact metal-basis/ECP equivalence to the OMol dataset has not been established. The original whole-reference specification explicitly states that PBE0 is “not a claim of reproducing the exact MACE training functional or spin representation.”

Thus a persistent MACE–PBE0 force difference can contain learned-model error, electronic-method/basis/state differences and reference numerical error. The present evidence does not uniquely assign their shares. PBE0 is an independent conditional comparator, not asserted physical truth. Adding the same saved ALPB term does not resolve those differences or independently validate solvent.

## Exact local evidence

Paths below are relative to the repository root. Hashes are SHA256 of the inspected files.

| File | SHA256 |
| --- | --- |
| `workspaces/mace_omol_20260917/software_v1/release.json` | `1855c8105133f34f98cb387c980c42abfacebb07a59505d8d09b0edaaf7d6fa5` |
| `workspaces/mace_omol_20260917/software_v1/software_manifest.json` | `7d2448b279b5191a529ef329a60442635601576794f9737efc64343c9a6e0168` |
| `workspaces/mace_hybrid_20260916/software_v1/venv/lib/python3.11/site-packages/mace_torch-0.3.16.dist-info/METADATA` (line317) | `ce3a831df793896eb140235a93dc0281e161db56bb89a0790b9d1cae3154935d` |
| `diagnostics/mace_omol_20260917/PLAN.md` (line32) | `f1c72cd1bde260170e949ce46a4ce0fcf64c6f48a45fa8664e8e0cd8138d0c83` |
| `diagnostics/lanm_coupled_mechanics_20261017/lbt_whole_reference/SPEC.md` (lines30,139) | `02b45065533355f06f53316e201e165971e0736f77130dde97dcc718bef2af75` |
| `diagnostics/lanm_coupled_mechanics_20261017/lbt_whole_7ccn_reference_v1/SPEC.md` (lines15–17) | `a8cb688c8d5f8f1f109a6c1436da24e68410cb92528f503d123cc83650307447` |

The manifest binds `workspaces/mace_omol_20260917/software_v1/MACE-omol-0-extra-large-1024.model`, SHA256 `9b64b4fd5153ca578c694abc57806d8111050de6ff652e695c9b525bc4d36469`. This review read its existing metadata; it did not reload or rehash the checkpoint.

## What remains to qualify

The inspected whole256 reference lineage has **no saved native finite-difference qualification**. Both original four-cell collectors retain `force_consistency_qualified=false`; their specifications explicitly exclude an independent derivative test. MACE/ALPB numerical checks, native EnGrad admission and the separate failed embedded-Mex derivative experiment do not supply it.

The three Dy force supplements15720 and saved localization15797 are already complete. They show substantial but incomplete cancellation of the broadly distributed force discrepancy between the metals; they are not pending tasks and do not qualify the reference derivative. See the [completed force report](../lbt_dy_force_completion_v1/FINAL_REPORT.md), SHA256 `41b204482f2da60e4e30afc430edcc6adfabf41983a6133cae6f12cbfac10014`, and [localization report](../lbt_force_error_localization_v1/FINAL_REPORT.md).

The next declared action is only the [two-single-point directional test](SPEC.md), SHA256 `c79d67370983f6100bc9e3ce3fea45a5641af26a4f924fdc3af41efd84f8c609`: reuse the step725 central gradient; evaluate exactly `x±0.01Å*u`, where `u` is the normalized full vacuum force discrepancy; retain the existing PBE0 method/state/settings and frozen numerical rule. This tests one sampled derivative and cannot validate every component, electronic state, solvent, a physical minimum or preference. It authorizes no matched-training-functional calculation, new optimization, retry or settings sweep. Any later Hamiltonian-attribution calculation would require a separate declared scope and verified method/state assets.
