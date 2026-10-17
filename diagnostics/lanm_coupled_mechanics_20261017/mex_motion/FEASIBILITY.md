# Mex A77–79: conditional local scaffold reference

Read-only feasibility review, 17 October 2026 UTC. No numerical geometry analysis, molecular calls, preparation or submission was performed for this memo.

## Decision

**A connected, chemically complete expansion is plausible, but no complete hydrogen-bond neighborhood has yet been established. Do not launch expanded-region DFT from this inventory.** It could test one local deformation response. It cannot validate the distributed whole-protein MM lowering, and it is not a prerequisite assumed to resolve that response. Root is retaining this option while measuring the existing 209-atom derivative scout.

Actual saved analysis 12983 selected **Glu77–Ala78–Asp79** using both metals, without reference energies or preference signs. This window accounts for **7.4198%** of outside-QM peptide-heavy squared motion. The reported effective participation is **49.27 residues**; 19/46/61 residues cover 50/80/90%. Its nearest heavy-atom separation from the existing QM region is about 5.3–5.4 Å, with four peptide edges to a QM residue and no overlap with existing QM atoms or cap anchors. These are saved descriptors, not newly calculated values. The large reported MM group changes are global energies; neither motion fractions nor endpoint force projections assign an energy to this window. Both searches stopped at 80 accepted steps with residual force, so their endpoints are not stationary reference structures.

## Actual source and existing partition

The normalized 8FNS model contains 2,049 physical atoms: protein, 167 retained waters, EF3 target, and three Nd spectators. The whole charge is +2. The current electronic region has 209 atoms: 205 real atoms including the target, plus four H caps. Its formal charge is 0 and the 1,840-site shifted MM field sums to +2; four omitted boundary CAs remain physical atoms despite carrying no external point charge.

| Existing real membership | Boundary/support |
| --- | --- |
| Complete A84–95 | A83 C/O and A96 N/H complete the bounding peptide units. |
| A66 C/O plus A67 N/H | Existing separate amide support. Saved origin contact: A66 O to A91 N, 2.8321 Å. |
| Whole waters A320, A331, A358, A373, A376, A394, A434 | Other represented waters remain in MM. |
| EF3 target, source A203/ND | La(III) or Dy(III); source identity is preserved separately from substitution. |

Existing caps follow retained/omitted pairs A83 C/CA, A96 N/CA, A66 C/CA and A67 N/CA. Their gradients map to **both physical anchors**. The old preparation's `full_hybrid:false` and spectator “no LJ” prose describes its historical electronic-only stage. The current Cartesian model explicitly retains three Nd(III) ions at EF1/EF2/EF4 (A201/A202/A204), with formal charges and the pinned pure 12–6 TIP3P IOD LJ parameters; their positions were movable. This review proposes no occupancy, spectator-state, proton or water-inventory change.

## Smallest connected whole-residue candidate under the existing boundary rule

Extend the main electronic peptide segment to **complete A77–95**, with **A76 C/O and A96 N/H** flanks. Keep the existing A66/A67 amide support, seven QM waters and target. This adds the complete A77–83 segment rather than making A77–79 a disconnected island. The actual exported graph contains every peptide C–N edge from A76 through A84; A80–83 is the bridge to the current region.

Move only the upstream main cap from **A83 C–CA to A76 C–CA**. The A96 and A66/A67 boundaries remain. A disconnected A77–79 fragment could use additional caps, but would leave its connecting peptide in MM and introduce another charged electronic fragment; it is not the preferred minimum for a coupled local-scaffold test. “Minimum” here means connected complete residues under this boundary convention, not a proven minimal converged electronic region.

The source has deprotonated Glu77 and Asp79 and protonated Lys81 and Lys82; their charged groups are complete in this candidate. Their formal charges balance, so the topology-only extension is consistent with retaining local charge 0. This is **bookkeeping inference**, not an admitted new charge vector or electron count. Any additional hydrogen-bond partners can change it. Actual preparation must independently derive formal charge, electron parity, source mapping, field closure and backend state. Preserve physical La/Dy multiplicities separately from the effective frozen-f reference state; do not reuse the old region's electron counts or cache identities.

## Missing hydrogen-bond closure

The existing `AUDIT.py` selects nearby O/N support within 3.4 Å of the **legacy EF3 core O/N seed**. It records distance-selected contacts, explicitly not established hydrogen bonds. Its only outside-loop protein support is the A66/A67 amide above. It does **not** audit the new Glu77/Asp79 carboxylates, Lys81/Lys82 ammonium groups, intervening backbone, or the new A76 boundary neighborhood. Thus the seven old QM waters cannot be asserted to complete the new network.

A future explicit region preparation would need one common, declared neighborhood rule applied to the saved before/after coordinates for both metals, taking the union before scoring. It must retain complete interacting carboxylate, ammonium, amide and water groups, and document contacts at the new boundary. This requires a bounded geometry calculation on Slurm; no such calculation was done here. If that closure expands substantially or reaches another metal's coordination shell, retain the local-test limitation rather than silently importing a new occupancy/electronic-state problem.

## Required accounting if the option is later activated

- Start from the original ff19SB/TIP3P export and regenerate the partition. Do not append atoms to the already shifted field: A83's old charge redistribution becomes obsolete, and the new A76 boundary needs its own source-derived redistribution. Preserve the stated residue-charge rule; it does not conserve the electrostatic dipole exactly.
- Rebuild retained bonded terms, including CMAP, MM Coulomb/LJ, QM–MM LJ and native exceptions/1–4 handling from the same physical graph. Newly electronic internal interactions must leave the classical ledger. Caps are electronic auxiliaries, not extra physical/MM atoms. Preserve both-anchor force projection.
- Keep the same complete physical coordinates, original FF-plus-formal-ion solvent charge vector and global ALPB cavity at each endpoint. Moving an existing water into QM changes representation, not physical membership. PBE0 embedding replaces the local MACE plus GFN2 embedding difference; no separate QM–MM Coulomb is added again.
- Use raw saved endpoint coordinates and current environments for energy evaluations. The Kabsch-aligned motion coordinates are diagnostic only. Requalify new boundary/field derivatives and state admission before relying on expanded-region forces.

## What the finite local test could establish

If later justified, the plan's four-cell before/after × La/Dy comparison can measure deformation work and restoring-force changes along the saved own-metal paths using one identical expanded-region rule. Compare **within each fixed partition** and then compare those work changes between models. Absolute totals from different partitions are not comparable. A discrepancy can identify sensitivity to the local electronic/MM treatment, boundary or local environment; it is not a unique residue-energy or stiffness decomposition. Agreement supports only this local path.

No outcome from a window representing 7.42% of peptide motion qualifies the global shared MM lowering, remote waters or the remaining scaffold. The saved motion is too distributed to claim that this is the decisive mechanics test. No additional optimizer, clamped release or new geometry is proposed here.

## Exact evidence and implementation pointers

Paths below are relative to the repository root. The actual result retains the complete La/Dy search and point pins.

| Evidence | Path / identity |
| --- | --- |
| Saved selection and descriptors, job12983 | `workspaces/lanm_coupled_mechanics_20261017/mex_motion_v1/execution/RESULT.json`; SHA256 `fee36f58d19a1067a2e18d683e5ec7c503990db7019518c1b620e49b649786c6` |
| Motion and global MM tables | Same directory: `REPORT.md`, `MECHANICS.json`, `WINDOWS.json`, `RESIDUES.json`, `ALIGNMENTS.json` |
| Actual repaired source/partition | `workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2/Mex8FNS/INPUTS.json`; SHA256 `4b3b98983af0cb5fcdbf7b8f52e35063308201ce9df777a30484acf5146b7930`; sibling `core_mapping_A.json`, `boundary_mapping.json`, `environment_atoms.json` |
| Native force-field source and graph | `workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1/Mex8FNS/export_v1/EXPORT.json`; SHA256 `b9005c860b2795807427c275c6dc9c4aff96352e2572db375bc65b819aa6275f`; sibling `atoms_charges.json`, `bonds.json`, `protein_water_system.xml` |
| Source proton inventory and IDs | `workspaces/lanm_global_occupancy_20260923/prepared_v1/sources/Mex_8FNS/atoms.json`; SHA256 `adfa05999245c5c9175df3d7c2ff33cb2166f59adcdc2810af5e370b36fbc4c4` |
| Existing contact selection and saved contacts | `diagnostics/metal_environment_response_20260926/lanm_preparation_feasibility/AUDIT.py` and `RESULT_v2.json`, Mex row `support_contacts` and `proposal` |
| Existing boundary/charge construction | `diagnostics/metal_environment_response_20260926/lanm_ef3_preparation/PREPARE_TRANSFER.py` |
| Actual current MM construction/state | `scripts/metal_environment_component_checks.py::build_system`; `diagnostics/lanm_coupled_mechanics_20261017/mex/MEX_CARTESIAN.py::MexEvaluator` |
| Existing embedded DFT assembly | `workspaces/lanm_coupled_response_20261016/mex_dft_origin_assembly_v1/execution/RESULT.json` |

This memo changes no source, parameters, prepared region or execution plan.
