# Exact-source scaffold feasibility — 28 September 2026

**The current Hans8DQ2 source already has complete standard-protein/water mechanics and exact atom maps. A new parameter build is unnecessary.** A covalent-preserving finite additive model is technically feasible, but boundary double counting, metal cross interactions and solvent treatment remain scientific qualifications. Do not start unconstrained whole-protein relaxation from this inventory alone.

## Reusable current artifacts

Base directory `workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1/Hans8DQ2/`:

- `export_v1/protein_water_system.xml`: native OpenMM ff19SB/TIP3P System,1,887 atoms, no constraints;1,820 bonds,3,106 angles,4,239 proper/improper torsions,108 CMAP interactions/16 maps,9,289 nonbonded exceptions. No ion particles/parameters are silently present.
- `export_v1/atoms_charges.json`, `bonds.json`, `EXPORT.json`: exact source IDs, normalized coordinates, native charges, topology and pinned force-field/source provenance. Read-only comparison confirms System topology order exactly matches this atom array. These are the current1,665 protein atoms plus74 complete waters, not a rebuilt tail or damaged relaxed candidate.
- `core_mapping_A.json`, `core_mapping_B.json`:190 real QM nonmetal atoms, target EF3 and four synthetic caps. Retained/omitted source coordinates and analytic cap Jacobians are present.
- `boundary_mapping.json`, `environment_atoms.json`, `INPUTS.json`: local charge redistribution, four omitted-MM1 CA charges with physical atoms retained, actual La/La/Na spectators, source receipts and exact A/B correspondence.

Existing `scripts/metal_environment_mechanics.py` and `metal_environment_components.py` implement an earlier1H4I ledger/components; reuse their term classification, exception handling and chain-rule patterns, **not their1H4I products, source-size assumptions or already chosen boundary policy**. No new engine was built in this review.

## Exact subtraction inventory

Read directly from the current XML and selected190 source IDs:

| Native term | Wholly QM | Cross boundary | Wholly MM |
|---|---:|---:|---:|
| Bonds |186|4|1630|
| Angles |323|20|2763|
| Proper/improper torsions |457|32|3750|
| CMAP |12|4|92|

Wholly-QM classical terms must be excluded when their physics is already in the electronic region. Keep MM-only terms with their original dimensional parameters and CMAP surfaces. Cross-boundary terms need a declared policy: the older1H4I implementation retained all such terms, whereas native ORCA by default removes certain link-associated terms containing only one MM atom (`DeleteLADoubleCounting`/`DeleteLABondDoubleCounting`). The cap already creates artificial local bonded energy. Therefore copying all cross terms and claiming exact cancellation would be unjustified. See the existing independently reviewed `lanm_preparation_feasibility/BOUNDARY_REVIEW.md`.

A candidate expression is embedded electronic energy(realQM+mappedcaps; shiftedMMcharges,DoEQfalse), plus retained MM/cross bonded terms, MM–MM Coulomb/LJ and realQM–MM LJ. No cap force-field particles, QM–QM classical terms, classical QM–MM Coulomb, C4 induction, second D4/gCP or duplicate solvent energy. The four zero-charge CA atoms retain real native bonded/LJ terms. Recompute shifted-charge MM exception products using the native graph/1–4 scale; preserve native exception LJ overrides. At moved coordinates regenerate cap Jacobians, and combine actual electronic nuclear and pointcharge derivatives before projecting to physical atoms. Pointcharge forces alone are not whole-scaffold forces.

## Metal parameters: available is not qualified

Installed file `/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/dat/leap/parm/frcmod.ions234lm_126_tip3p` contains **pure12–6** IOD/TIP3P candidates (not C4-stripped12–6–4):

| Ion | Rmin/2 Å | epsilon kcal/mol |
|---|---:|---:|
| LaIII |1.718|0.15060822|
| DyIII |1.609|0.08389240|

The adjacent `frcmod.ions1lm_126_tip3p` contains NaI1.475Å/.03171494kcal/mol (HFE/TIP3P). These exact file values were inspected; no parameters were assigned or evaluated. Standard real QM protein atoms/waters already have source-native LJ, so no PQQ/cofactor typing is needed in this LanM region.

These ion values supply a declared candidate for target–exterior and spectator repulsion/dispersion, but hydration fits do not establish compatibility with the new embedded frozen-f electronic surface or accurate metal–protein cross response. No validated La/Dy-specific QM/MM cross-LJ set was found in the current prepared ledger. Do not use zero terms or invent covalent ion bonds. Existing legacy12–6–4 Hans trajectories used different topology/occupancy and failed their discrimination test; their C4 induction cannot be added on top of explicitly represented electronic polarization as a shortcut.

If spectators remain fixed, their classical interactions with moving protein atoms still matter. If allowed to move, a fixed monopole without qualified short-range interactions is particularly inadequate. Keep the actual EF4Na identity and the declared occupancy; no conversion to four La sites.

## Smallest useful next test

After the active electronic force gate, first implement only a **classical source-decomposition check** using this exact XML and existing A/B coordinates, with all ion terms explicitly unavailable until their chosen candidate is declared. Compare full native protein/water mechanics with the sum of classified terms to verify subtraction/exception bookkeeping and inspect source strain. This is a numerical ledger test, not a La/Dy selectivity gain. Current Asp85 A/B motion is entirely inside QM and does not move the distant caps; it cannot qualify scaffold compliance by itself.

For the first actual compliance question, select one covalent-preserving small peptide torsion adjacent to EF3 from the real source graph, include attached hydrogens, and move the same physical atoms for both targets. Compare a finite set of small signed displacements under the declared parent mechanical energy and combined electronic energy, retaining affine source forces rather than assuming equilibrium. Establish the chosen boundary policy and cross-LJ response on this finite test before any coupled optimization. No arbitrary springs, elastic-network energy scale, whole-protein Hessian or entropy term is needed.

The74 retained crystal waters are not bulk solvation or a periodic equilibrated solvent bath. A successful finite-field mechanical test would still be a local conditional scaffold response. A defensible solvent-consistent full-protein relaxation model, dimerization/occupancy thermodynamics and effective stiffness remain unqualified. Source native covalent parameters provide a useful starting stiffness, not proof that the capped frozen-f force surface transfers.

This review read existing source artifacts and native parameter files only. No force/energy calls, new structures, optimization, reserved outcomes or Slurm submissions occurred.
