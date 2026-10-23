# Paired explicit-water exchange legs: prepared, awaiting a qualified Hamiltonian

The aqueous geometry is prepared and the existing bound geometry is preserved. **Both legs are neutral, and each has an exact La/Dy atom mapping. No energies, forces, minimization or dynamics were run.**

## Prepared systems

| Quantity | Bound LBT3 | Aqueous ion |
| --- | ---: | ---: |
| Peptide atoms | 257 | 0 |
| Metal ions | 1 | 1 |
| Waters | 3,894 | 3,980 |
| Na+ / Cl− | 8 / 7 | 7 / 10 |
| Total atoms | 11,955 | 11,958 |
| Net charge | 0 | 0 |
| Cubic starting cell | 49.324 Å | 49.324 Å |
| Native metal index | 258 | 1 |

The bound peptide is the deposited amidated I10 LBT3 construct, NH3+–FIDTNNDGWIEGDELLA–CONH2, with the five Asp/Glu side chains deprotonated: peptide −4, metal +3, complex −1. One extra sodium neutralizes that complex. Three extra chlorides neutralize the aqueous +3 ion. Each leg also contains seven background NaCl pairs, corresponding to **0.0968662 M background salt pairs** in the starting cell. This does **not** equate total ionic strength, ion activity or chemical potential. The mobile monovalent ions alone contribute 0.103785 M and 0.117623 M to ionic strength, respectively; these numbers exclude the metal and peptide.

The aqueous fill was rebuilt from the same pinned 500-water packing and tiling recipe, with only the metal excluded from the water volume. It does not retain a peptide-shaped cavity. Three of the 4,000 tiled waters overlap the metal exclusion region; seventeen further complete waters are replaced by ions. Deterministic farthest-site placement and the original 5 Å ion–ion separation rule were retained. This deliberately dispersed initial ion atmosphere requires equilibration.

The smallest aqueous metal–water-O separation is 2.89176 Å; ion–ion separation is 5.12751 Å. Minimum aqueous Na/Cl–water-O distances are 2.56969/2.56279 Å; the corresponding water-H minima are 1.70918/1.73060 Å. The declared gross-overlap gates all pass. The bound fill has no water O within 5.06324 Å of the metal; this initial packing does not establish equilibrium coordination or water occupancy.

An independent saved-coordinate audit passes 111,392 source/topology/charge/identity checks. Water and ion coordinates reconstruct from the pinned source with maximum residual 7.1×10−15 Å. All waters still have the packing geometry, approximately 0.95720 Å O–H and 104.520° H–O–H. The assigned flexible AMOEBA water has equilibrium angle **108.50°**, so these are unequilibrated starting coordinates and need a declared initial relaxation.

## Exchange definition and atom mapping

The intended direction is La → Dy in both legs:

`ΔΔG = ΔG(bound, La→Dy) − ΔG(aqueous, La→Dy)`.

A negative result favors binding Dy relative to La under the eventual model and solution convention. Within each leg, every nonmetal atom, atom order, bond, charge and coordinate is identical between the initial La and Dy files. The metal is index258 in the bound leg and index1 in the aqueous leg. Separate thermodynamic legs need no atomwise map from peptide atoms to aqueous solvent. `ATOM_MAPPING.json` records every native index, peptide source join, and water/ion source molecule/tile.

Both boxes remain net-neutral at every λ; the metal charge remains +3. **No charge-changing alchemical correction is being applied.** Net neutrality and unchanged metal charge do not cancel all periodic electrostatic effects. Ion-atmosphere, polarization, box-size and finite-concentration effects need explicit later assessment. Seven background pairs and unlike neutralizing ions are a finite-box convention, not a proof of matching experimental salt chemical potential.

Use the same fixed metal mass, 138.91 Da, throughout both legs; retain native nonmetal masses and no hydrogen mass repartitioning. For an initial fully flexible classical calculation, use identical unconstrained treatment and global translational-motion removal in every endpoint/window. The resulting kinetic factors are λ-independent. This statement does not cover isotope effects, nuclear quantum effects, different constraints, or mass-dependent restraint conventions. With all atoms active, no constraints and three removed translational degrees of freedom, the intended thermostat counts are 35,862 bound and 35,871 aqueous; verify actual native counts before dynamics. Do not remove global periodic rotations.

## Concrete missing pieces before simulation

1. **Qualified protein metal parameters.** Both input files deliberately retain undefined reserved metal types900/901. The small-donor repair is not automatically a protein repair. Protein carboxylate O143/161 (class31, polarizability0.837 Å³) differs from small acetate O403 (class403, 1.931 Å³). Protein amide O5/11/153 (class5, 0.837/0.834 Å³) differs from small amide O105 (class55, 0.914 Å³), including permanent multipoles. Small acetate type403 also collides with protein water H403. A new explicit donor mapping and transfer qualification must precede protein energies. The fixed-charge alchemical path must be regenerated and checked for the mapped protein and ion classes, including Na/Cl pairs.
2. **A complete periodic Hamiltonian.** The existing `protein.key` specifies only the parameter file, one thread and print precision. It does not select Ewald, cutoffs, grids, polarization convergence or a metal model. No executable periodic key was prepared here. Water types402/403 (classes73/74), Na405/class76 and Cl414/class85 are fixed to the pinned `amoebabio09.prm`; all permanent multipoles, induced polarization and intramolecular terms must remain matched across the two legs.
3. **Equilibrated solvent, ions and protein.** Packing is not an equilibrium sample. Initial water strain, ion mixing, bound water exchange, donor coordination and peptide conformational populations require sampling. An initial minimization or short stable trajectory cannot settle these quantities.
4. **Observable and finite-box mapping.** MES is absent. The prior observable contract still separates free-ion exchange, buffer/speciation and finite-concentration assay effects. A 283.15 K development calculation corresponds to the recorded 2020 ITC temperature; it must not be presented as the 2025 amidated assay temperature, which remains unverified. Protonation and the different assay targets remain separate uncertainties.

## Proposed first stability pilot — finite and conditional

Execute only after one candidate passes the predeclared physical chemistry and protein donor-transfer gates, its exact pair parameters are bound to both legs, and the GPU/native numerical checks pass. This is a proposed pilot, not a submission or a new qualified model.

- **Four endpoints:** bound-La, bound-Dy, aqueous-La, aqueous-Dy. Same four pinned starting geometries, water/ion force field and fixed metal mass. No metal–donor bonds, positional restraints, frozen atoms, HMR or bond constraints.
- **Explicit periodic settings to bind before execution:** 49.324 Å fixed cube; conducting/tinfoil electrostatics, 9 Å real-space Ewald cutoff, permanent and polarization Ewald coefficients0.4 Å−1, 64³ PME grid and order5 for both; mutual polarization tolerance10−7 and at most100 iterations. Buffered14–7 VDW cutoff9 Å, taper start8.1 Å and no long-range VDW correction for this finite stability test. These choices need later cutoff/grid/size checks before production free energies. In the pinned native code, `EWALD-BOUNDARY` selects vacuum regardless of following text; tinfoil is the pinned default, so omit that keyword and verify the resolved boundary. Do not write `EWALD-BOUNDARY TINFOIL`.
- **Initial preparation per endpoint:** one periodic energy/gradient evaluation; then one all-atom minimization capped at500 iterations, target RMS gradient1 kcal mol−1 Å−1. If the cap or native failure occurs before that target, stop that endpoint and preserve the result; do not silently extend. Record every actual energy/gradient evaluation and final coordinates.
- **Trajectory per admitted endpoint:** 4,000 velocity-Verlet steps of0.25 fs, fixed-volume Bussi thermostat at283.15 K with0.1 ps coupling time; total1 ps per endpoint. Use declared seeds2101/2102/2103/2104 in the order above, remove total translational motion every100 steps, and save coordinates/energies every100 steps. Verify thermostat degrees of freedom and runtime settings from actual native output. No exchange estimator is computed from this pilot.
- **Finite cost and stop:** at most12 scientific program invocations: four initial evaluations, four capped minimizations, four finite trajectories; later stages depend on each endpoint's prior pass. Record actual costs and internal evaluation counts. Require finite outputs, converged mutual polarization at every evaluated state, actual successful minimizer/trajectory termination, and all4,000 requested steps. Water O–H outside0.75–1.25 Å or H–O–H outside50–170° is a declared gross-instability stop. Record metal coordination, donor identity, ion association, energy and temperature; stable bond lengths or successful completion alone do not establish a physical binding basin or equilibrium.

After this pilot, any free-energy calculation still needs equilibration, λ-overlap, independent replicas, forward/reverse agreement, uncertainty and convergence checks, and the finite-size/assay mapping above. The existing full strategy remains in `../PLAN.md`.

## Evidence and cost

Prepared inputs and complete atom mapping: `workspaces/lanm_selectivity_free_energy_20261021/paired_legs_v1/prepared_v1/`. The parent box and source pins are preserved in `RESULT.json`. Independent audit: `paired_legs_v1/AUDIT_v1.json`. Preparation took7.50 s and saved audit1.79 s, each on one local thread; **zero native calls or scientific allocations**.

Native source used for the proposed settings is the already pinned `gpu_engine_feasibility_v1/tinker9_source`: `tinker/source/kewald.f` (PME/boundary), `cutoffs.f` (cutoffs/taper), `prmkey.f` (polarization controls), `mdinit.f` and `mdrest.f` (degrees of freedom/inertia), and `lbfgs.f` (iteration cap). These source files and the report are pinned in `READY.json`; no implementation or force-field source was changed.
