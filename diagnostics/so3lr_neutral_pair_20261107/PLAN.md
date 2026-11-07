# Neutral cis/trans pair and bounded solvent feasibility

Declared7 November2026 UTC, before new potential evaluations. Jacob's active
goal is to proceed until there is a good reason to stop. The previous charged
prototype passed L's energy test, while force sensitivity remained. Its result
does not decide behavior of the intended neutral, solvated ligand system.

Use cis599/trans600 only. No new experimental pair search, reserved outcomes,
DFT, force-field fitting, model training or replacement of M1. Keep original
results and failures. Pinned L-2606/source70707f5/GLP2294f4f/environment unchanged.
M was the prior sensitivity model and is not substituted for L.

Preparation:
- Exact stereochemical structures from existing curated rows and geometry traces.
  Fully deprotonated carboxylates, neutral tertiary amines: each ligand L4-.
  Verify stereochemistry from3D, atom order, formal charge and covalent distances.
- Reuse conformers0/1/2, assignment0,k20: original La/Dy/La size-prepared
  coordinates respectively for BOTH ligands. These are initial preparations,
  not measured structures or Tauri-optimized minima. No selection by Tauri energy.
- Site0=trans600,site1=cis599; metal positions(11,21,21)/(31,21,21) A in a42-A
  cubic cell. Apply predetermined rigid rotations to vary orientations.
  A=La-trans+Dy-cis; B=Dy-trans+La-cis. Only the metal identities change.
  Two K+ counterions at(21,10,21)/(21,32,21): total charge0, multiplicity6.
- Three independently packed arrangements, seeds44100/44101/44102,2400 real
  OHH water molecules each. Use archived water geometry, no virtual sites.
  Packmol20.15.0 commit7ddbd6562301aebf2a4980599ebb3ee4191020cd is built in the
  new workspace because installed20.14.3 lacks PBC packing. Existing tools/envs
  remain unchanged. Periodic42-A packing,tolerance2A; check actual minimum-image
  contacts and restore only rounded fixed-solute coordinates from exact inputs.
- This is a zero-added-background-salt DEVELOPMENT system with neutralizing K+,
  not the unverified1-M experimental medium. No experimental agreement claim
  until the exact electrolyte and original source records are established.

Checks before dynamics:
1. Neutral finite first-shell extracts from the prepared boxes, including both
   ligands/metals/K and every water whose O is within3.5A of either metal.
   Preserve variable denticity; never force six donors to coordinate. Repeat the
   spectator control at200/400A in two orientations on all three preparations.
   Also translate the whole trans ligand +/-0.1A along its nearest carboxylate-O
   radial direction to test cancellation under a local motion. A/B share exact
   coordinates. Primary max|spectator-induced change in swap energy|<=0.5kcal/mol;
   forces recorded separately. A failure in this extract does not itself reject
   the full solvated system, whose global conditioning differs.
2. Full periodic42-A boxes: L,float64,head5,Q0/five unpaired,standard engine,
   SR5A,LR12A/dispersion damping2A,PME smearing4A/grid1A/order4. Compare endpoint
   energies and forces on all starts. Check rigid translation/wrapping invariance.
   Remove the same three independently selected remote OHH molecules one at a
   time from A/B (farthest from all solute heavy atoms, outside12A if available).
   Record reaction-energy/force changes. Use0.5kcal/mol as the same practical
   robustness target; diagnose legitimate electric interactions rather than
   silently treating every change as model error. No fragment-charge overrides.

Feasibility, only when the actual system is evaluable without a demonstrated
large spurious response: measure H200 throughput on that real box first.
Up to six unrestrained endpoint starts (3A/3B),298.15K,initial timestep0.25fs.
No association walls initially: observe whether bound basins sustain themselves.
Monitor metal/donor distances, water occupancy/exchange, coordination, ligand
stereochemistry/bonds, proton transfers, ion partitioning and independent-start
agreement. A trapped trajectory or dissociation is not automatically a potential
failure, and endpoint mean energies are not the requested free energy.
No automatic lambda ladder, binding-free-energy estimate or protein transfer.
Any later alchemical run requires explicit bound domains and restraint-release
assessment; unrestricted metal permutation must not be mistaken for selectivity.

Resource envelope: Slurm standard_expand, root owns submissions. Initial
preparation/static CPU budget4 allocated CPU-hours; first solvent allocation
at most one H200 GPU-day TOTAL including benchmarks/failures. Record measured
costs, use internal finite timeouts; no scheduler memory/time overrides. Respect
the live dnv1_fold allocation. Scheduler waiting is not a scientific failure.
Declare any necessary preparation repair, preserve its failed attempt and do
not change criteria to rescue a selectivity result.
