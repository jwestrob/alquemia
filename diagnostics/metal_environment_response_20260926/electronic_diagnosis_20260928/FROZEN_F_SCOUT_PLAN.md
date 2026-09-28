# Frozen-f origin force scout — declared before execution

Question: does the audited ECP55 valence model yield a converged energy and analytic
gradient on the exact consumed Hans EF3 50-atom Dy core? One origin endpoint only.
Protocol nikasha_DyIII_frozen4f_pbe0_force_scout_v1. PBE0-D4/def2-TZVP CHNO,
lcecp-1-TZVP Dy matching ECP55 plus actual ORCA6.1.1 AutoAuxJ, RIJCOSX,
DefGrid3 VeryTightSCF EnGrad, PModel, no solvent/embedding, no gCP, no iteration
increase. Charge−1; physical DyIII4f9 sextet, effective restricted valence singlet,
208explicit electrons. Existing ECP28 adapter is unchanged.

Use exact pinned archived source and basis assets; no geometry optimization or
unconverged-wavefunction restart. Native convergence, ECP55, electron count,
matching engrad energy/coordinates and analytic SCF/ECP/D4 components required.
Convergence alone does not qualify force consistency, 4f approximation or affinity.
Absent results stay absent; no historical score/reference reuse.

After admissible origin, prepare separately declared central differences along
one actual Dy-to-nearest-donor axis with fixed ligand geometry (metal motion is a
physical coordinate; caps do not move). Steps0.005 and0.0025 angstrom; compare
analytic directional derivative to both central differences, tolerance0.1kcal/mol/A
and refinement change0.05kcal/mol/A. These are development force checks on a small
real coordinate perturbation, not a full numerical gradient or a label fit.
Do not launch these dependent endpoints until origin proof is inspected.

Resource policy: one40MPI endpoint on currently available standard-shared CPUs,
--mem=0 with the existing observed-available-memory renderer. Normal priority,
no PQQ dependency/interference. Live SCF health and terminal root wake armed.
No automatic retry, iteration escalation, whole-protein run or production change.
