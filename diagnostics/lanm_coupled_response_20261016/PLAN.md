# LanM coupled response: finite preparation and mechanical control

16 October 2026. Root owns integration and submissions; Leibniz owns PREPARE.py.

## Goal and immediate question

Develop a useful La/Dy profiler combining local electronic chemistry with
whole-protein mechanical response, recovering supported relative preferences
across credible structures with explicit occupancy and hydration. Completing
this preparation does not complete that goal.

The saved Hans source-force diagnostic (Slurm11788) identifies Asp85/Asp87
differential loads but also substantial common strain. Can a small common set
of covalent-preserving motions generate admissible own-metal downhill proposals
with the full existing scaffold mechanics included? This stage checks actual
geometry and classical energy response; it does not calculate new QM energies.

Jacob's latest steering favors experimentally profiled monomers as the first
predictive benchmark. Recover those before committing to a larger Hans/dimer
campaign. Hans remains a useful consumed source-response control. The real
partner Arg100 contact network is a separate physical hypothesis, not an
assumed explanation or prerequisite for the monomer benchmark.

## Exact inputs and physical state

Pin the complete saved PROFILE.json and four ASSEMBLED.json outputs of11788,
their exact repaired Hans8DQ2/Hans8FNR INPUTS.json and existing full classical
ledgers. Preserve all real atoms, source-specific waters, spectators, proton
inventory, charges, caps and parameters. EF3 is the exchanged La/Dy target.
Source compositions differ; no cross-source minimization or population estimate.
No new force-field parameters, protonation, solvent term or quantum Hamiltonian.

Reuse the existing native protein graph, physical/cap maps, ff19SB/TIP3P and
pure12-6 spectator/cross-region ledger. Electronic embedding and classical
electrostatics retain their existing exclusions; no duplicate QM-MM Coulomb.
The existing force/state limitations remain explicit, including frozen-f Dy.

## Frozen candidate construction

Use the SAME eight physical directions for both metals within each source:
Asp85 CA-CB and CB-CG torsions, Asp87 CA-CB and CB-CG torsions, the N83-CA83
downstream protein hinge, and three EF3 target translations. Native graph cuts
define complete moving components and their attached hydrogens. Waters and
spectators remain fixed. Synthetic caps follow their physical anchors and never
become independent coordinates. Context membership remains fixed.

Normalize Jacobian columns by physical Cartesian displacement. Remove redundant
directions with a normalized-Gram relative pseudoinverse cutoff of1e-10. Project
each metal's own full electronic-plus-classical gradient into that common span;
never minimize the metal difference. Report individual, common and differential
loads, rank and conditioning.

Construct one deterministic proposal per metal, with a single conservative scale
limiting every physical atom's displacement to0.02angstrom and each torsion to
0.5degree. Use the triangle bound for composed rotations. No label-dependent
direction, line search, retries or halving. Failures remain unavailable.

## Admission and actual work

Source graph torsions preserve native bond lengths/angles within1e-10angstrom/
radian and tetrahedral signed volume within1e-9angstrom^3, without inversion.
Retain the existing source-relative severe-contact rule (heavy1.0angstrom,
H0.55angstrom): reject a nonexcluded contact below min(cutoff, source distance)
minus1e-7angstrom. Exclude native1-2/1-3 pairs. Preserve and report pre-existing
short contacts; do not weaken gates after inspecting a result. Check peptide,
chirality, state, inventory, finite coordinates and cap/boundary reconstruction.

Evaluate the existing complete classical ledger at both candidates for both
metals: two sources x two proposals x two metals = eight new classical cells
(four component queries each). Reuse four origins. Compare finite classical
work with its saved-gradient prediction. Electronic linear work is explicitly
a local diagnostic, never a new energy or affinity score.

Emit explicit candidate INPUTS with actual moved physical exports, point-charge
coordinates and cap anchors. Each singleton candidate may use identical A/B
aliases, but the original geometry remains separately pinned. Do not label a
multimode displacement as the old single-torsion perturbation.

## Execution and next decision

One standard Slurm CPU, existing lanm_qmmm environment; no memory/time overrides.
Pin plan, code and inputs in a fresh workspace, retain all failures, record costs
and arm the actual terminal wake. No login computation or new molecular QM/MACE
call. No PQQ/default/calibration changes, and no reserved numerical labels read.

After inspecting results, choose the smallest scientifically useful matched
response experiment on the monomer benchmark or a justified Hans control.
Quantum evaluation of this eight-cell proposal matrix is not automatic: first
resolve the relevance of the source state and the available monomer targets.
This is an experiment-selection decision, not a per-analysis permission gate.
