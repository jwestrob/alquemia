# Native solvent-force consistency after fixed self continuation

## Approved scope, frozen before molecular calls

Jacob, September23: “keep on goin ... launch subagents to do what experiments you
like ... goin to bed, email when cool stuff”. Parent assigned this finite test:
4MAE and Q88JH5 original complete contexts; both metals; vacuum and ALPB(water).
No new geometry optimization, MACE, DFT, biological labels, classifier calibration,
production change or extension of the unsuccessful protein-only scaffold model.

Existing analytic-gradient work (`accommodation_response_20260920`) already showed
that the installed native GFN2 driver can match energy derivatives. The later
full-composite search had failed numerical curvature checks. Confirmed native
restarts now remove one 4.805 kcal electronic discontinuity. This experiment asks
whether that **same fixed continuation policy** produces consistent physical
solvent derivatives; it does not claim analytic capability as a new discovery.

## Physical coordinate and exact states

Use existing source-backed Kinematics and its coupled source/cap Jacobian:
4MAE anchor Glu172 chi3 and Q88JH5 anchor Glu221 chi3. This rotates the complete
terminal carboxylate about CG–CD. No independent cap coordinate, donor identity,
metal, water, protonation, heavy-bond length or context-membership changes occur.
The former is an existing eligible mode omitted by the later adaptive four-angle
selector; this is explicitly a force diagnostic in the existing full map.

For each source and metal, q=0 and q in {-0.001,-0.0005,+0.0005,+0.001} radian.
All other existing physical coordinates stay zero. Geometry replay, Jacobian,
source bonds, fixed atoms, no new severe overlaps and displacement <0.8Å must
pass before execution. Preserve paired nuclear coordinates/charge policy.

The eight q0 energy/gradient cells are **reused from stage2 of the coordinated
native_pool_continuation_20260923 experiment**. They are not rerun here. This
branch executes exactly32 displaced cells, each initial NoAutostart followed by
two own-result continuations, at most96 calls. The final continuation always
supplies the reported energy, never the lowest/favorable stage. Failure prevents
its dependent continuation; every fixed cell remains in the denominator.

## Native recipe and electronic qualification

ORCA6.1.1 native GFN2, native mixer true, 300K, MaxIter500, existing native
parameters, exact charge/singlet and solvent switches. No ordinary SCF or claim
that generic TightSCF tightens the special mixer. Initial runs have fresh folders
and NoAutostart. Continuations copy only their own preceding matching GBW/xtbw to
the fresh runtime basename. Actual XTBRESTART, native mixer, parameters, charge,
electron count and normal SCF termination must be confirmed from the output.
Each state and stage retains original seed hashes/receipts; no adaptive retries.

Stage1→2 energy change must be ≤0.1kcal/mol per cell. Within each geometry the
Ca−La difference of solvent transfers must change by ≤0.2kcal/mol. Report
individual solvent-transfer changes too. These are separate numerical gates;
normal termination alone is insufficient. No ground-state guarantee is inferred.

## Force comparison, frozen acceptance

For each vacuum/ALPB component, each metal's ALPB−vacuum transfer, and the
Ca−La transfer contrast, calculate centered derivatives at h=0.001 and h/2.
Use raw unrounded energies and convert Hartree to kcal/mol exactly once.
Project the actual q0 Cartesian analytic gradient through the existing context
Jacobian, including cap chain rule. Hartree/bohr becomes kcal/mol/Å; derivative
is dE/dq, not negative force. Compare the gradient's own energy with the printed
final energy, exact atom order and q0 coordinates. Require native analytic SCF,
CN and correct solvent-driver evidence; numerical gradients are unsupported.

Both |D_h−D_h/2| and |g_projected−D_h/2| must be no larger than
max(0.2kcal/mol/radian,0.05*|D_h/2|), separately for every declared quantity.
Record kcal/mol/Å by dividing by the RMS q0 physical-heavy Jacobian speed over
actually moving heavy atoms (norm >1e−10Å/radian). Also retain the speed over all
mapped heavy atoms, moving-atom IDs and actual displacements, so normalization
cannot hide a large angular discrepancy. These are engineering force checks,
not biological accuracy or entropy estimates.

If inconsistency occurs, inspect retained electronic components, occupations,
finite-temperature information and derivative energy identity. Do not add an
invented entropy term, select another root, tune the step or extend the solver.
An unqualified derivative ends this branch's accommodation-enabling claim.

## Execution and reporting

Existing pinned native executor, eight concurrent eight-rank workers, one64CPU
128GiB allocation using the CPU-sharing GPU partition route, no GPU requested.
Finite manifests and stage receipts; preserve all errors, actual starts and costs.
No project wall-clock/compute cap beyond the finite agreed scope and scheduler
policy. Outputs under workspaces/native_solvent_force_20260923. Report tests,
all40cells including eight separately owned q0 cells, source/implementation pins,
numerical and derivative gates independently, and a vault note. No production
promotion follows from this diagnostic.

Parent coordination before execution: approved Glu172 chi3 / Glu221 chi3 as an
eligible-mode diagnostic despite the former's absence from the adaptive four.
The fixed32cells/96calls and eight externally supplied q0 cells are unchanged.
