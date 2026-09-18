# Physical donor coordinates from saved matched-hybrid gradients

Declared2026-09-18 under the active autonomous MACE goal. The preceding
metal-only response failed structural transfer (4/12 raw directions; no new
success on2FW0/2FVY). Do not expand/tune that candidate or reinterpret its gates.

Question: do the actual saved hybrid gradients support testing coupled donor
motions, and how representation-sensitive are those physical derivatives?
The earlier coupled experiment used original H, POLAR+GB and one peptide angle;
it supplied no paired correction. This is a new coordinate interface for the
existing normalized-H vacuum DFT/masked-OMOL hybrid, not a rerun of that model.

## Fixed inputs and coordinate rules

Use all eight prepared representations: GGR1GLG/2FW0/2FVY, both extended and
connected; alpha1F6S and6IP9 extended cores. Reuse only their actual center
native DFT and core/full masked analytic gradients. Two biological groups,
all consumed development; no new affinity label or calibrated decision.

Keep source atoms, protonation, normalized H, core membership, assembly,
explicit waters, caps and charge states unchanged. Reuse actual topology and
forcefield connectivity to define the following physical motions uniformly:

- Three Cartesian translations of the selected metal.
- Every side-chain chi bond on the CA-to-donor path of each declared Asp/Asn
  (CA-CB,CB-CG) or Glu/Gln (CA-CB,CB-CG,CG-CD) donor. Rotate the distal connected
  component, including its attached hydrogens, about that actual bonded axis.
  Reject rings, missing atoms or unsupported donor chemistry explicitly.
- A peptide crankshaft for every declared backbone-carbonyl donor: rigidly
  rotate its C/O and the actual bonded next N/H around the two source CA anchors.
  Use covalent connectivity, not residue_number+1. Require a nonproline peptide
  with exactly one attached N-H; reject unsupported chemistry.

Actual dataset therefore has16coordinates/GGR physical system and11/alpha.
No independent fragment/cap translations; water and exterior atoms remain fixed.
Both GGR representations use identical physical coordinates and measures.

Map source atoms directly and link caps through both retained/omitted anchors'
analytic Jacobians. Use the same physical Jacobian for DFT(core), T(core) and
T(full): g_H = J_core^T(g_DFT - g_Tcore) + J_full^T g_Tfull.
Record each component and g_R = g_Ca - g_La, with native units kcal/A for metal
translation and kcal/radian for angular motions. Do not combine their raw norms
or invent a coordinate covariance. For display only, also report the linearized
change for a fixed probe whose maximum physical-heavy displacement is0.02A
(the angular probe is0.02/max_atom_speed radians). These are first-order
predictions, not calculated energies, stiffnesses or affinity scores.

## Runs and checks

Zero new DFT, MACE, solver, fold, trajectory or optimization calls. Reuse16
matched electronic states. This is actual-data projection plus pure geometry
arithmetic. Verify source/core pairing and gradient atom ordering; finite-check
the coordinate Jacobians at1e-6A or1e-6radian (not numerical DFT gradients), with
absolute derivative tolerance1e-7A per coordinate unit. Check all physical bond
lengths to1e-9A under each probe and exact unchanged water/exterior inventories.
Report any failure, do not silently drop a motion. Source/core replay tolerance
5e-10A is inherited serialized-cap precision. Require common physical Jacobians
across the two GGR representations. No displaced scientific energy is implied.

Deliver reusable physical-coordinate/mapping code, real-fixture tests, source-
pinned gradients and report/vault note. Keep all response corrections unavailable
(`response_model_not_validated`). Choose any subsequent response potential,
optimization rule, trust limits and endpoint inventory in a separate frozen plan
before new scientific outputs. This projection does not itself improve accuracy
or qualify a mechanical model. Existing baseline and experimental records stay
unchanged; no push or default promotion.
