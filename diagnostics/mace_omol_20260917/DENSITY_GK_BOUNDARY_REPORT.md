# Full physical AMOEBA/GK source boundary: preparation checks pass

All16native initializations and12paired environment/cavity checks pass. These
are four real representations times Ca/La source and zero-source reference.
No energy, Born calculation, induction solve, DFT or MACE ran during preparation.

The full physical systems retain4698GGR,1932alpha1F6S and1898alpha6IP9 atoms,
including their real source metal, declared assembly, protonation and0/2/3waters.
No cap nucleus enters the physical cavity. Existing source-projected quantum
charges provide the reaction-field proxy; printed charge residuals are retained.

AMOEBA residue sums independently reproduce their chemical formal charges.
Remove force-field moments on source-support atoms and restore each intersected
residue's exterior formal charge through its distinct same-residue exterior
neighbors directly bonded to source support. Retain complete per-residue charge
and bond ledgers. No global neutralization, endpoint-specific redistribution,
duplicate caps or copied ff19SB increments. Exterior dipoles/quadrupoles stay
unchanged. Exterior charges close to -3e forGGR and -4e foralpha; source ligand
formal charge is -3e, with Ca-source -1e and La-source0e.

Both endpoints use the explicit Ca2018-derived source-metal cavity specification:
native GK radius1.82485A, descreen radius1.795A, source polarizability0, internal
dielectric1, native bulk78.3, Grycuk neck/tanh andGKC2.455. This is a common
boundary convention, **not an AMOEBA La force-field parameterization**. Actual
zero-source Ca/La preparations prove environment and cavity identity. Native
`kpolar` initially removes the neutral nonresponsive metal from its multipole
index; the isolated frontend restores all physical source indices after moment
replacement. It introduces no artificial initialization charge and changes no
native library routine.

## Actual cost and tests

Preparation65.3894131966s wall/47.460063526processCPU;16native initializations
sum16.2488231175s wall/16.064919CPU. Build0.5353315s wall/0.498959CPU.
No GPU or cluster allocation for initialization-only work. Three real-artifact
tests pass39.789s, no skips: residue ledger/explicitly corrupted formal charge,
all16native states with retained source index, paired cavity/reference identity
and an explicitly corrupted radius copy.

Detailed compact record: [DENSITY_GK_BOUNDARY_RESULT.json](DENSITY_GK_BOUNDARY_RESULT.json).
Protocol `density_direct_common_Ca2018_GK_source_boundary_v1`.
Manifestdf84a6612ba2df87f822f854289cdbf08fb43929726056018094bc7d30c715fe.
Products: `workspaces/mace_omol_20260917/density_gk_boundary_v1`.
This qualifies preparation only; the separate full hybrid pilot must establish
numerical stability and predictive usefulness. Production baseline unchanged.
