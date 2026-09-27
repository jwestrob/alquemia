# Hans8DQ2 EF3: explicit-state local response scout

Declared 27 September 2026 before any new LanM reference energy. This is stage5
of Jacob's 26September handoff, after the completed three-core Ca/La response
experiment. It does not restart whole-protein vacuum subtraction or relaxation.

## Question and gate

Can the installed native r2SCAN-3c model execute the physically declared LaIII
singlet and DyIII sextet on a chemically complete EF3 region, and what restoring
response does it predict for one small real donor reorientation? This probes a
surface needed for future specificity modeling, not protein-level affinity.
The paper-matched MACEPOL-EF model remains unavailable; no substitute is called
that candidate. Existing MACE results on other regions are not matched references.

Start with consumed Hans8DQ2 only, chosen as the handoff's first source before
new results. Four endpoints: La/A, La/B, Dy/A, Dy/B. Measure this first scout
before extending to Hans8FNR and Mex8FNS. The latter are not optional favorable
source selectors: any subsequent source comparison must retain all three and
all missing/failed cells. No reserved numerical outcomes are opened.

## Fixed preparation

Use the audited whole-chain source's already normalized hydrogen coordinates,
not old unprojected donor XYZs. Preserve original source occupancy, including
actual EF4 sodium. Target EF3 alone changes LaIII↔DyIII; spectator EF1/EF2 LaIII
and EF4 NaI remain frozen formal monopoles (+3,+3,+1) in both endpoints.
This is an explicit electrostatic approximation to neighboring occupancy,
not a spectator force field or treatment of intersite quantum spin coupling.

The declared region in lanm_preparation_feasibility/RESULT_v2.json contains
195 electronic atoms including four mapped caps: complete EF3 loop83–94,
flanking peptide units and nonlocal65–66 hydrogen-bond support, with the
source-supported selected waters. Keep all region identities/protons/waters
fixed. Record every retained/cut bond and cap Jacobian. Heavy geometry is not
optimized or repaired to match an expected preference.

A is the source. B rotates the actual Asp85 side-chain CG/OD1/OD2 and HB2/HB3
by +2degrees about CA→CB using the declared atom identities. The motion preserves
bond lengths and angles at CB; require maximum physical displacement0.15Å
and unchanged covalent graph before admission. It represents a modeled alternative, not an
experimental structure or population. Other physical atoms remain fixed.

MM protein/water charges come from the exact existing topology with ff19SB/TIP3P,
transferred by atom identity to normalized coordinates. Do not invoke hydrogen
addition or metal parameterization. Selected real QM atoms have no MM charges.
Zero cap-omitted CA charges and restore each affected partial residue's native
formal charge minus its selected neutral peptide fragment on that residue's
actual bonded exterior heavy neighbors. Record every amount and recipient;
this local partial-residue charge restoration is NOT native ORCA CS and is not
dipole-preserving. No global neutralization. For Hans8DQ2: QM−1, MM+7, total+6,
including the unchanged spectator monopoles. Preparation must verify these values
and inventory before submission. Preserve a common finite potential gauge
(zero at infinity) and identical fields for paired metals.

## Electronic treatment and evidence

Use native `r2SCAN-3c NoAutostart DefGrid3 TightSCF EnGrad`, finite point charges,
DoEQfalse, no added CPCM or duplicate MM Coulomb. Installed def2-mTZVPP export
confirms native LaECP46 and DyECP28: Dy4f remains explicit. No basis/ECP substitution,
spin-independent xTB, effective-singlet shortcut, numerical gradient or SOC claim.
Declare La multiplicity1 and Dy multiplicity6 with closed-shell ligand hypothesis;
verify actual charge/electron counts, ECP, RHF/UHF and analytic SCF/ECP/D4/gCP/PC
components from outputs. Retain S², per-atom spin populations and any actual
stability evidence; missing stability/alternate-solution evidence stays missing.
Physical sextet is not proof of the ground state or a complete spin-orbit model.

Primary output: W_M=E_M(B)−E_M(A), and W_Dy−W_La. Report the individual works,
actual gradients and projection onto the exact source-derived torsion (including
mapped-cap chain rule wherever relevant). No raw cross-element energy or Ca/La
band becomes affinity. Both configurations have identical stoichiometry.

The earlier DFT rigid-transform test failed. Do not promote this recipe to
qualified hybrid forces through this scout. The existing0.05kcal/mol numerical
response target is an intended resolution, not a proven error bar; contrasts
at or below that scale are unresolved. No blind validation claim, new fitted
threshold, source averaging or adaptive sign rescue. A failed/unsupported Dy
state stops dependent expansion; inspect evidence rather than retrying unchanged.

## Execution

Reuse the manifested native runner, allocation-derived memory renderer, strict
explicit-state parser and immutable receipts. Four concurrent independent native
endpoints share the actual full-node MPI slots and RAM, --mem=0, normal priority,
no dependency on another session's PQQ jobs and no test queue. Select a currently
available supported node/layout and record it before submission. No fixed project
compute/time cap, no blanket panel: these four tasks define the initial work.
Record failed attempts and actual allocatedCPU time/RSS. New high-level work is
explicitly4endpoint evaluations; no GPU, MD, optimizer or quantum Hessian.
Actual completion wake must be armed. Preserve restart notes ahead of Monday
28September10:00America/Los_Angeles shutdown. No molecular stage starts after
that cutoff. No production/default/push/email action follows automatically.
