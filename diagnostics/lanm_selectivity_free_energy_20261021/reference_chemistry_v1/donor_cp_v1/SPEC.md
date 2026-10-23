# Matched dry donor counterpoise qualification

Question: does the0.4678kcal/mol PBE0-D4 versus canonical CCSD(T) discrepancy
for the paired water interaction cancel or change in balanced water-to-donor
exchange? This finite experiment extends the same frozen-f Hamiltonian, basis
and CC settings; no parameter fit, optimization, new geometry or force call.

## Fixed inventory

Use exactly the acetamide and bidentate-acetate central coordinates from
matched_water_removal_v1/prepared_v1/NATIVE_PANEL.json. These are the metal and
donor atoms retained unchanged after deleting their surrounding waters. They
are not silently replaced with the earlier published dry placements. Preserve
La/Dy paired coordinates and their physical1/6 multiplicity metadata; the
frozen-f calculations remain effective singlets.

Twelve canonical CCSD(T) cells:2donors x2metals x3roles (whole complex, real ion
with ghost ligand, real ligand with ghost ion). Eight matching electronic-only
PBE0 ghost cells, plus two real-ligand PBE0-D4 cells, shared between metals.
Reuse the four newly completed matched-dry PBE0-D4 complexes; never rerun them.
Retain the completed canonical/DFT water CP comparison as the water baseline.
Exactly22 new single points at most; no automatic expansion or retry.

| Donor | Complex atoms/AOs | Complex charge/electrons | Ion fragment | Ligand fragment |
|---|---|---|---|---|
| Acetamide |10/201|+3/40|+3/8|0/32|
| Bidentate acetate |8/189|+2/40|+3/8|−1/32|

Every CP role keeps the complete donor-complex basis. Real ligand-only cells
have9atoms/154AOs (acetamide) and7atoms/142AOs (acetate). Both real metals have
effective nuclear charge11; LaECP46 and DyECP55 remain unchanged. A ghost metal
uses DelECP, retains the exact metal GTO, and contributes zero nuclear charge,
electrons and ECP potential. Ghost ligand centers retain their bases. The
nominal Dy sextet is not modeled by this effective closed-shell Hamiltonian.

## Methods and scout

Canonical RHF/CCSD(T), NoFrozenCore, def2-TZVP light atoms and exact archived
metal GTO/ECP; no D4, RI/COSX, solvent or gradient. VeryTightSCF, PModel,
NoAutostart, SCF MaxIter300; MDCI MaxIter100/STol1e-8/PrintLevel3. NoFrozenCore
correlates all remaining40/8/32 electrons, including light-atom cores.

PBE0 ghost cells retain original RIJCOSX/def2-J/DefGrid3 and archived metal
auxiliary bases, PModel/VeryTightSCF/MaxIter300/NoAutostart, with no D4. The
two real ligands use the same PBE0 recipe plus D4, with their actual0/−1 charge.
Dispersion is evaluated only on real molecules; isolated-ion D4 is the exact
zero empty-pair/triple sum. Free-ligand D4 is subtracted in each interaction,
even though it cancels in a paired metal contrast.

All native cells use4MPI/one thread each, standard_expand/default memory/time,
MaxCore=floor(75% actual memory grant/ranks). First submit only the largest
La acetamide complex, cc__acetamide__La__complex, to establish actual canonical
cost and admission at201AOs/40electrons. Root decides remaining release from
the measured result. No molecular preflight repeats already qualified water.

## Admission and balanced contrast

Generalize the qualified water parser to the actual number of centers; check
every basis assignment, total AO count, coordinates, effective nuclear charges,
ECP assignment, charge, spin and electron inventory. CC must report RHF,
NoFrozenCore, all expected correlated electrons, native SCF/CC convergence,
STol and final residual evidence, and CC/final energy closure within1e-8Eh.
DFT must show PBE/PBE with25% exact exchange, RIJCOSX, actual dispersion scope
and SCF/final/D4 additive closure. Preserve all raw outputs and SCF stopping
signals; do not reinterpret inactive printed criteria as an undeclared failure.

For each method, donor L and metal M:

    I(M,L) = E(ML) - E(M|ghost L) - E(L|ghost M)
    Exchange(M,L) = I(M,L) - I(M,water)
    Paired(L) = Exchange(Dy,L) - Exchange(La,L)

This is the CP-corrected fixed electronic interaction form of
M-water + L -> M-L + water. Charges and elements balance for neutral amide
and anionic acetate. PBE0-D4 adds only the matched real-molecule dispersion
interaction to its electronic CP value; CC receives no added D4. Report
electronic-only and D4-added results separately, all individual terms and the
metal-dependent donor contrast relative to the retained water baseline.

Before outcomes, use the same0.25kcal/mol exploratory resolution for the primary
paired PBE0-D4-versus-CC exchange discrepancy, separately for each donor. Also
report whether its signed difference is reduced or increased relative to the
water-only discrepancy. Do not select a Hamiltonian or threshold after results.

The acetate anion's missing diffuse-basis test, remaining semicore basis,
metal-basis convergence, frozen-f/spin-orbit physics and numerical precision
remain separate uncertainties. This same-basis comparison is not correlated
physical truth, force qualification, hydration thermodynamics or selectivity.
