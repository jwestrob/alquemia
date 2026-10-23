# Matched PBE0 water counterpoise comparator

Question: does the paired La/Dy water interaction differ between canonical
CCSD(T) and the original PBE0 model after matching the full counterpoise cycle?
Keep the original geometries, frozen-f ECPs, orbital bases, effective singlets,
PBE0/RIJCOSX/DefGrid3/VeryTightSCF recipe and archived auxiliary metal bases.
This is a five-call diagnostic; root owns submissions. No new geometry or force
calculation, electronic-state assertion, or method-accuracy admission follows.

## Exact inventory and reuse

Reuse both admitted original complex outputs. Their printed TOTAL SCF ENERGY
is the electronic PBE0 energy; verify it plus printed D4 equals the final energy
within 1e-8 Eh. Do not use final PBE0-D4 energy as an electronic-only quantity.

Four new no-D4 single points: La/Dy real ion with ghost water, and real water
with ghost ion. Use the exact four-center canonical counterpoise decks' basis,
coordinates, charge, effective state and ghost mask. Each has 90 orbital basis
functions. Ion roles have charge+3/electrons8; water roles charge0/electrons10.
Real metal keeps ECP46/55; ghost metal has DelECP and the exact retained GTO.
Keep all auxiliary basis centers, including metal auxiliary functions on ghosts.
No D4 term is requested or admitted for a ghost calculation.

One separately declared real-water PBE0-D4 single point supplies its dispersion
correction, with the exact common three ligand atoms and no ghost metal. This
is neutral, effective singlet, 10 electrons and 43 orbital basis functions.
La/Dy water coordinates must match exactly before reuse. No compatible cached
real-water DFT/D4 result was found in the active 21October reference workspace.
The isolated real ion's D4 dispersion is exactly zero because it has no atom
pairs or triples; record this algebraic empty-sum term, not a calculated result.

For each metal:

    electronic_CP = SCF(complex) - SCF(ion|ghost water) - SCF(water|ghost ion)
    dispersion_interaction = D4(real complex) - 0 - D4(real water)
    PBE0_D4_CP = electronic_CP + dispersion_interaction

Report both electronic-only and D4-added values and their Dy-minus-La contrasts.
The common real-water D4 cancels in the paired contrast, but is retained in each
individual interaction. Compare either matched CP Hamiltonian explicitly with
canonical CCSD(T), which receives no added D4. PBE0 RIJCOSX versus conventional
CC also retains an integral-approximation difference; no unique-correlation or
physical truth claim follows. Do not compare original uncorrected complex gaps
as counterpoise interaction energies.

## Runtime and admission

ORCA6.1.1, four MPI ranks, one thread/rank, standard_expand, scheduler default
memory/time; MaxCore=floor(75% actual grant/ranks). PModel, NoAutostart, SCF
MaxIter300. No automatic fallback, solver retry, optimization or gradient.
Native output must admit normal termination and SCF convergence, closed-shell
PBE0 with 25% exact exchange, correct charge/electrons/multiplicity, all orbital
basis functions/centers and their native effective nuclear charges, exact ECP
placement, finite electronic/final energies and additive D4 closure. Native
settings must show RIJ-COSX. Ghost outputs must have no D4 correction or module.
Retain all printed SCF residuals and raw output even when native admission passes.
The existing protocol's SCF stopping behavior remains an uncertainty, not a
reason to silently tighten/relax this matched comparison after seeing results.

First inspect actual ghost-ion and ghost-water native output admission. The
separate real-water cell must prove its three-center/no-ECP/43-function inventory.
Collection requires all five new results, both pinned saved complexes and the
complete six-cell canonical collection; partial comparisons are not final.
