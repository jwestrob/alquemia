# Six-cell canonical water counterpoise pilot

Version2 preserves the original version1 source and prepared bundle. Scientific
input decks are unchanged. Parser hardening isolates orbital from auxiliary
basis summaries, requires the actual CC-iteration convergence marker and rejects
explicit CCSD nonconvergence, and verifies native NO frozen core treatment.

Question: what is the paired La/Dy water interaction at canonical CCSD(T), with
the initial panel's frozen-f ECP and exact orbital basis? This is the first
bounded tier in CORRELATED_QUALIFICATION_PLAN.md. No ligand exchange or affinity
claim follows from this water-only pilot.

## Inventory and method

Exactly six cells: La/Dy times complex, real ion with ghost water, and real water
with ghost ion. Coordinates are the exact original central water dimers. Every
role retains four centers and90 spherical orbital basis functions: metal47,
oxygen31, hydrogen6 each. Complex charge+3/electrons18; ion fragment+3/electrons8;
water fragment0/electrons10. All are effective singlets. Physical Dy sextet is
retained as provenance; the removed4f shell remains outside this Hamiltonian.

Use ORCA6.1.1 conventional RHF/CCSD(T), def2-TZVP on O/H and the exact archived
metal NewGTO coefficients. Use LaECP46/DyECP55 only on a real metal. NoFrozenCore
correlates every electron remaining outside the ECP. No DFT, D4, RI/COSX,
continuum solvent, optimization or gradient is requested. VeryTightSCF,
PModel, NoAutostart, SCF MaxIter300; MDCI MaxIter100, STol1e-8, PrintLevel3.
The CC residual target is declared before any result and is not fitted to data.

Ghost atoms use ORCA's colon notation. For a ghost metal remove the entire
NewECP block and specify DelECP on that element while keeping its exact NewGTO
block. A real metal retains the complete archived ECP/basis block. This prevents
reliance on an assumed implicit ghost-ECP conversion. Actual output admission
still checks zero ghost nuclear charges and no ghost ECP assignment.

Four MPI ranks, one thread per rank, standard_expand, scheduler default memory
and time. MaxCore=floor(0.75×actual memory grant/ranks). No repeated or fallback
calculation. A technical failure preserves output and requires a separately
recorded recovery decision.

## Admission and evidence

Static preflight renders all six decks and verifies input pins, fragment charge,
electron bookkeeping, exactly4centers, ghost masks, unchanged geometry and
unchanged metal orbital coefficients. It makes zero native calls. Native input
acceptance and actual ghost behavior are not proven by static preflight.

Native output must show normal termination, SCF and coupled-cluster convergence,
RHF, multiplicity1, the expected charge, total and correlated electron counts,
all90 basis functions and the expected per-center contracted basis counts.
The atomic-unit coordinate table must contain all four centers, correct symbols
and coordinates within1e-6Å (printed precision), and effective nuclear charges
exactly [11,8,1,1], [11,0,0,0], or [0,8,1,1] for the three respective roles.
ECP assignment must occur only on real center0 and have the expected46/55 core;
the water-with-ghost-ion role must have no ECP assignment or replacement count.
The CCSD(T) total must agree with FINAL SINGLE POINT ENERGY within1e-8Eh.

Save actual output, stderr, input, orbital file metadata, wall time and allocated
CPU time. Extract numeric printed SCF residuals and retain native correlation,
amplitude, triples and diagnostic lines. Missing diagnostic fields stay missing;
T1 alone is not a multi-reference or accuracy certificate. Numerical admission
does not establish basis-limit, physical-state, force or method accuracy.

## Collection

Require the complete six-cell denominator. For each metal:
Einteraction_CP=Ecomplex−Eion_with_ghost_water−Ewater_with_ghost_ion.
Then report Einteraction_CP(Dy)−Einteraction_CP(La), with all term pins. This is
a fixed water interaction contrast, not the water-to-donor exchange contrast.
The corresponding ligand CP calculations, basis/force tests, hydrated cluster
tests and electronic-state/SOC uncertainties remain separate future work.

Official implementation references: [CCSD(T), frozen-core and diagnostics](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/mdci.html),
[ghost-center counterpoise](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/counterpoise.html),
[native correlated-electron output example](https://www.faccts.de/docs/orca/5.0/tutorials/prop/corren.html).
