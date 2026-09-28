# One repaired embedded Dy origin with explicit TRAH

Declared before evaluation. The original embedded four-cell job1220312 has spent
about an hour with Dy oscillating in SOSCF: recent maximum orbital gradients
~4e-4–5e-3, substantial reversals and no accepted solution. La residuals improve,
so preserve that job while a distinct solver capability scout is tested. Do not
extend iteration limits, accept unfinished energies or restart the old explicit-f
campaign. Save actual residual evidence and account both attempts.

Use exactly one Dy_A from reviewed C-bound-H repaired Hans8DQ2 full-regionv2.
Keep195atoms, fixed inventory/protonation/occupancy, charge−1,806explicit electrons,
physicalDyIIIsextet represented by ECP55restrictedvalence1. Same PBE0-D4,
light def2TZVP/def2J, literal Dy lcecp1TZ/AuxJ, RIJCOSX,DefGrid3,VeryTightSCF,
NoAutostart,PModel,DoEQfalse and finite repaired permanent field. Add only the
explicit TRAH solver keyword to this repaired-source Hamiltonian. No old energy
can satisfy this new-coordinate task; no old wavefunction guess copied live.

This changes hydrogen coordinates/field relative to the ongoing old matrix AND
uses a different SCF algorithm. Therefore success cannot establish which change
resolved convergence. The purpose is to obtain a usable physically improved
reference origin, not a causal solver-vs-H experiment or favorable metal score.
No A/B or La/Dy contrast exists from this single endpoint. Native gradients and
pointcharge gradients must pass existing state/coordinate/output checks; require
actual TRAH/NR macro-iteration evidence, not just an echoed keyword.

ORCA6.1 documents TRAH for restricted/unrestricted Hartree–Fock/Kohn–Sham and
RIJCOSX. Its robust orbital optimization is distinct from nuclear Hessian work;
no numerical DFT forces or nuclear Hessian is requested. A converged local orbital
minimum is not a global electronic-state or biological validation.
Reference: https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html

One endpoint on an idle112CPU standard node,112MPI,mem0, full allocation-aware
memory policy, normal schedulerpriority, no GPU. Fresh terminal/health watches.
Costs are unknown before this scout and will be measured; do not extrapolate a
speedup or automatically expand to the remaining sources. If it converges,
collect real forces before deciding the finite repaired A/B continuation. If it
fails/stagnates, diagnose actual residuals instead of increasing limits blindly.
Original and repaired matrices, source-H defects and full-hybrid/solvent limitations
remain explicit; no production/default change.
