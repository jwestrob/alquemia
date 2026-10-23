# One bounded electronic-solver recovery for the Dy–water reference

Question: can the exact failed17260 all-electron Dy–water endpoint converge when
the trust-region solver is retained through the observed late oscillations?
The saved assessment `ae_water_recovery_assessment_v1/REPORT.md` (9cb284b3)
and actual native output establish electronic nonconvergence, not a reporting
defect. No old energy or gradient is admitted. This recovery is independent of
the frozen-water orientation test and does not alter its model or reference.

Run exactly one fresh four-atom EnGrad endpoint. Preserve the original geometry,
charge+3, physical/effective sextet, unrestricted PBE0-D4, scalar ZORA(MP),
SARC-ZORA-TZVP Dy and ZORA-def2-TZVP water bases, SARC/J, RIJCOSX, DefGrid3,
PModel guess, NoAutostart, VeryTightSCF and MaxIter300. Add only the nested
`trah / QuadRegionStart 1e-6 / end` setting inside the original SCF block.
The documented NR-entry control is separate from unchanged TolG2e-6. This tests
a solver hypothesis; it does not change the Hamiltonian, spin sector or accuracy
thresholds. No unconverged-wavefunction restart, smearing, extra state, automatic
iteration extension or fallback is included.

Singer owns a small versioned adapter/preparation proving the exact sole deck
delta and preserving native source/input/output provenance. Root owns review and
submission; Arendt reviews source and preparation independently. Original failed
inputs/results stay immutable. Any reused native helper remains pinned; recovery
metadata must explicitly identify the wrapper and changed solver setting, never
label this as an unchanged baseline input. Preserve raw native result separately
from the additional recovery admission result.

One four-MPI Slurm job, one thread/rank, standard_expand and scheduler default
memory/time. Derive MaxCore from the actual grant under the existing75% rule.
Native startup must print the intended NR-entry threshold when TRAH activates.
If SCF converges before TRAH is used, report that and do not claim the new control
caused the success. Retain native SCF criteria, convergence/normal ending,
energy–gradient closure, basis/electron/state checks and finite gradient admission.
Inspect final orbital residual against the actual active threshold. Preserve spin
expectation, occupation and population output when available; unavailable state
analysis cannot become a physical-state pass.

A converged result would fill one missing reference endpoint. A different local
determinant remains possible, and frozen-f versus explicit-f comparisons also
change basis/relativistic treatment. This one endpoint cannot establish global
electronic state, spin–orbit cancellation, force accuracy or protein selectivity.
The balanced comparison still needs its other admitted endpoints. Failed recovery
remains failed, with no automatic second attempt.

Preparation/review make zero native calls and use at most one local CPU each.
Freeze actual configuration/release before the single molecular submission;
record completion/health watchers, actual cost and independent admission.
