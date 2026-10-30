# Bounded native ALMO quartet

## Question and release

At these four archived fixed geometries, does suppressing metal–whole-shell
delocalization materially change the La/Dy water-to-acetate double difference?
The active user goal authorizes one La-aquo compatibility pilot, followed only
on admission by the other three systems concurrently. All previous calculations,
software attempts and other users' jobs remain preserved. This is not a fitting
or protein campaign.

`DECK_PROVENANCE.json` binds the four original XYZ hashes, exact translated
coordinate tokens, electronic states and deck hashes. Each deck has two groups:
M3+ with eight explicit electrons, and the entire shell with 80 (aquo) or 102
(mixed) electrons. All fragments and complexes are effective singlets. Dy's
physical sextet metadata is not the runtime spin state.

## Single electronic protocol

The pinned isolated CP2K2026.2 build implements native molecular-guess ALMO,
block-diagonal polarization, then fully delocalized `FULL_SCF`. The only source
change is the separately reviewed seven-line Dy55 input-reader repair.
PBE0, literal archived La46/Dy55 semilocal ECPs and metal orbital contractions
are preserved. H/C/O retain def2-TZVP with the documented named-library
equivalence limitation in `CONVERSION.md`. This is not a claim of bitwise
equivalence to ORCA's internal library or RIJCOSX implementation.

All four decks use the same prospective numerical realization: 30 Å isolated
MT cell without centering, GAPW cutoff1200 Ry/relative100 Ry, five grids,
200 radial/590 angular points, exact full-range25% HF, explicit tight Gaussian
screening, fragment SCF threshold1e-10, block/full ALMO thresholds1e-9 and
300 iterations each. See the decks and `PILOT_SETTINGS_REVIEW.md` for every
control. No ADMM, truncated exchange, alternate potential, solvent, optimization,
wavefunction reuse, smearing or protocol sweep is authorized.

CP2K performs the electronic chain only. The exact archived per-system D4 scalar
is added once in the analysis to frozen, polarized and full energies equally.
It leaves POL and CT unchanged. Native isolated-fragment electronic references
are retained. Unavailable absolute whole-shell D4 references are not invented;
they cancel between matched metal pairs in the quartet double difference.

## Prospective admission gates

Before release of the remaining three, root reviews actual pilot output:

1. Exact input/binary/assets, correct atom/AO/electron counts, two fragments and
   requested effective states. Fresh fragment, block and full SCFs must converge
   at their declared native thresholds; native process exits successfully.
2. All R/F/P/U energies are available: native fragment sum R, frozen initial F,
   final converged block P, final fully delocalized U. Preserve the raw trace.
3. `abs(U - archived electronic energy) <= 1e-4 Eh`, equivalent to0.0627509474
   kcal/mol. Adding the same archived D4 must give the same full-total agreement.
   This selected cross-code gate is not an estimated uncertainty bound. Archived
   individual density/residual flags and final COSX shifts remain explicit.
4. Telescoping arithmetic closure <=1e-8 Eh, with actual printed-digit bounds
   reported separately; P<=F and U<=P within1e-8 Eh. These checks cannot substitute
   for native convergence or reference agreement.
5. No unresolved fragment/state/convention discrepancy that contaminates the
   double difference. Missing output is unavailable evidence, not assumed success.

All four cells must meet these gates for an admitted quartet. The gate's worst
case sum over four cells is0.2510037896 kcal/mol; report the actual full-energy
double-difference disagreement rather than treating this as a universal error
bar on CT. No failed gate triggers a silent tolerance change or numerical retry.

## Execution and finite cost

One shared Slurm node per system: banfield/standard, `--qos=standard_expand`,
24 allocated CPUs,12 native image-MPICH ranks ×2 OpenMP threads. BLAS uses one
thread. Runtime and libraries come from the pinned image; no host MPI mixing.
Scheduler memory default is4000 MB/CPU (96,000 MB for24 CPUs). No memory or
walltime override is added. HFX integral storage512 MB/rank is included within
that grant; peak actual memory is recorded when available.

Each native parent call has a **7200-second process cap**, no restart/extension.
This is a declared development budget for internal fragment/block/full SCFs;
archived four-rank ORCA parent timings were374–653 seconds including gradients.
The cap is not a CP2K performance prediction. At most691,200 allocated CPU-seconds
are budgeted for the four native calls, plus short startup/analysis overhead.
The runner enforces this limit only on its own new process group, without
cancelling any existing Slurm job. Syntax validation precedes the molecular call.
All calls retain raw output, launcher logs, source pins and receipts. There is
no automatic retry, automatic remaining release or repeat of completed work.

## Analysis and stopping rule

Compute FRZ=F−R, POL=P−F and CT=U−P. Native FULL_SCF supports combined FRZ;
an electrostatic/Pauli split is unavailable in this implementation path and
must be marked unavailable rather than reconstructed from unrelated terms.
For each quantity X, use C_X=(X_Dy,mixed−X_La,mixed)−(X_Dy,aquo−X_La,aquo).
Negative means relative Dy stabilization on substitution. Retain numerical
fragment-reference cancellation residuals; never force them to zero.

The final report compares full and CT-suppressed contrasts, all component
contrasts, archived D4 and reference terms, the original model discrepancy
(−6.810868 kcal/mol), actual numerical residuals and uncertainty limitations.
A large CT contrast establishes a contribution in this declared electronic
comparison, not the identity of a force-field defect or a binding free energy.
On pilot incompatibility, report the exact blocker and stop this quartet.
On completion, give one evidence-supported recommendation and stop.
