# Bounded Tauri screen — 6 November 2026

**Decision: unresolved reference; no solvent pilot. End automatic expansion.**

1. **Improved selectivity prediction was not demonstrated.** L/M inference and
   its numerical checks succeeded, but no balanced contrast has a qualified
   physical-spin reference. The historical donor-count/charge M1 remains the
   baseline; its [earlier 0.854-kcal/mol family-held-out result](../ln_structure_increment_20261030/REPORT.md) was not rerun here.
2. **The potential shows material sensitivity.** L and M disagree by several
   kcal/mol on the balanced contrasts. A neutral water 100 Å away changes L's
   combined-system energy by −4.25 kcal/mol beyond the isolated-fragment sum and
   changes a base-atom force component by 1.75 kcal/mol/Å. These observations
   require attention before assuming cancellation in a periodic identity swap.
3. **The quantum screen is unresolved, rather than a Tauri success or failure.**
   Three La references finished, but each misses the frozen Tight density/DIIS
   criteria. The three unfinished Dy calls were stopped once that prerequisite
   had failed. No Dy final energy, gradient or spin qualification is inferred.
4. **Neither experimental development pair was admitted.** The eight native
   records match the stated within-ligand conditions, but shared electrolyte and
   original-study provenance remain unresolved.
5. **The next action is to stop before MD and preserve M1.** Further simulation
   is not released by this diagnostic. A future continuation needs a separately
   justified reference qualification and a verified experimental medium.

## Frozen experiment and implementation

[PLAN.md](PLAN.md) was frozen before inference. The checkout was f1181b4; newer
uncommitted work and other-session executors were preserved. The complete panel
contains four W8/W9 aquo shells, two dry monoacetate centers, sixteen archived
radial/angular displacements, and two hydrated AcW7 configurations. All positions
and atom orders were copied from their archived records. No structure was
optimized. All twelve La/Dy coordinate pairs match exactly apart from the metal.
These are existing development fixtures, not prospective selectivity validation.

Source: [SO3LR 70707f5](https://github.com/general-molecular-simulations/so3lr/tree/70707f5bb0c333ecec707e5a12b1efdf1c83603a),
GLP2294f4f, L/M-2606 weights, separate Python3.12.0/JAX0.10.2 environment.
[SOFTWARE_PINS.json](SOFTWARE_PINS.json), [environment lock](environment-lock.txt),
[input pins](INPUT_PINS.json) and [artifact inventory](ARTIFACTS.json) retain
source, weights, inputs and raw-output hashes. Only the required two models were
retrieved; no potential code or embeddings were fitted or changed.

All numerical work used Slurm `standard_expand`. Static inference used twelve
allocated CPUs because H200 scout49306 was blocked by QOSMaxGRESPerUser and
cancelled before allocation. Physical state inputs were explicit: aquo Q=+3,
acetate Q=+2; La multiplicity1 and Dy6. The native loader's charge, unpaired0/5,
head5 and exact coordinates were asserted. Standard native inference used float64,
x64 enabled before loading, highest matmul precision, no PBC, SR5 Å and LR1000 Å.
Raw device arrays were saved before formatting. The no-cutoff control uses
`lr_cutoff=None` and its required `dispersion...damping=None`.

Both models pass the six checks in [TAURI_RESULTS.json](TAURI_RESULTS.json):
intact LaW8 ASE/native agreement, an asserted metal(+3)/acetate(−1) translated
fragment graph, and a real-coordinate collective-shell finite difference.
Maximum ASE/native discrepancies are below 8e-14 eV and 2e-14 eV/Å. Shell
finite-difference residuals are below 5e-7 eV/Å at the frozen 1e-4-Å step.
No field argument, artificial field-energy response or Flash fragment substitute
was used. Fragment charge constraints were confined to the translated diagnostic.

## Fixed-coordinate results

These are **total electronic configuration contrasts**, in kcal/mol:

| Model and LR cutoff | C_W | C_Ac |
|---|---:|---:|
| Primary L, 1000 Å | −1.113 | +3.963 |
| L, 100 Å | −1.095 | +4.260 |
| L, no cutoff | −1.115 | +3.930 |
| Sensitivity M, 1000 Å | +2.560 | +10.769 |
| M, 100 Å | +2.718 | +10.810 |
| M, no cutoff | +2.542 | +10.765 |

C_W=(DyW9−DyW8)−(LaW9−LaW8);
C_Ac=(DyAcW7−DyW8)−(LaAcW7−LaW8).
W8/W9 also differ in geometry. Neither contrast is a solution free energy or a
La/Dy ligand affinity. Elemental offsets cancel algebraically; global charge/spin
conditioning and size consistency do not thereby become validated. M is the
predeclared capacity sensitivity, not a replacement selected for agreement.

[Force readouts](force_readouts.tsv) contain individual oxygen radial projections,
metal projections and rigid-water collective shell derivatives. Increasing shell
R translates every OHH triplet outward along its metal–O direction. The acetate
translation/rotation axes are the original geometry-defined axes; rotation uses
q=theta×1 Å. AcW7's seven water triplets and the unchanged acetate atom group are
reported separately. Oxygen projections do not assert that every oxygen binds.
[Dry curves](dry_curves.tsv) retain both original ±0.005/±0.010-Å steps.
[Exact-coordinate delta forces](../../workspaces/so3lr_tauri_screen_20261106/DELTA_FORCES.json)
report Dy6 minus La1 on all twelve pairs for each model.

Selected collective shell dE/dR values, kcal/mol/Å:

| Fixture | L | M | New native La reference | Qualification |
|---|---:|---:|---:|---|
| LaW8 | −117.13 | −119.75 | −125.98 | unresolved SCF residuals |
| DyW8 | +45.74 | +70.92 | unavailable | stopped before final gradient |
| LaW9 | −260.80 | −248.80 | −258.32 | unresolved SCF residuals |
| DyW9 | −79.12 | −51.34 | unavailable | stopped before final gradient |
| LaAcW7 waters | −112.92 | −98.76 | −115.28 | unresolved SCF residuals |

The raw La slopes resemble the native results descriptively. That observation
does not qualify the reference or establish improved accuracy. All complete
comparisons and qualification flags are in
[reference_force_comparison.tsv](reference_force_comparison.tsv).
The [historical force comparison](historical_force_comparison.tsv) is explicitly
descriptive: its Dy uses effective-singlet ECP55, whereas Tauri uses physical Dy6.
No CP, correlated-fragment or old force-field operator replaced the original
PBE0-D4 targets.

The spectator energy excess is −4.29/−4.25 kcal/mol for L at50/100 Å and
−3.68/−3.61 for M. Its persistence and corresponding force changes show sensitivity
to the combined-system representation; they do not isolate a unique mechanism.
Explicit fragment charge treatment changes the translated dry frame by
18.12/17.71 kcal/mol for L/M relative to plain conditioning. No correction was
fitted to remove these effects from intact systems.

## Reference bridge and stopping decision

The [old raw-output audit](../../workspaces/so3lr_tauri_screen_20261106/OLD_REFERENCE_AUDIT.json)
reproduces DyW9→DyW8+water as +18.51456 kcal/mol with ECP55/effective singlet and
−9.03750 with ECP28/physical sextet: a 27.55206-kcal/mol discrepancy. The latter
W8/W9 calculations took5859/513 SCF cycles, with S²=9.272561/9.244627. Their
printed density/DIIS residuals fail, and no orbital stability qualification was
available. Their optimized geometries differ from this panel, and they provide
no matched gradients. The discrepancy is reproduced, not established as physics.

Exactly six fresh energy/gradient invocations49309–49314 were started under the
frozen plan. The pinned [FAIRchem OMol recipe](https://github.com/FAIR-Chem/fairchem/blob/d4e10d0f4e87593922fda1d1ffb7dd28e504cadd/src/fairchem/data/omol/orca/calc.py)
and external def2-TZVPD basis were used with wB97M-V, RIJCOSX, def2/J and DEFGRID3.
This is an ORCA6.1.1 port of the6.0 recipe, with predeclared PModel initialization
and non-restarting stability diagnostics; La retains the20-degree beta rotation.
Actual native output confirms La ECP46, Dy ECP28, UKS and expected electron
counts88/98/110 and115/125/137. No D4 or SCNL keyword was appended. EnGrad
automatically selected self-consistent NL and disabled AutoTRAH; it is not an
exact mixed-training-Hamiltonian identity claim.

| Completed La reference | MAX density (tol1e-7) | RMS density (tol5e-9) | DIIS (tol5e-7) | S² | Lowest stability eigenvalue, Eh |
|---|---:|---:|---:|---:|---:|
| W8 | 3.3242e-6 | 2.5060e-7 | 2.1343e-6 | 0 | +0.312950 |
| W9 | 6.1388e-6 | 2.8757e-7 | 3.5111e-6 | 0 | +0.298345 |
| AcW7 | 1.5973e-5 | 5.4155e-7 | 7.9191e-6 | 0 | +0.155629 |

Each terminated normally after energy-based SCF convergence, with finite matched
gradients and stable wavefunctions as tested, but failed three mandatory printed
Tight residuals. Both balanced contrasts require LaW8. Therefore, finishing the
remaining Dy jobs could not produce a qualified contrast under this protocol.
[SCIENCE_STOP.json](SCIENCE_STOP.json) records the decision before cancelling
only49311/49312/49314, after2686 seconds each. Their partial outputs remain
preserved; they were not run to a convergence failure or assigned final energies.
No restart, tighter solver, new reference method or changed acceptance criterion
followed. The provisional1-kcal/mol contrast budget remains **unassessed**.

[REFERENCE_RESULTS.json](REFERENCE_RESULTS.json) retains method/state checks,
occupations, S², stability, residuals, gradients and complete native file hashes.
Coordinates are checked against the twelve-decimal native property file and
the lower-precision engrad export using
[ORCA's documented units](https://www.faccts.de/docs/orca/6.0/manual/contents/general.html#units-and-conversion-factors).
Two output-writer failures and export precision/unit repairs are disclosed in
[EXECUTION_NOTES.md](EXECUTION_NOTES.md), without deleting their artifacts.

## Experimental admission, cost and reproduction

[EXPERIMENTAL_ADMISSION.json](EXPERIMENTAL_ADMISSION.json) records the two frozen
development choices. Native decoding verifies LaIII/DyIII K1, ML=M.L, 1:1,
25°C and the stated ionic strength for all eight selected records.

| Pair | Conditional D, kcal/mol | Admission |
|---|---:|---|
| cis599/trans600, I=1.0 M | −3.66982 | K+ background noted; anion/shared medium and original provenance unresolved |
| EDTA576/trans-cyclopent597, I=0.1 M | +1.41882 | medium unresolved; adjusted cyclopent values versus unannotated EDTA |

The trans **alternate1.0-M** row was used; its primary0.5-M row was not substituted.
These D values are conditional development references, not simulation results or
prospective validation. NIST-derived rows remain in the ignored workspace under
the existing redistribution policy. TMDTA/TTHA and reserved outcomes stay closed.

Warm29-atom native energy/force evaluations took **0.1262 s for L** and
**0.0644 s for M**, ten repeats each on twelve allocated CPUs. Compilation
is separate in the timing artifact. This does not estimate bulk-solvent or H200
throughput, mixing, overlap, a converged D or its uncertainty.

Actual cost: **23 Slurm jobs,150,814 allocated CPU-seconds (41.89 CPU-hours),
146,636.412 scheduler-measured CPU-seconds (40.73 CPU-hours), zero GPU time**.
This includes failed writers, cancelled references, metadata checks and the
saved-data reproduction. Source preparation/downloads/environment setup were
unmetered. [COMPUTE_COST.json](COMPUTE_COST.json) and
[full accounting](../../workspaces/so3lr_tauri_screen_20261106/SCHEDULER_ACCOUNTING.tsv)
retain the individual records. All owned jobs are terminal; unrelated jobs remain.

On biotite, reproduce the saved analysis into a fresh directory:

```bash
bash diagnostics/so3lr_tauri_screen_20261106/reproduce.sh --saved
```

`--inference` additionally reruns the fixed L/M panel and its controls in a fresh
attempt directory using pinned existing assets. Neither mode launches ORCA.
The saved-data replay reproduced numerical tables and qualification decisions;
see [REPRODUCTION_CHECKS.json](REPRODUCTION_CHECKS.json).

No MD, alchemical ladder, peptide/protein simulation, potential training, fitted
correction or production change followed this screen. **Unresolved reference;
preserve M1 and stop before solvent expansion.**
