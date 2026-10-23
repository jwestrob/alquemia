# Canonical method check of the consumed dry mono-acetate rotation

## Fixed question and inventory

Does the opposite model/reference rotational energy slope persist with canonical RHF/CCSD(T), or is it sensitive to replacing PBE0? This is a fixed-geometry method comparison, not a fit, a CP-corrected force, or a physical-truth claim.

Exactly **eight new complex single points**: La and Dy × the existing rotation q=−.010, −.005, +.005, +.010 Å from the admitted dry-acetate derivative panel. The original coordinate is q=(1 Å)θ; preserve the axis, ligand centroid, all literal coordinates and atom order. No center, ghost-fragment, gradient, optimization, new displacement or electronic field calculation is included. Existing canonical bidentate-acetate energies have different geometries and provide recipe/cost evidence only. They do not replace these mono-acetate endpoint energies.

## Unchanged canonical recipe

Reuse donor_cp_v1's qualified deck and parser: `RHF CCSD(T) def2-TZVP NoFrozenCore VeryTightSCF NoAutostart`, SCF Guess PModel / MaxIter300; MDCI MaxIter100 / STol1e−8 / PrintLevel3. Use the exact original metal GTO/ECP files and light-atom basis. All eight centers are real; no auxiliary correlation fit, local-correlation approximation, RIJCOSX, D4, solvent or EnGrad is introduced. Preserve the original ORCA6.1.1 executable.

Each system has charge+2, effective singlet RHF, 40 explicit and correlated electrons, 189 orbital functions with per-center counts [47,31,31,31,31,6,6,6]. La ECP46 and Dy ECP55 leave effective metal nuclear charge11. Physical parent La singlet / Dy sextet metadata remain distinct from the effective closed-shell Hamiltonian. NoFrozenCore correlates all remaining electrons, including light-atom cores; it does not restore the removed ECP electrons.

Each endpoint uses **four MPI ranks on one node**, one thread per rank, standard_expand with default memory/time. MaxCore=floor(.75 × actual memory grant /4); 16,000 MiB gives3000MiB/rank. No memory/time override or hidden resource-driven scientific input change. Runtime verifies the actual single-node grant and retains the exact submitted deck, argv, START, return code, stdout/stderr and elapsed allocation cost.

## Staging and native admission

The fixed first scout is `cc__La__rotation__minus_010__complex`. Root reviews its actual native admission and measured cost before releasing exactly the other seven. The scout is reused in the full eight denominator. No duplicated endpoint, automatic retry or alternative solver is allowed. An existing output directory is never reused. Source preparation performs zero native calls.

Retain the donor parser's actual coordinate, per-center orbital basis, effective charges, ECP assignments, RHF, charge/electron/multiplicity, NoFrozenCore, normal termination, native SCF/CC convergence, CC/final-energy closure and declared final max-residual≤1e−8 gates. Retain SCF stopping mode and every printed criterion, T1, amplitudes, triples contribution and CC trace. Failed printed criteria remain visible under the original native stopping rule; no new state/accuracy guarantee or post-outcome diagnostic cutoff is introduced.

## Saved comparison and prospective numerical resolution

For each metal and q∈{.010,.005} Å, compute

`dE/dq = [E(+q)−E(−q)]/(2q)`, in kcal/mol/Å,

and the Dy−La difference. Use the native final CCSD(T) energy token after the unchanged CC/final closure gate, retain its decimal precision and separately retain printed CCSD/triples values. Report both steps, their signed difference and propagated half-last-digit bounds. Apply the existing **0.1 kcal/mol/Å** resolution only to the two-step consistency of the La, Dy and paired CCSD(T) derivatives: three required checks, all unavailable if the full eight-cell denominator is incomplete. A failed two-step check remains failed; no Richardson extrapolation, changed steps, tolerance or added points.

Compare both CC slopes with the already admitted PBE0 electronic-only, PBE0-D4, separate D4 derivative and fixed selected17621 model slopes/components. Comparisons of signs and magnitudes are **descriptive**; there is no new cross-method physical-accuracy tolerance. D4 is not added to CC. Existing PBE0 analytic center projections remain contextual evidence; no CC analytic gradient or CP correction is claimed.

If the canonical slope supports the PBE0 sign and broad size, it independently supports a model defect for this consumed rotational direction within the shared frozen-f/basis approximation. A substantially different canonical slope establishes method sensitivity, not automatic replacement of targets. Unresolved numerical checks prevent a qualified directional conclusion. Missing diffuse functions for acetate, basis/semicore convergence, physical f-electron/spin-orbit state, single-reference adequacy and physical force accuracy remain open. Canonical CC and AMOEBA components must not be given a one-to-one mechanistic interpretation.

## Measured cost and retained evidence

The admitted same-size donor controls at189AOs/40correlated electrons took La1750.337s and Dy1856.244s, four ranks each (jobs17552/17555). Four analogous endpoints per metal estimate57,705 allocated CPU-seconds, **16.03 CPU-hours**, for this panel. Their geometries are bidentate and timings are estimates, not guarantees; the actual mono-acetate scout controls continuation. The donor201AO scout approached its default memory grant, so actual native/memory behavior remains a limitation rather than a reason to change settings silently.

`PANEL.json` binds the original derivative panel and admitted quantum/model comparison, exact canonical method sources/assets, two cost controls and their strict review. All original failed controls and results remain untouched. Root owns release and submissions; preparation and collection invoke no new molecular program.
