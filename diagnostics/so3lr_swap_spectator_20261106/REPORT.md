# Does the distant-water effect cancel in the metal swap?

6 November 2026. **L passes the predeclared energy test. M misses it.**

Adding the same distant water to both arrangements changed L's predicted swap
energy by at most **0.068 kcal/mol**, below the **0.5-kcal/mol limit**. This is
new evidence that the earlier large single-complex effect mostly cancels for L
in this prototype. It supports keeping L as a candidate. It does not establish
chemical accuracy, reliable forces throughout an ensemble, or a binding free
energy.

**What we tested.** The two arrangements were LaAcW7+DyW8 and DyAcW7+LaW8.
Both contain the same54 atoms, total charge+5 and multiplicity6. Only the metal
identities swap. The internal coordinates came unchanged from the earlier real
fixtures; the complexes were translated40 A apart. The same archived water was
added to both arrangements at100,200 and400 A, in two opposite orientations.

This is a charged isolated acetate/aquo prototype. The intended ligand experiment
would use a neutral periodic box. Passing or failing this prototype alone cannot
validate or reject that later system.

The measured quantity was the change in the swap energy caused by the extra
water: [E(B+water)-E(A+water)]-[E(B)-E(A)]. No isolated-water reference energy is
needed for that difference. Global composition, charge and spin match between
the two endpoints of every comparison.

| Model | Maximum effect at200/400 A, LR1000 A | No-cutoff repeat | Declared energy criterion |
|---|---:|---:|---|
| **L, primary** | **0.0670 kcal/mol** | **0.0677 kcal/mol** | **Pass** |
| M, fixed sensitivity | 0.5821 kcal/mol | 0.5827 kcal/mol | Misses0.5 limit |

The energy effect persists with distance, rather than vanishing completely. For
L, changing200 to400 A changes it by about0.002 kcal/mol. The estimated ordinary
remote electric interaction changes the swap energy by at most0.000153 kcal/mol
at the scored distances. Reversing the water and removing the cutoff preserve
the conclusion. M was not substituted for the preselected L model.

**An important warning: energy cancellation did not imply force cancellation.**
The largest change in a component of the endpoint force difference was about
**4.18 kcal/mol/A for L**; the acetate-metal force-difference vector changed by
about0.50 kcal/mol/A. These were recorded diagnostics, not an added acceptance
criterion. The declared energy test still passes. However, the small energy
effect at this fixed geometry need not remain small as the complexes move.

The bare joint-system swap energy was+11.857 kcal/mol for L and+10.480 for M.
For context, the earlier separate-cluster C_Ac values were+3.963/+10.769. Joining
the charged complexes introduces real interactions and changes the global model
conditioning. Without a qualified reference, neither set of absolute predictions
is accuracy evidence. Partial-charge totals on the complexes and spectator,
arithmetic component changes, and all force diagnostics are retained in
[swap_changes.tsv](swap_changes.tsv) and [RESULTS.json](RESULTS.json).

**Checks and reproducibility.** The frozen [plan](PLAN.md) was declared before
new values. We reused unchanged source70707f5, GLP2294f4f, L/M-2606 weights and the
existing Python3.12/JAX0.10.2 environment. Float64, highest matmul precision,
head5, explicit charge/spin, native loader and standard inference were retained.
The neighbor list was5 A; each model retains its native radial cutoff, L5/M4.5 A.
The controls used LR1000 A/damping2 A and no cutoff/dampingNone. No fragment
charge constraints, embedding overrides, fitting or optimization were applied.

All geometry/state/head, complete-pair, charge-conservation, finite-value and
energy-component checks passed. Real LaW8 replays agree with the previous
predictions within3e-14 eV and6e-15 eV/A. Raw float64 arrays were saved before
formatting. [Artifact hashes](ARTIFACTS.json) retain the executed source copies,
inputs and all62 model outputs. A saved-array replay reproduces both numerical
tables byte for byte and reproduces the decision; see
[REPRODUCTION_CHECKS.json](REPRODUCTION_CHECKS.json).

**Cost:** inference job52169 took53 seconds on12 allocated CPUs; saved replay
52184 took4 seconds on one CPU. Total **640 allocated CPU-seconds, or10.7 CPU
minutes**, within the approved60-minute budget. Scheduler-measured CPU use was
175.149 seconds. No GPU, DFT or MD was used. Both jobs completed. Source
preparation was unmetered. [Accounting](COMPUTE_COST.json).

Reproduce the analysis into a fresh workspace directory on biotite:

```bash
bash diagnostics/so3lr_swap_spectator_20261106/reproduce.sh
```

This command reads saved arrays; it does not make new model or quantum calls.

**Decision:** retain L for the next limited investigation. Prepare the intended
neutral ligand system and check energy/force sensitivity there before starting
trajectories. This result does not yet justify a selectivity free-energy
calculation. Test2/MD was not started. M1, historical results, reserved outcomes,
other-session jobs and unrelated edits remain unchanged.
