# Second-shell test of the La/Dy aqua-leg deficit — PLAN (frozen 2026-11-10, before any job)

Origin: `diagnostics/lady_acetate_forensics_20261109/` found the static cluster-continuum aqua leg
ΔΔG_hyd(Dy−La) = −59.21 against −66.9 experimentally (error +7.69), with the acetate exchange error
(−1.88) equal to the 24% of that error the complex leg leaves uncancelled. This test asks whether an
explicit second hydration shell carries the aqua deficit (stage A) and, only then, whether it moves acetate
(stage B). Jacob authorized the staged test on 2026-11-09. Root (the coordinating agent) owns submission.

## Stage A — aqua gate

**Systems (2 jobs).** `Dy_ssaq_w21_s1` and `Dy_ssaq_w21_s2`. Each start geometry is the final optimized
geometry (`job.xyz`) of the finished `workspaces/lady_relaxed_exchange_20261024/run_ss/La_ssaq_w21_s{1,2}`
(La + 21 waters: first shell 9 from the conformer-best La_aq9 core `conf/La_aq9/dft_c3`, plus 12 second-shell
waters), with the metal label La replaced by Dy and every coordinate unchanged. The La partners are those two
finished runs, used as they are (no rerun).

**Method.** The Dy `job.inp` is byte-identical to the La `job.inp` except (1) the metal's `%basis` block
(NewECP/NewGTO/NewAuxJGTO, Dy ECP55 effective singlet, the block the project already used for the
cancelled `Dy_ssaq_w21_s*` runs and from `diagnostics/lady_relaxed_exchange_20261024/basis_Dy.inc`), and
(2) the literal `nprocs 24` written as the placeholder `@NPROCS@`, which the job script replaces with
`$SLURM_CPUS_ON_NODE` (24) at run time, as the La runs' generator did. Method line:
`PBE0 D4 def2-TZVP def2/J RIJCOSX DefGrid3 TightSCF Opt Freq CPCM(Water)`, `Guess PModel`, SCF and geometry
MaxIter 500, charge 3, multiplicity 1 — all unchanged.

**Quantity per seed s.**
ΔΔG_hyd,ss(s) = [G(Dy_ssaq_w21_s) − E(Dy³⁺)] − [G(La_ssaq_w21_s) − E(La³⁺)], in kcal/mol (627.509474 per Eh).
G = ORCA "Final Gibbs free energy" (298.15 K, ORCA default thermochemistry, the same quantity `analyze.py`
parses). E(M³⁺) = bare-ion large-core gas single point from `workspaces/lady_ingredient_tests_20261026/gas/{La,Dy}_ion`,
exactly as `lady_ingredient_tests_20261026/collect.py`. Both clusters hold 21 waters, so no water
standard-state term enters.

**First-shell reference, chosen by the final Dy coordination number (O < 2.95 Å).**
- Dy CN 9 → reference **−57.54** (quoted earlier as −57.55, a rounding of the same value) = first-shell La_aq9 (conf/dft_c3) and Dy_aq9 (primary), both conformer-best,
  same formula, no water term (9 waters each).
- Dy CN 8 (one first-shell water moved to the second shell; composition still Dy + 21 waters) → reference
  **−59.21** = La_aq9 and Dy_aq8 (conf/dft_c2) plus one water at G_water + RT ln 24.46 + RT ln 55.34, the
  value `collect.py` reports.
- Shift(s) = ΔΔG_hyd,ss(s) − reference(CN of that seed). Negative = toward −66.9.
`second_shell_analyze.py` recomputes both references from saved outputs and must reproduce them to 0.01;
otherwise it stops.

**Seed admissibility.** A seed counts only if: ORCA terminated normally with Opt converged and Freq done;
no imaginary mode below −30 cm⁻¹ (the `conformer_best` rule); Dy CN is 8 or 9; no O lies in the
ambiguous 2.95–3.10 Å band; the La partner passes the same checks. Any other outcome makes that seed
inadmissible. A Dy CN change is therefore handled by the reference rule above, never by rerunning.

**Gate (needs both seeds admissible).** m = mean of the two shifts, d = |shift(s1) − shift(s2)|.
- PASS if m ≤ −4.0 and (d ≤ 1.0 or both shifts ≤ −4.0).
- FAIL if m > −4.0 and (d ≤ 1.0 or both shifts > −4.0).
- Otherwise (seeds disagree by more than 1.0 and straddle −4.0), or either seed inadmissible: INCONCLUSIVE.

**Stop rule.** FAIL or INCONCLUSIVE ends the test: stage B is not released, no third seed, no rerun, no
change of method. FAIL is reported as: an explicit second shell at this level does not carry the aqua deficit.

**Cost cap and runtime action.** Hard cap 1,000 allocated CPU-h for stage A, summed over all its job IDs
(main and restart; `sacct` CPUTimeRAW). Expected 18–20 h wall per job on 24 CPUs (the La partners: 19.9 h,
221 cycles; 18.2 h, 187 cycles) ≈ 430–480 CPU-h each. If any stage-A job passes 30 h wall, or the summed
allocation reaches 1,000 CPU-h first, root touches `STOP` in that run directory, cancels its restart job by
exact ID, then the main job by exact ID; the seed is recorded as inadmissible (unconverged) and the gate is
INCONCLUSIVE. The afterany restart job runs only when the main job ended without a Gibbs energy and no
`STOP` file exists; its CPU-h count toward the cap.

## Stage B — acetate (released only if the gate passes)

Matched CN9 monodentate quadruple, first-shell value −1.97 (La_acmono_w8, Dy_acmono_w8, La_aq9, Dy_aq9; n = 0;
DLPNO-CCSD(T) − PBE0-D4 correction +0.09). Two matched seed pairs of `La_ssacmono_w20` and `Dy_ssacmono_w20`
(acetate monodentate, 8 first-shell + 12 second-shell waters; Dy built by swap into the finished La geometry,
same method rules as stage A); 4 jobs, ≤ 2,000 CPU-h. Exchange per (complex seed i, aqua seed j) =
[G(Dy_ssacmono_w20_i) − G(La_ssacmono_w20_i)] − [G(Dy_ssaq_w21_j) − G(La_ssaq_w21_j)] (n = 0); report the mean
over the four (i, j) combinations; acetate shift = mean − (−1.97).
- Under (a) the acetate shift is about 0.24 × the stage-A aqua shift, so ≥ +1.0. Under (b) or a
  complex-intrinsic error it is about 0.
- **Pass for (a):** shift ≥ +0.7 and shift / |aqua shift| ≥ 0.15.
- **Against (a):** shift < +0.3.
- Between: inconclusive.

## Files

`build_stage_a.py` (writes the two inputs) · `orca_job.sh`, `restart_job.sh` (copied from
`lady_relaxed_exchange_20261024`, changed only to `-c 24`, `$SLURM_CPUS_ON_NODE` ranks, and a `STOP` check)
· `submit_stage_a.sh` (prints and runs the sbatch chain; root runs it) · `second_shell_analyze.py` (gate as
JSON). Runs in `workspaces/lady_second_shell_20261109/runs/`.
