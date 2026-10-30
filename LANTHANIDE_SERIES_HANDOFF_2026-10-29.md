# Lanthanide-series calibration: session handoff (2026-10-29 12:00)

## Goal

A scalable, physics-based prediction of lanthanide selectivity across the whole series (La–Lu), usable as a discriminator for protein sites. La/Dy is the worked example, chosen because of Hans-LanM; it is not the goal. The route:

1. **Calibrate (failed 2026-10-29; see the 19:30 update).** Test whether the careful cluster-continuum DFT protocol can be calibrated against NIST 46 data on La/Dy. This is the panel in `diagnostics/lady_panel_20261028/`.
2. **Extend to Nd and Lu (not launched; its gate failed).** If a calibration form passes, run Nd and Lu on the same panel, giving one calibration per metal. Staged in `diagnostics/lady_series_20261028/`.
3. **Reconstruct the series.** Map computed La, Nd, Dy and Lu onto all 13 metals with the NIST two-shape basis.
   - Leave-one-ligand-out error is 0.12 kcal/mol, about experimental noise.
   - Protein and peptide series fit the same two shapes. Mex-LanM's 14-metal series fits to 0.12 kcal/mol, and predicting it from La/Nd/Dy/Lu alone is off by 0.17.
   - Details are in the series README.

## Update 2026-10-29 13:00 (next session)

- **Watcher re-armed** as task `bjczmry32` (12:46, after the Dy TTHA restart). The session was `/clear`ed, not restarted, so the old watcher survived. It was stopped by its exact PID.
- **TMDTA:** both La geometries converged.
  - La_tmdta_w3 converged at **CN 8**: the third water left. It is excluded by rule.
  - La_tmdta_w2 (CN 8, as intended) carries La-TMDTA. Both are in Freq.
- **La_cdta_w3** is at TightOpt threshold level (cycle 160) and still has 3 waters.
- **Dy_ttha_w0 restarted** from La_ttha_w0's DFT geometry with Dy swapped in (jobs **30037/30038**; JOBS.txt updated).
  - The Dy MACE start's strained backbone made TightOpt oscillate: +0.26 Eh steps, and max gradient still 0.14 at cycle 7.
  - It is the same conformer as La. Details are in the panel README. The old run is in `abandoned/`.
- **Nd/Lu code changes are done and dry-run tested.** Details are in the series README. TTHA is added, and the H3 form plus the panel-selected whole-series form are implemented.
  - The four-ligand calibration is reported as a secondary check. For Nd, acetate and EDTA are only 1.2 kcal/mol apart experimentally.
  - Four dead-end states were dropped, leaving 33 states.
- **La_tmdta_w3 cancelled** (28689/28690, 13:20) during its Freq step. It had already converged to CN 8, so the rule excludes it whatever the Freq gives.
- **The restarted Dy_ttha_w0 is smooth:** energy falls steadily, with max gradient 0.017 → 0.003 by cycle 3.
- **Recommendation (13:20): finish the panel, don't launch Nd/Lu as designed.**
  - The panel's remaining cost is small, about 1.5k CPU-h, and it decides whether the method can be calibrated at all.
  - The series costs about 12k CPU-h and buys only 0.28 → 0.12 kcal/mol in reconstruction error, below the 0.8 calibration error being tested. It would pay off only if the calibrated error came in near 0.2–0.3.
- **Launch still waits for the panel.** The leave-one-out RMS of the two forms on the four fitted ligands is 0.77 and 0.91, right at the 0.8 bar. A pass is far from certain, so launching early would be an expensive gamble.

## Update 2026-10-29 19:30: panel decided, it fails

- **Verdict:** CDTA and TMDTA finished at calc −11.7 and −10.5, against exp −5.7 and −4.9.
  - Every frozen form misses them by 1.8–4.1 kcal/mol.
  - Even with malonate and TTHA perfect, the best held-out RMS is 1.30 (H3), so the pre-registered test fails.
  - Ranking (CDTA < TMDTA < EDTA) is right.
- **Cancelled the remaining runs** at 19:20, restart job first: La_mal_w7 (28737/28736) and La_ttha_w0 (29780/29779), then Dy_ttha_w0 (30038/30037).
  - Watcher `bjczmry32` exits on its next poll and reruns `panel_analyze.py`; the verdict is unchanged.
- **Root cause** (panel README, "Result"): CDTA and EDTA are near-identical ligands, but their calculated values differ by 3.7 kcal/mol more than experiment.
  - The difference is an 11 kcal/mol swing in the CPCM term, tied to the La-CDTA waters H-bonding back onto the carboxylates (3 such H-bonds, against 1 in La-EDTA).
  - Static, single-structure cluster-continuum is therefore sensitive to the water arrangement at the several-kcal/mol level.
- **Post hoc, not a pass:** a stretch refit on all six ligands has a leave-one-out error of 0.79, but its slope moved from 1.28 to 1.62.
- **The Nd/Lu series is not launched** (its gate failed).
- **Open decision for Jacob:**
  - stop the static protocol (recommended); the next method has to average over the first-shell waters;
  - or first run one water-orientation check on La-CDTA (~1k CPU-h, ~1 day).

## Update 2026-10-29 23:45: next direction proposed, waiting on outside review

- **Static route stopped.** I advised against the La-CDTA water check. Even if it rescued CDTA, the fix (a water search for every complex) does not scale.
- **Proposal: a NIST-trained predictor** of the one-number-per-ligand selectivity.
  - It uses chelate-ring and backbone features, plus one computed descriptor: force-field size-fit strain at La versus Dy metal–donor distances.
  - The full write-up, with pre-registered pass/fail, is `LANTHANIDE_NIST_PREDICTOR_REVIEW_2026-10-29.md`. It is self-contained and limited to small-molecule chemistry, so Jacob can hand it to other models.
  - Jacob is collecting outside opinions before deciding.
- **Why it could work** (`diagnostics/ln_selectivity_calibration_20261027/DATA_AUDIT_20261029.txt`):
  - 7 ligands carry 63% of the donor-count model's squared error; the other 94 are fit to 0.54 kcal/mol.
  - 15 diamine tetraacetates with the same donors span 4.3 kcal/mol, and 11 of them are missing from the current model.
  - The data can grow from 101 to 208 ligands: 58 La/Dy pairs were lost to the PubChem lookup, 12 to the I ≤ 0.5 filter, and 37 more ligands can enter through the series shape.
  - The noise floor is 0.29 kcal/mol, or 0.57 for strong chelators.
- **The NIST pipeline is now tracked** in `diagnostics/ln_selectivity_calibration_20261027/` (scripts, README, audit). The NIST-derived tables stay in the untracked workspace because the repo is public.
- **Committed and pushed:** the La/Dy record since 2026-10-23 (the `lady_*` and `ln_surrogate_20261027` diagnostics, this handoff, `ALMO_SESSION_HANDOFF_2026-10-24.txt`). Other projects' uncommitted changes in the tree were left alone.
- **Nothing is running and no watcher is needed.** Nothing new launches until Jacob decides.

## Restarting the session

Nothing is running for this work (2026-10-29 23:45), so a restart loses nothing and no watcher needs re-arming.

## First steps in the new session

1. Read `LANTHANIDE_NIST_PREDICTOR_REVIEW_2026-10-29.md`, then ask Jacob for his decision and for any outside opinions.
2. If he says go, follow section 5 of that document in this order:
   - structure curation;
   - 2D features;
   - strain descriptor;
   - fits, scored by the pre-registered criteria.

   All of it is single-threaded and runs on the login node.
3. Everything below this point describes the finished panel run. Rerunning `diagnostics/lady_panel_20261028/panel_analyze.py` reproduces the verdict.

## Jobs in flight (2026-10-29 11:57)

Each state has an Opt+Freq job (`pn_`) and an afterany restart (`pnr_`). The restart exits "already done" if TightOpt finished; otherwise it restarts from the last geometry with normal Opt + Freq.

| State | Jobs | Status |
|---|---|---|
| La_cdta_w3 | 28669 / 28670 | Running 22.6 h, more than 153 opt cycles. Still moving (max step 0.06): the floppy waters. |
| La_tmdta_w3 | 28689 / 28690 (standard, 12 ranks) | Running 22.6 h. Gradients converged, only max step left, so close to done. Its third water sat at 3.08 Å, just inside the 3.1 Å CN cutoff. If that water leaves, the state is excluded and w2 carries TMDTA. |
| La_tmdta_w2 | 29036 / 29037 (standard, 12 ranks) | Running 18.5 h. Gradients converged. |
| Dy_nta_w4 | 28744 / 28745 | Running 3.8 h. |
| La_mal_w7 | 28736 / 28737 | Running 1.1 h. Unchained from La_cdta_w3 today. |
| La_ttha_w0 | 29779 / 29780 | Running 1 h. New today. |
| Dy_ttha_w0 | ~~29781 / 29782~~ → 30037 / 30038 | Restarted 12:45 from La's geometry (see the update above). |

- **Finished, terminated normally:** Dy_cdta_w2, Dy_dtpa_w1, La_dtpa_w1, Dy_mal_w6, Dy_nta_w5, La_nta_w5, Dy_tmdta_w2.
- **Jacob's jobs** (bmeth_ipdb, bmeth_agg, fs245L, rrna_cap, rcnv_x21): don't touch. rrna_cap and rcnv_x21 were user-held at handoff.

## Results so far

La→Dy exchange ΔΔG in kcal/mol; negative means Dy is favoured. Acetate and EDTA states are conformer_best.

| Ligand | Calculated | NIST | Bias | Waters displaced |
|---|---|---|---|---|
| acetate | −1.95 | −0.07 | −1.88 | 1 |
| NTA | −4.26 | −1.58 | −2.68 | 3.5 |
| EDTA | −6.35 | −4.01 | −2.34 | 6 |
| DTPA | −8.33 | −4.56 | −3.78 | 7.5 |

**H0 (constant offset from acetate) has failed, with certainty.**
- EDTA, NTA and DTPA alone give RMS ≥ 0.86 over the six test ligands. The threshold was 0.8.
- The NTA residual can only get more negative. Under the state rule, an extra state can only lower its family's G.

**The calculated order matches experiment for all four ligands.**

**Two post hoc two-parameter forms fit the four ligands.** Both were frozen 2026-10-29, before any held-out ligand finished:
- **H2, per-water bias:** bias = −1.642 − 0.2282·n_disp. In-sample RMS 0.41.
- **H3, stretched selectivity:** calc = 1.2792·exp − 1.956. In-sample RMS 0.48.
- The primary-only fits have the same slopes, with intercepts +0.141.

**Leave-one-out RMS on these four:** constant 0.93, +q_c² 1.35, +n_disp 0.77, +exp 0.91.

## Pre-registered tests still pending

These are coded in `panel_analyze.py` (the `FROZEN` block and TTHA call) and described in the panel README.

1. **Frozen forms on the held-out ligands (malonate, CDTA, TMDTA, TTHA).** A form passes if its RMS is below 0.8 kcal/mol, measured in experimental units (H3 residuals divided by 1.279).
2. **The TTHA call.**
   - Why TTHA: it displaces every water (n_disp 8.5) but is nearly unselective (NIST −0.68). That is the protein-site regime, where H2 and H3 extrapolate differently.
   - Call: calc ≤ −3.90 → per-water form; calc ≥ −3.20 → stretch or constant form; in between → inconclusive.
   - These cutoffs are for conformer_best. Shift them by +0.141 for primary_only.
3. **Order test:** CDTA < TMDTA < EDTA.

Frozen predictions of the calculated ΔΔG (conformer_best):

| Ligand | H0 | H2 | H3 |
|---|---|---|---|
| TTHA | −2.56 | −4.26 | −2.83 |
| CDTA | −7.55 | −8.69 | −9.21 |
| TMDTA | −6.75 | −7.88 | −8.18 |
| malonate | −2.53 | −2.75 | −2.79 |

**State rule:** use the lowest-G, water-balanced hydration state per metal. A state whose CN (donors within 3.1 Å) changes during DFT is excluded. TTHA is exempt from the CN check: it has no water to lose, so 9 vs 10 bound donors are isomers of the same state.

## Decision tree when the watcher reports

- **At least one frozen form passes on the held-out ligands:** launch Nd/Lu (next section).
  - The states are the same whichever form wins. The form only matters at analysis, where the TTHA call picks it.
- **Borderline:** run the tier-2 alternate hydration states. Their inputs already exist in `workspaces/lady_panel_20261028/runs/`, with no job.out yet: La_cdta_w2, Dy_cdta_w3, Dy_tmdta_w3, La_nta_w4, La_mal_w6, Dy_mal_w7.
  - Each can only lower its family's G, so the direction of each effect is known in advance.
- **Both forms fail:** static cluster-continuum can't be calibrated to 0.8. The series run is moot; report that directly, without spin.
  - The physical lead is a per-water bias. CPCM misses the second hydration shell, which is stronger (ΔH) and more ordered (TΔS) around the smaller Dy. That would explain the ΔH and TΔS errors in `diagnostics/lady_relaxed_exchange_20261024/EXPERIMENT_THERMO.md`: both grow with waters displaced and partly cancel in ΔG.
  - That idea is untested. Discuss it with Jacob before spending compute on it.

## Launching Nd/Lu (staged in `diagnostics/lady_series_20261028/`; README there)

**Code changes required before launch (DONE 2026-10-29 13:00; see the update above):**
- **`build_series.py`:**
  - add `ttha_w0` to `PANEL_STATES` (Nd from La_ttha_w0, Lu from Dy_ttha_w0);
  - set `N_LIG_DONORS["ttha"] = 9`;
  - skip the CN-match test in `source()` for TTHA.
  - Fallback starts if needed: `workspaces/ln_surrogate_20261027/nist/678/{Nd,Lu}_cn8_mace.xyz`.
- **`series_analyze.py`:**
  - add `LIG_ID["ttha"] = "678"`, `Q_C["ttha"] = -3` and `DONORS["ttha"] = 9`, plus the CN exemption;
  - add the H3 form: per metal, fit calc = α·exp + β on acetate + EDTA, then score the rest;
  - make the whole-series block use the panel-selected form. It currently subtracts only the acetate offset.
- **`submit_series.sh`:** add `ttha_w0` near the front of `ORDER`, since it is a large state.

**Steps:**
1. Rebuild with the final panel geometries: `/groups/banfield/users/jwestrob/envs/lnsel/bin/python build_series.py`.
2. Check `squeue --me`.
3. Run `./submit_series.sh N_SE N_STD`.
4. Arm `workspaces/lady_series_20261028/watch_series.sh` in the background.

**Size:** about 37 states at 10–20 h each on 24 ranks. The pinned Nd/Lu basis files are verified; regenerating La reproduces the production AuxJ exactly.

## Translating this into a discriminator (discussed with Jacob 2026-10-28)

- **Precision:** a per-metal bias cancels when ranking sites. An error of 0.8 kcal/mol is about 4× in separation factor, enough for three bins.
- **Known gaps:**
  - Amide and backbone donors: there are none among the 75 NIST La+Dy ligands.
  - Rigidity: compare frozen vs free attachment atoms, plus an MLIP/MM elastic term.
  - Cost: about 2 days per state for 60–80-atom clusters, at 4 states per site. That is fine for about 10 confirmed families.
- **Bridge test:** LBT3, where the ITC value is −2.36 kcal/mol (Hatanaka 2020). Pass if |prediction − (−2.36)| < 1.0.
- **Structure checks:** crystal pairs Hans-LanM 8DQ2/8FNR and Mex-LanM 6MI5/8FNS.
- **Out of scope:** Hans-LanM dimer allostery, by construction.
- **Data:** protein affinities are in `workspaces/ln_selectivity_calibration_20261027/protein_ln_affinities.tsv`.

## Jacob's standing instructions (beyond CLAUDE.md and memory)

- **Slurm:** never put `--mem` or `--time` in sbatch. Cancel exact job IDs only. Arm a watcher after every submission.
- **Login node:** at most 16 threads in total; run single-threaded analysis inline.
- **Python:** `/groups/banfield/users/jwestrob/envs/lnsel/bin/python`.
- **Scope:** stay on lanthanide physics. Don't run the GPU fold.
- **Autonomy:** no permission needed ("yolo"). Still no expensive checks whose result would be inconclusive.
- **Queue use:** throttle only while Jacob's queue is busy.
- **Reporting:** keep reports short and high-level. Jacob objected to getting "bogged down in the weeds"; the target is the whole series.
- **Never** offer "flag for a wet-lab experiment" as the solution.
- **Email** jacobwestroberts@gmail.com (`/usr/sbin/sendmail -t`) only for genuine blockers.

## Gotchas

- Running `panel_analyze.py` by hand rewrites `PANEL_RESULT.json`. That is harmless: the watcher rewrites it at the end.
- `runs/*/job.inp` can still show `@NPROCS@` after a run, because the inputs were regenerated after the jobs started. The ORCA echo in job.out shows the real count (24).
- ORCA sbatch needs `-N 1`. `orca_job.sh` and `restart_job.sh` in `diagnostics/lady_relaxed_exchange_20261024/` already have it. Use `< /dev/null` inside while-read pools.
- `squeue -j <list>` silently drops purged IDs. Use sacct for final states.

## Key paths (repo = `alchemical_bvs/`)

| What | Path |
|---|---|
| Panel README, analysis, TTHA build | `diagnostics/lady_panel_20261028/{README.md,panel_analyze.py,build_ttha.py}` |
| Panel runs, JOBS.txt, watcher, result | `workspaces/lady_panel_20261028/` |
| Series staging (build, analyze, submit, basis_Nd/Lu.inc) | `diagnostics/lady_series_20261028/` |
| Series runs (35 dry-run inputs) and watcher | `workspaces/lady_series_20261028/` |
| Original La/Dy acetate/EDTA runs and shared helpers | `diagnostics/lady_relaxed_exchange_20261024/` (`analyze.py`, `EXPERIMENT_THERMO.md`) |
| NIST ligand table and MACE starts | `workspaces/ln_surrogate_20261027/{surrogate_nist.tsv,nist/<id>/}` |
| Session transcript | `/home/jwestrob/jwestrob/.claude/projects/-groups-banfield-projects-environmental-sr-srvp2020-Jacob-lanthanide-binding-on-density-scanner-alchemical-bvs/39867d9b-3f9a-4f8f-86e2-21cbaede90a8.jsonl` |
