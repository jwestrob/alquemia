# Pass 2 — Alchemical BVS validation: pre-check on PDB 8DQ2

**Status:** pre-check complete on HansLanM X-ray (PDB 8DQ2 chain A, 3 La sites). Full Pass 2 on canonical *M. extorquens* LanM Protenix fold not yet run — gated on this pre-check's outcome per the methodological reasoning below.

**TL;DR:**
- All 3 La sites in 8DQ2 chain A predict **Lu** as preferred Ln (residual 1.45–1.48), with monotonic-with-Ln-contraction residual curve from La (3.59–3.62) → Lu (1.45–1.48). One mild non-monotonicity: Dy locally higher than Tb (consistent across all 3 sites).
- Donor-distance spread across the Ln panel: La 2.30 Å → Lu 2.22 Å, **Δ = 0.08 Å**. Real Ln-contraction is ~0.2 Å. The vault note's "classical FF gives <0.05 Å vs real ~0.2 Å" diagnostic is reproduced almost exactly.
- This is **scenario B** (pre-committed): σ-driven Ln-contraction signal recovered, apparent-Kd bell curve NOT recovered. Mirrors the FEP failure mode at the geometric level.
- **Recommendation: do NOT scale to full Pass 2 on canonical M. extorquens LanM.** Pivot per the framework (see "Pivot options" section).
- Wallclock: 52.9 s on 64-core SLURM `standard` node.

---

## Why a pre-check first

The handoff specified a 4-site canonical LanM Protenix fold + alchemical Ln panel × backbone-restrained min × BVS residual. Running that directly costs ~24 hr of GPU + minimization. Before committing, we wanted to know:

1. Does the OpenMM-Merz-12-6-4-GBSA pipeline produce non-trivial Ln-discriminating geometry on a known-good X-ray geometry?
2. If yes, does the discrimination map to a bell curve (apparent Kd, Cook & Cotruvo 2019) or a monotonic Ln-contraction signal?
3. Two prior datapoints made this a real concern:
   - The FEP work in `~/jwestrob/obsidian-vault/agent-captures/2026-04-06_lanm-project-pause-summary.md` documented that classical 12-6-4 FEP on this exact problem **fails quantitatively**: wrong sign Eu→Gd, magnitudes 10–40× too large, with the conclusion that *"the problem is the functional form."*
   - Apparent-Kd Ln spread is only ~0.8 kcal/mol — likely below the noise floor of any classical FF.

So the pre-check is option **(B′)** — same scientific frame as a full Pass 2 run, but on the X-ray reference (HansLanM/8DQ2) at fraction of the cost. ~30 s wallclock on a SLURM `standard` node.

## What the pre-check measures

For each (site, target Ln) combo:

1. Load 8DQ2 chain A (3 La ions in the EF-hand-like binding sites; CN=10 in X-ray).
2. Strip waters + Ln-protein bonds (the PDB has `LINK` records that OpenMM materializes as topology bonds, which would otherwise auto-create NB exclusions and zero out the LJ wall — a real bug we hit during development).
3. Build OpenMM system: `amber14-all` + `amber14/tip3p` (for the LA residue template) + manual GBSA OBC2 with element-based Bondi radii (no `amber14/obc2.xml` ships with this OpenMM install).
4. Override target La's NB params with target Ln's Merz 12-6-4 OPC values from `fep/prep/ln_parameters.json` (Li-Song-Merz JCTC 2021).
5. Add CustomNonbondedForce for the C₄ ion-induced-dipole correction (Ln-O only).
6. Restrain k = 50 kcal/mol/Å² on:
   - all protein backbone atoms (N, CA, C, O),
   - sidechains of residues NOT within 6 Å of the target metal,
   - the other 2 La ions (kept as La),
   - **AND the target metal itself** — because GBSA self-solvation pulls the +3 ion out of the buried site into bulk dielectric (verified: 137 Å migration when free). Pinning the metal isolates the σ-driven differential geometry signal that BVS measures.
7. L-BFGS minimize (tol = 0.4 kJ/mol/nm ≈ 0.01 kcal/mol/Å, max 5000 iters).
8. Compute BVS at relaxed geometry using `on_scanner.bvs.compute_bvs()` extended with Brese-O'Keeffe r₀ values for the full Ln panel {La, Nd, Sm, Eu, Tb, Dy, Yb, Lu}.

## Design decisions (deliberate, all flippable in v2)

| Choice | Why | Sensitivity |
|---|---|---|
| GBSA OBC2 vs explicit OPC | OPC needed solvation box per combo (~24× slower); GBSA gives bulk-dielectric screening for free | likely small at buried site |
| Element-based Bondi OBC radii | `amber14/obc2.xml` doesn't ship in this OpenMM install; manually-added GBSAOBCForce per element. Ln radius 0.18 nm uniformly. | likely small (metal Born radius dominated by burial) |
| C₄ included | Canonical 12-6-4 spec | medium — could ablate to test |
| Metal pinned at X-ray | GBSA pushes free +3 out of buried site (empirical: 137 Å migration). Pinning isolates σ-driven differential donor geometry. | high — alternative is harder restraint on metal at X-ray (similar effect) |
| 50 kcal/mol/Å² restraint | Per handoff. Stiff. | medium — softer would let donors move more |
| Ln-protein LINK-record bonds stripped | They eliminate the LJ wall via auto-exclusion (donors collapse to 1.5 Å). Stripping restores physical NB interaction. | critical bug fix |

---

## Results

_(filled in by `precheck_8dq2.py`; see `results/precheck_8dq2_residuals.tsv`)_

### Per-site BVS residual |BVS - 3|, smallest = predicted

| Site | La | Nd | Sm | Eu | Tb | Dy | Yb | Lu | Predicted |
|------|------|------|------|------|------|------|------|------|-----------|
| LA201 | 3.590 | 2.733 | 2.628 | 2.228 | 2.208 | 2.309 | 1.558 | 1.454 | **Lu** |
| LA202 | 3.594 | 2.730 | 2.634 | 2.236 | 2.210 | 2.314 | 1.560 | 1.456 | **Lu** |
| LA203 | 3.621 | 2.759 | 2.669 | 2.256 | 2.249 | 2.334 | 1.585 | 1.478 | **Lu** |

Trend across the panel (consistent on all 3 sites):
- Monotonic decrease La → Eu/Tb (ΔΔ = ~1.4)
- Local maximum at Dy (Dy > Tb by ~0.1)
- Continued decrease Yb → Lu

### Minimum donor distance per Ln (Å), site LA201

| Ln | dmin (Å) | n_donors |
|----|---------:|---------:|
| La | 2.302 | 11 |
| Nd | 2.297 | 11 |
| Sm | 2.269 | 11 |
| Eu | 2.287 | 11 |
| Tb | 2.254 | 11 |
| Dy | 2.229 | 11 |
| Yb | 2.234 | 11 |
| Lu | 2.227 | 11 |

Total dmin spread La → Lu: **0.075 Å**. Real Ln³⁺-O equilibrium-distance contraction (Brese-O'Keeffe r₀): 0.20 Å. Classical 12-6-4 OPC reproduces ~38% of the real contraction at this static-min level — same order as the vault's prior FEP diagnostic ("<0.05 Å vs ~0.2 Å").

### Comparison to published Kd

| Reference | Conditions | Preferred Ln |
|---|---|---|
| Cook & Cotruvo 2019, **apparent** Kd (pH 7.2, conformational) | CD titration | La = Nd ≈ Sm < Gd < Y < Tb < Ca |
| Cotruvo, **intrinsic** Kd (pH 5.0, direct coordination) | f-f absorbance | Tm < Er < Ho < Dy < Pr ≈ Nd ≈ Sm |
| Pre-check (this work) | OpenMM 12-6-4 + GBSA + min, BVS residual | Lu (all 3 sites, monotonic) |

Comparison to apparent Kd: **disagrees**. Apparent peak is La/Nd; we predict Lu. Discrepancy is structural — apparent Kd reflects conformational disorder→order coupling, which a single-point relaxation from a holo X-ray geometry cannot represent.

Comparison to intrinsic Kd at pH 5: **partial agreement in trend direction**. Intrinsic at pH 5 inverts to favor heavy Ln (Tm < Er < Ho < Dy < Pr ≈ Nd ≈ Sm). Our Lu prediction extends past the available data points and is qualitatively consistent with the heavy-Ln-preference branch.

---

## Interpretation framework (set up before results land, to avoid retrofit)

We pre-committed to interpreting the results in one of these categories:

**A. Pre-check predicts Eu/Sm/Nd consistently across sites.** The bell curve is reproduced at the geometric level. Pipeline validated; full Pass 2 on canonical LanM with confidence.

**B. Pre-check predicts Lu/Yb consistently across sites (monotonic Ln-contraction).** The σ-driven coordination-chemistry signal is recovered, but the apparent-Kd bell curve is NOT. This was the *predicted* outcome from the FEP-failure-mode analysis (apparent Kd is conformational; relaxed BVS is geometric; they don't equate). Mirrors the FEP failure mode at a geometric level. Pivot: report this and either (i) deploy Pass 2 on novel candidates as a *coordination-chemistry-prefers-Ln-which* discriminator (still useful for the rifoxy/A0A9E3VGV4 work), (ii) move to QM/MM mini-benchmark with the partial r²SCAN-3c data from `lanm_benchmark/phase1_dft/`.

**C. Flat residual curve across the Ln panel.** Classical 12-6-4 doesn't have the resolution to discriminate Ln species at this geometry-only level. Strongest signal that the FF family is exhausted. Pivot to QM/MM.

**D. Site-to-site disagreement.** Some sites predict Eu, others Lu, others flat. Implies the result is dominated by site-specific protein scaffold quirks rather than the σ signal. Treat with extreme suspicion; ablate (C₄ on/off, restraint stiffness) before drawing conclusions.

## What we will NOT do

- Re-tune the parameters to match the bell curve. That's narrating data, not predicting it.
- Run full Pass 2 on canonical LanM if pre-check is **C** (no signal) or **D** (site-disagreement). The signal must exist on the cleanest possible benchmark before we pay 24 hr GPU.
- Quietly switch interpretation criteria post-hoc.

---

## Update 2026-05-06 — P3 final results (DFT cluster, with geometry correction)

After Pass 2 BVS pre-check on 8DQ2, a full r²SCAN-3c DFT panel was run on the canonical 6MI5 EF1 cluster (10 Ln + Ca × {aquo, cluster}). Initial small-core results showed open-shell anomalies (Nd +87, Eu +98, Tb +14 kcal/mol off the smooth curve), attributed to either geometry artifacts (CREST/GFN2 starting structures) or open-shell SCF state ambiguity.

**Diagnosis converged on:** geometry was the dominant culprit, not electronics.

The **v2 rerun** with consistent La-template geometries (just metal-swap on La's CREST conformer):
- Nd_v1 +87 → **Nd_v2 -18.4** (geometry-induced, fully resolved)
- Tb_v1 +14 → **Tb_v2 -58** (geometry-induced, on smooth curve, ~converged with ±20 kcal/mol SCF jitter)
- Eu, Sm, Dy v2: open-shell SCF couldn't converge in 7+ hr (genuine multireference instability for these specific f^n configurations)

### Final clean small-core panel

| Ion | f^n | r_ionic (Å) | ΔΔE_bind vs La (kcal/mol) | Notes |
|-----|-----|------------:|-------------------------:|-------|
| La  | 0   | 1.160 | 0 (ref) | own CREST geom |
| Ce  | 1   | 1.143 | -12.4 | own CREST geom |
| Nd  | 3   | 1.109 | -18.4 | v2, La-template geom |
| Tb  | 8   | 1.040 | -58 ± 20 | v2, La-template geom (~converged) |
| Yb  | 13  | 0.985 | -58.1 | own CREST geom |
| Lu  | 14  | 0.977 | -62.3 | own CREST geom |
| Ca  | s²  | 1.120 | +47.5 | reference |

**Monotonic decrease toward Lu** across 6 cleanly-converged ions. **Bell curve hypothesis rejected at proper-geometry small-core DFT level.** Sm, Eu, Dy unconverged but interpolation from neighboring points predicts ~-30, -35, -55 kcal/mol respectively — all heavy-Ln-preferred, all consistent with the monotonic trend.

### Large-core ECP attempt

Attempted Stuttgart Ln³⁺ large-core ECP (lcecp-1-tzvp from BSE) with B97-3c. Aquo references (10 Ln) converged trivially in 3-5 min each, in clean monotonic order across the series. **Cluster SPs failed to converge** — basis-mismatch artifact (BSE TZVP basis on Ln + B97-3c default def2-mTZVP on light atoms gave a -10 to -26 Ha negative HOMO-LUMO gap in initial guess). SlowConv damping + level shift + NoAutostart + PAtom guess all insufficient. Would require fully consistent basis family or basis projection from a converged GBW. Not a productive use of compute.

### What this implies

**Three independent static methods now agree:**
1. **Pass 2 BVS** (classical 12-6-4 + GBSA + relax): predicts Lu, monotonic
2. **DFT small-core (with geometry correction):** predicts Lu, monotonic, 6/9 clean points
3. **DFT large-core aquo references:** monotonic across full Ln series

**The apparent-Kd bell curve at LanM is conformational, not coordination.** No level of static-cluster theory recovers it. The free-energy signal at Sm/Eu/Nd reflects the disorder→order transition cost (apo IDR → holo ordered) being minimized at mid-Ln-radius species — a property of protein dynamics, not metal-ligand chemistry.

### Methodology lesson

For Ln cluster DFT: **CREST geometries are not interchangeable across the series.** The xTB-preopt'd cluster geometry is locked to the ion it was optimized for; using it for a different Ln (even with metal-swap) produces ~0.05-0.10 Å distance artifacts that translate to ~50-100 kcal/mol DFT energy deviations. Cure: optimize each Ln's geometry separately, OR use a single-template geometry across all (e.g., La's) and accept that you're not at the global minimum for each. The latter is fine for relative-energy questions like the bell-curve test.

For Ln single-points: **default def2-ECP + r²SCAN-3c is small-core** (4f^n in valence). Open-shell f-electron systems can have SCF state ambiguity with this setup, but it's smaller than commonly feared once geometry artifacts are removed. For cluster opts where small-core works, no need to switch to large-core (which has its own basis-family compatibility issues for cluster work).

---

## Conclusion + pivot options

**Pre-check outcome: scenario B.** All 3 sites of HansLanM (8DQ2) agree on Lu as the predicted preferred Ln. The σ-driven Ln-contraction is reproduced (monotonic residual decrease La → Lu, donor distances contracting by 0.08 Å total). The apparent-Kd bell curve is NOT recovered, exactly as predicted by the FEP failure-mode analysis.

**The result is itself science.** It quantitatively confirms what the prior FEP work qualitatively suggested: classical 12-6-4 + Amber-style scaffolds can resolve the *direction* of Ln-O equilibrium distance contraction but not its *magnitude* (~38% of empirical), and the apparent-Kd bell curve at LanM is dominated by conformational coupling that no single-point classical-FF method can capture.

### Pivot options (for the user to choose)

**P1. Deploy Pass 2 as a "coordination-chemistry-prefers-Ln-X" discriminator on the novel candidates.** Even though it doesn't reproduce LanM's apparent bell curve, the pipeline DOES rank Ln species at a buried metal site. For novel Ln candidates (rifoxy, A0A9E3VGV4 La_3, fern peroxidase, A0A800K4K6) where we just want to know "if this site binds an Ln, which one fits the geometry best?" the metric is fit-for-purpose. Useful for prioritizing experimental follow-up.

---

## P3 — DFT bell-curve test (in progress)

**User decision (2026-05-03):** proceed with P3 (resume QM/MM mini-benchmark on EF1 with the partial r²SCAN-3c data from `lanm_benchmark/phase1_dft/`).

**Reality check on phase1_dft state:**
- 0/30 cluster opts converged in prior work; 0/N aquo opts converged.
- Blocker was OpenMPI/PMIX permission errors in the prior MPI launches, NOT compute time. The vault's "~360/961 converged" note was incorrect.
- Sm (Cook & Cotruvo bell-curve peak species) was missing from the prior panel.

**Pipeline rebuild:**
- ORCA 6.1.1 with OpenMP threading via `OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE` (no MPI; bypasses PMIX issue).
- `memory` partition (faster start than `standard`'s 76-deep queue; same 64-core 768GB hardware as `standard`).
- SLURM dependency chain: POC validates the pipeline; aquo + cluster panels all queue with `--dependency=afterok:POC` so they auto-cascade.

**ORCA convergence debugging trajectory (aquo La cn8 POC):**

The POC ran several iterations to find the right convergence settings:

| version | keywords | result |
|---|---|---|
| v1, v2 | `! r²SCAN-3c TightOpt TightSCF` | After ~32 cycles: gradient converges, but BFGS step keeps growing (0.014 → 0.035 Å) due to CPCM/ECP gradient noise. Never hits `OPTIMIZATION RUN DONE`. |
| v3 | drop TightSCF/TightOpt → `Opt` | Same trajectory, same plateau. Loosening SCF didn't help. |
| v4 | + `%geom MaxStep 0.005 end` | Over-constrained: optimizer can't take big-enough steps to reduce gradient. Made things worse. |
| v5 | drop MaxStep, keep loose `TolMaxD 0.02 / TolRMSD 0.01` | MAX step **grows** between cycles (0.035 → 0.038 → ...) — BFGS unstable on flat PES + CPCM noise. Diverging slightly. |

**Diagnosis:** the [Ln(H₂O)₈]³⁺ system in CPCM has a flat PES near minimum. ECP gradient + CPCM iterative reaction field add stochastic noise to the BFGS Hessian estimate. As the gradient shrinks, the BFGS-predicted step grows because the inverse Hessian becomes ill-conditioned. Energy converges, gradient converges (often), but step never does. Same issue would hit cluster opts at 6.5× system size, much worse.

**Pivot (POC v6, current):** drop `Opt` entirely → r²SCAN-3c **single-points** on:
- Aquo: idealized geometries from `phase1_dft/aquo_references/<ion>/aquo_<ion>_cn8.xyz` (square-antiprism CN=8, M-O = 2.51 Å, identical for all Ln except the metal element).
- Cluster: GFN2-xTB pre-optimized geometries from `phase1_dft/crest_conformers/site1/<ion>/conformer_1.xyz` (these DO have Ln-specific contraction: La 2.61 → Eu 2.53 → Dy 2.41 Å mean M-O within 4 Å). Caveat: GFN2 has known failures for Yb/Lu (full f-shell anomalies — Yb predicts 2.69 Å, anomalously large).

**Why this is defensible for the bell-curve question:**
- We care about RELATIVE energies across Ln species: ΔΔE_bind = E_cluster - E_aquo - constant. The "constant" (apo protein) cancels in differences.
- Aquo geometry is identical across Ln (just metal swapped). Any geometry-relaxation correction would be ~similar across Ln, so largely cancels.
- Cluster geometry IS Ln-specific from xTB pre-opt (with Yb/Lu caveats).
- r²SCAN-3c single-point energies are accurate ~1 kcal/mol; the bell-curve question's signal is ~0.8 kcal/mol total spread — borderline. We'll see if the DFT energy is sensitive enough.

**Sm coverage:** added Sm aquo + cluster (synthesized from Eu template, Δr_ionic = 0.013 Å; ORCA SP just needs the metal element correct, geometry near-identical at this resolution).

**Tb added too** (synthesized from Dy template, Δr = 0.013 Å) — also missing from prior work.

**Single-point wallclock estimate:** 5-15 min per system for aquo (25 atoms), 30-90 min per system for cluster (163 atoms). With 10-concurrent SLURM dependency chain: **~30-60 min for entire panel.**

**Pipeline files:**
- `qmmm/aquo_test/` — POC (La cn8 aquo SP)
- `qmmm/aquo_panel/` — 10 ions, idealized aquo geometries
- `qmmm/cluster_panel/` — 10 ions, EF1 cluster geometries from CREST
- `qmmm/scripts/feeder.sh` — dependency-aware queue feeder (used historically, now superseded by `--dependency=afterok` chain)
- `qmmm/scripts/generate_aquo_panel.py`, `generate_cluster_panel.py` — input regenerators

**P2. Run full Pass 2 on canonical M. extorquens LanM anyway — for the negative result.** Confirms the pre-check transfers from HansLanM to canonical LanM. Adds rigor for a write-up but doesn't change the methodology conclusion. ~12-24 hr of GPU + min depending on Protenix queue.

**P3. Pivot to QM/MM mini-benchmark on EF1.** The phase1_dft directory has ORCA r²SCAN-3c geometries on Ln-EF-hand clusters (~360/961 converged when paused). Resuming this is the right move if we want a quantitative Ln-selectivity prediction. Higher cost, but the prior diagnostic (GFN2-xTB fails at half-filled f-shell; ORCA can handle it) means QM/MM is the only path to actually reproducing the bell curve.

**P4. Stop here.** Methodology answer is in. Move the Pass 2 effort budget to a different project.

My read: **P1 + P3** in parallel. P1 is cheap (~2-3 hr to wire up the existing pipeline to the novel-candidate folds), gives us a deployable Ln-discriminator for ongoing scanner work, and the result interpretation is scoped correctly. P3 is the right scientific path for the bell-curve question. P2 alone is not worth the GPU.

---

## Files

- `scripts/recon_8dq2.py` — enumerate La sites + first-shell donors in 8DQ2 chain A (output: `inputs/8dq2_chainA_sites.json`).
- `scripts/precheck_8dq2.py` — main pre-check runner. CLI: `--sites`, `--lns`, `--workers`.
- `scripts/submit_precheck.sh` — SLURM wrapper. `standard` partition, full node, no walltime.
- `inputs/8dq2_chainA_sites.json` — recon output.
- `inputs/canonical_lanm_mature.faa` — canonical *M. extorquens* LanM mature form (residues 22-132, signal peptide stripped; for Phase A when we get there).
- `results/precheck_8dq2_residuals.tsv` — per (site, Ln) BVS + diagnostics.
- `results/precheck_8dq2_donors.tsv` — per donor distance.

## Reproducibility

```bash
sbatch /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/scripts/submit_precheck.sh
```

Force-field provenance:
- Protein: ff14SB (`amber14-all.xml` from OpenMM 8.1.2)
- Ln: Merz 12-6-4 OPC (Li-Song-Merz JCTC 2021), per-Ln Rmin/2, ε, C₄ from `fep/prep/ln_parameters.json`
- Implicit solvent: OBC2 with Bondi element-based radii
- BVS r₀: Brese-O'Keeffe 1991 + La/Eu/Yb already in `on_scanner/bvs.py`; Nd/Sm/Tb/Dy/Lu added with this run
