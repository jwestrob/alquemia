# Acetate La→Dy forensics — PLAN (frozen 2026-11-10 00:50 UTC, before any recomputation)

Question: why do static cluster-continuum DFT (−1.95) and classical 12-6-4 FEP (−2.29) both make
acetate ~2 kcal/mol too Dy-selective against NIST (−0.07, I = 0.1), and which single cheapest
calculation would decide between the leading explanations?

Scope: saved outputs only (login node, single thread). No Slurm, no new QM/MD/MLP calculation.
Outputs: scripts and REPORT.md here; tables in `workspaces/lady_acetate_forensics_20261109/`.

## Step 0 — positive control (stop if it fails)

Using the existing code (`lady_relaxed_exchange_20261024/analyze.py` functions, imported, never run as
`main()` because main() overwrites saved JSON):
- acetate pooled `ac`, conformer_best, lowest-G water-balanced rule: must equal −1.95 ± 0.01.
- ΔΔG_hyd(Dy−La) with the `collect.py` formula (conformer-best La aq9 / Dy aq8 G, gas bare-ion E,
  one water + WATER_CORR): must equal the saved `RESULT.json` −59.21 ± 0.01. The brief's −59.4 appears only
  in the 10-29 review document; record which value reproduces.

## Hypotheses and pre-stated evidence criteria

Sign convention: exchange = complex leg − aqua leg, legs are Dy − La; error = calc − exp. Negative error =
too Dy-selective. Experimental aqua leg −66.9 (Marcus via Li/Song/Merz, from EXPERIMENT_THERMO.md).

**(a) Aqua-ion deficit that fails to cancel.**
Support: each model's aqua-leg error has the sign that produces Dy over-selectivity (calc ΔΔG_hyd less negative
than −66.9) in BOTH models; a knob that moves ΔΔG_hyd (the saved CPCM metal-radius variants) moves the acetate
exchange in the same direction, and extrapolating that knob to ΔΔG_hyd = −66.9 brings acetate within 0.5 of −0.07.
Against: the two models' aqua-leg errors have opposite signs (then a shared aqua cause is excluded), or the
radius knob moves ΔΔG_hyd several kcal/mol while acetate moves < 0.5.

**(b) Carboxylate–Ln overbinding that grows with charge density.**
Support: the acetate error sits in the complex leg in both models; at fixed geometry DLPNO-CCSD(T) − PBE0-D4
(saved `dlpno`/`gas` single points, n = 0 combinations only, because the DLPNO water job failed) is ≤ −0.5
kcal/mol on the acetate exchange (PBE0 over-favours Dy); the error grows with the number of bound carboxylate O.
Against: DLPNO − PBE0 within ±0.3 or positive; FEP complex leg error near zero.

**(c) Binding mode and hydration-state treatment.**
Support: the acetate ΔΔG spread across the treatments below exceeds 1 kcal/mol AND physically reasonable
treatments (Boltzmann, matched sets) move toward −0.07 by ≥ 1 kcal/mol; or the La/Dy candidate-set asymmetry
can produce ≥ 0.5 kcal/mol of selection bias given the measured conformer noise.
Against: all ensemble-type treatments stay within ±0.5 of −1.95.

**(d) Experimental side (ionic strength, what K1 contains, inner/outer sphere).**
Support: ΔΔG differs by ≥ 1 kcal/mol between ionic strengths (I = 0.1 NIST vs 2 M Grenthe), or an outer-sphere
share large enough that the inner-sphere ratio could be ~e^(2/RT) while the total ratio is ~1.
Against: the two ionic strengths agree within 0.3, and a Fuoss outer-sphere estimate gives K_os ≪ K1
(then K1 is inner-sphere-dominated and the measured ratio is the inner-sphere ratio).

**(e) Thermochemistry / standard-state conventions.**
Support: the selected acetate combination releases n ≠ 0 waters (so the 1 atm→1 M and 55.34 M terms enter),
or the thermal (G − E) part of the exchange is ≥ 1 kcal/mol.
Against: n = 0 for the selected combination and |G − E| exchange < 0.5. A secondary check: shift the per-water
reference free energy by δ ∈ [−4, +4] kcal/mol (covers the monomer- vs cluster-cycle conventions) and record how
the selected states and the Boltzmann exchange move; this is where (e) couples to (c).

## Analyses

1. **State table** (`state_table.py`): every acetate and aqua state for both metals, primary and conformer-best:
   E, CPCM dielectric term, D4, H, G, −TS (= G − H), coordination number (donors < 2.95 Å), bound carboxylate O,
   imaginary modes, selected or not. Exchange under:
   (i) lowest-G rule; (ii) Boltzmann sum over water-balanced states per metal and family (same WATER_CORR);
   (iii) matched state sets (only state labels present for both metals); (iv) every fixed state quadruple.
   Report min/max/spread. Selection effect: per-state conformer noise σ from the saved conformer candidates
   (G spread among valid candidates), number of candidates per metal, and the expected-minimum bias
   σ·[c(k_Dy) − c(k_La)] with c(k) = E[max of k standard normals].
2. **Leg decomposition** (`legs.py`): for ac, edtax, nta, dtpa, cdta, tmdta (malonate/TTHA unfinished):
   aqua-leg error e_aq = calc ΔΔG_hyd − (−66.9); complex-leg error = exchange error + e_aq.
   Inheritance models, form fixed here, fitted by least squares, LOO RMS reported (descriptive, n = 6):
   - M0: bias = a
   - M1: bias = a + b·n_disp (mean over metals of aqua waters − complex waters, selected states)
   - M2: bias = a + b·(9 − q_c²) (continuum Born scaling between aqua ion and complex)
   - M3: bias = a + b·n_Ocarb (mean over metals of bound carboxylate O, < 2.95 Å, from the xyz)
   Zero-parameter inheritance of the aqua error: Z1 bias = −e_aq·n_disp/n_aq (n_aq = mean aqua waters, 8.5);
   Z2 bias = −e_aq·(1 − q_c²/9). Pass for "inheritance explains it": Z RMS < 1.0.
   Plus exchange decomposition per ligand (E, D4, CPCM, rest, G − E) and the radius knob (acetate + aqua).
3. **Thermo split** for acetate: model ΔΔH, TΔΔS vs Grenthe 2 M; per-leg ΔΔH_hyd vs −69.8.
4. **FEP**: legs with and without the analytic metal–OW tail (the 12-6-4 fit convention for that tail is not in
   the saved outputs; report both), binding-mode populations, hydration numbers, where the error sits.
5. Rank (a)–(e); state where the saved outputs cannot discriminate.

## Design (not launched)

Cost, from actual workspace timings, at least: (1) explicit second-shell cluster convergence (aqua ions +
acetate); (2) matched-pair vertical swap at fixed geometry + harmonic relaxation from existing Hessians vs the
independently optimized endpoints; (3) small-box DFT-MD thermodynamic integration. Each with: what, which
systems, predicted outcome per hypothesis, pass/fail threshold, CPU/GPU-h, and how it could be uninformative.
