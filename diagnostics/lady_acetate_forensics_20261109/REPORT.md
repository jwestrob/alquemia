# Acetate forensics: why two models miss La→Dy acetate — 9 November 2026

**Result.** The shared sign comes from two different mechanisms. The static DFT misses acetate mainly
through its aqua leg; the force-field FEP misses it in its complex leg. Coupled cluster agrees with
PBE0-D4 on the metal swap, so the functional is cleared.

Written by the lead agent from the analysis subagent's reply and the saved tables (see
DEVIATIONS.md). `PLAN.md` sha256 `8bf3c02a…` was frozen first. Convention: error = calc − exp, negative
= too Dy-selective; kcal/mol.

## Positive control

Acetate reproduces at −1.947 (La_acmono_w8, Dy_acbi_w7, La_aq9, Dy_aq8; zero net waters released).
ΔΔG_hyd(Dy−La) is −59.21 from conformer-best structures and −59.35 from primary structures.

## Acetate across state treatments

| Treatment | ΔΔG |
|---|---|
| Lowest-G rule | −1.95 |
| Boltzmann over states | −2.02 |
| Matched state labels | −1.95 / −2.02 |
| Monodentate only | −1.44 |
| Bidentate only | −2.58 |
| Per-water reference shifted ± 4 | −1.38 to −2.09 |

- Same-type quadruples range −7.07 to +3.38; the spread comes from mismatched aqua coordination
  numbers.
- La_acmono_w8 is the only selected state never conformer-searched. Searches on the other complex
  states lowered G by 0 to 0.52, so the expected bias is about +0.15.
- La's 10-coordinate bidentate state (La_acbi_w8) was never computed. The FEP puts La at
  10-coordinate about 1% of the time.
- DLPNO-CCSD(T) − PBE0-D4 at fixed geometry (from saved single points) is about 0 for matched-mode
  exchanges (−0.16, +0.09, −0.12). It stabilises bidentate over monodentate by 2.30 (La) and 2.34
  (Dy), a metal-independent shift. Re-selecting states with that correction makes La bidentate and
  gives about −2.7 ± 0.4, assuming a common per-water offset (observed spread 0.64).

## Leg decomposition (n = 6; malonate and TTHA unfinished)

DFT aqua-leg error: **+7.69**.

| Ligand | ΔΔG calc | NIST | Error | Complex-leg error |
|---|---|---|---|---|
| acetate | −1.95 | −0.07 | −1.88 | +5.81 |
| EDTA | −6.35 | −4.01 | −2.34 | +5.35 |
| NTA | −4.60 | −1.58 | −3.02 | +4.67 |
| DTPA | −8.33 | −4.56 | −3.78 | +3.92 |
| TMDTA | −10.45 | −4.87 | −5.58 | +2.11 |
| CDTA | −11.66 | −5.67 | −5.98 | +1.71 |

- Acetate's error is the 24% of the aqua error that its complex leg does not cancel.
- Leave-one-out RMS of the inheritance models: constant 1.86; per water displaced 1.61;
  9 − q_c² 2.11; bound carboxylate oxygens 1.79. The zero-parameter forms score 1.84 and 2.90. None
  explains the pattern. Waters displaced and bound carboxylate oxygens correlate at r = 0.99, so the
  two cannot be separated.
- Acetate exchange terms: ΔE −3.27, D4 −0.26, CPCM +19.13, rest −22.14, thermal +1.32.
- CPCM metal-radius variants move acetate by −0.53 to +1.02 and ΔΔG_hyd by −0.91 to +0.35, so the two
  are not coupled through that setting.

## Thermodynamics

Against Grenthe at 2 M, the acetate ΔΔH error is −3.51 and the TΔΔS error is −1.71; they partly cancel
to −1.78. Enthalpy carries the error. The aqua-leg deficit is also enthalpic: ΔΔH_hyd is about −54.9
calculated against −69.8 (approximate water reference).

## FEP

- With the analytic tail the aqua leg is −70.16 (error −3.26); without it, −66.92, matching −66.9.
  That suggests the 12-6-4 fit excludes the tail (unverified).
- In both conventions the acetate error sits in the complex leg (−2.25 without tail, −5.48 with).
- Structure: La 95% monodentate, coordination 9.01, 1.04 waters displaced. Dy 89% bidentate,
  coordination 8.93 (aqua 8.71), 1.67 waters displaced. The DFT shows the same La-monodentate /
  Dy-bidentate split.

## Ranking of hypotheses

1. **(a) Uncancelled aqua-ion deficit: the DFT route.** The aqua error is large and enthalpic, acetate
   is its incomplete cancellation, and DLPNO clears the functional. The leg reasoning is inference; the
   radius setting is evidence against.
2. **(b) Carboxylate overbinding: the FEP route.** The error is in the complex leg with a fitted aqua
   leg. For DFT it is refuted electronically (≤ 0.2 on matched exchanges).
3. **(c) Binding mode and hydration state: secondary.** Physical treatments span −1.4 to −2.7; the
   binding-mode mismatch contributes about 0.5.
4. **(e) Standard-state conventions: excluded.** n = 0, so the terms cancel.
5. **(d) Experimental side: excluded.** 0.1 M versus 2 M differ by 0.10. A Fuoss estimate gives
   K_os ≈ 5 M⁻¹ against K1(La) ≈ 63 (inference), an effect of ≤ 0.1.

The saved outputs leave four things open: whether the DFT complex leg should keep more of the aqua
deficit than the 76% it keeps; the FEP tail convention; whether Dy-acetate is bidentate in reality;
and La_acbi_w8.

## Proposed decisive test (designed, not launched)

**Matched-arrangement second shell, staged**, on the matched 9-coordinate monodentate quadruple
(first-shell value −1.97, DLPNO correction +0.09, n = 0).

- **Stage A, aqua gate.** Build Dy_ssaq_w21 by swapping La→Dy in the finished La_ssaq_w21 s1 and s2,
  then Opt+Freq. The La originals took 18.2–19.9 h on 24 CPUs: budget ≤ 1,000 CPU-h. Gate: ΔΔG_hyd
  moves at least 4 kcal/mol toward −66.9 from the first-shell value of −57.55. Below that, stop.
- **Stage B.** Two matched seed pairs for La_ssacmono_w20 and Dy_ssacmono_w20 (Dy built by swap), 4
  jobs, ≤ 2,000 CPU-h. Under (a) acetate shifts by about 0.24 × the aqua shift, at least +1.0. Under
  (b) the shift is about 0. Pass for (a): shift ≥ +0.7 and shift / |aqua shift| ≥ 0.15. Against (a):
  shift < +0.3. Between: inconclusive.
- **Failure modes.** Seed-to-seed spread above 0.5; the third shell and continuum remain.
- **Total.** ≤ 3,000 CPU-h.

Other candidates:
- **Vertical swap plus harmonic relaxation.** Hessians exist; 8 swapped gradients cost about 10–20
  CPU-h. It tests whether −1.95 is an arrangement artefact and leaves (a) versus (b) open. A cheap
  control.
- **DFT-MD thermodynamic integration.** GGA in CP2K is about 13,000 CPU-h for 5 windows × 10 ps × 2
  systems, against nanosecond water exchange. Expected inconclusive.
- **La_acbi_w8.** One Opt+Freq (10–200 CPU-h) closes the missing-state question.

## Side observation

The DTPA and TMDTA exchanges carry thermal terms near +11 (DTPA +11.05), against +1.32 for acetate.
Their panel states' frequencies deserve a check.

## Limits

Unverified: the FEP tail convention; the Fuoss inputs; the water enthalpy reference in ΔΔH_hyd; whether
Dy-acetate is bidentate in reality.

## Files

- `diagnostics/lady_acetate_forensics_20261109/`: `PLAN.md`, `state_table.py`, `legs.py`,
  `dlpno_isomers.py`.
- `workspaces/lady_acetate_forensics_20261109/`: `state_table.tsv`, `state_treatments.json`,
  `legs.json`, `dlpno_isomers.json`, logs.
