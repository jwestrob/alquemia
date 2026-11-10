# Second-shell test, stage A — REPORT (10 November 2026)

**Verdict: FAIL under the frozen gate. Stage B stays unreleased and the test ends here.**

An explicit second hydration shell moves ΔΔG_hyd(Dy−La) by −2.38 kcal/mol (seeds −1.35 and −3.41). The gate
required −4.0. The second-shell value is −59.92 against −66.9 experimentally, so 6.98 of the 9.36 kcal/mol
matched nine-coordinate deficit remains. All values below are in kcal/mol.

## Gate (frozen analysis)

`second_shell_analyze.py`, output `workspaces/lady_second_shell_20261109/GATE.json`.

| Seed | La CN | Dy CN | ΔΔG_hyd, second shell | First-shell reference (CN 9) | Shift |
|---|---|---|---|---|---|
| s1 | 9 | 9 | −58.89 | −57.54 | −1.35 |
| s2 | 9 | 9 | −60.95 | −57.54 | −3.41 |
| mean | | | −59.92 | | −2.38 |

- Both seeds are admissible: optimization converged, frequencies done, Dy nine-coordinate, the 2.95–3.10 Å
  band empty. Dy s1 has one imaginary mode at −9.9 cm⁻¹ (the rule excludes modes below −30).
- Rule: mean −2.38 > −4.0 and both shifts > −4.0, so FAIL. Seed difference 2.06.
- `PLAN.md` and `second_shell_analyze.py` match `FREEZE_PLAN.sha256` and `FREEZE_ANALYZE.sha256`, checked
  after the run.
- Jobs: 69800 (6 h 31 min, 40 optimization cycles) and 69802 (12 h 19 min, 120 cycles) on 24 CPUs each;
  the restart jobs 69801 and 69803 had nothing to do. Total 452 CPU-h.

## Where the shift sits (descriptive, post hoc)

`stage_a_descriptive.py`, output `workspaces/lady_second_shell_20261109/STAGE_A_DESCRIPTIVE.json`. Written
after the verdict was read; it takes no part in the gate.

Dy-minus-La hydration difference by term of the ORCA thermochemistry (E_el includes CPCM and D4; the
bare-ion energies enter E_el; enthalpy = E_el + ZPE + thermal):

| | E_el | ZPE | thermal | enthalpy | −TS | G |
|---|---|---|---|---|---|---|
| First shell, CN 9 pair | −63.48 | +2.78 | −2.25 | −62.94 | +5.40 | −57.54 |
| Second shell, seed 1 | −61.66 | +1.32 | −1.07 | −61.40 | +2.51 | −58.89 |
| Second shell, seed 2 | −63.59 | +1.63 | −0.78 | −62.73 | +1.78 | −60.95 |
| Second shell, mean | −62.62 | +1.48 | −0.92 | −62.07 | +2.15 | −59.92 |
| Shift, mean of seeds | +0.85 | −1.30 | +1.32 | +0.87 | −3.25 | −2.38 |
| Experiment | | | | −69.8 | +2.9 | −66.9 |

The experimental enthalpy difference −69.8 is the approximate value quoted in
`diagnostics/lady_acetate_forensics_20261109/REPORT.md`; +2.9 is its difference from −66.9.

- **The whole shift is in the entropy term.** −TS moves from +5.40 to +2.15, close to the experimental +2.9
  (error +2.50 with the first shell, −0.75 with the second).
- **The enthalpy stays where it was.** It moves by +0.87 (seeds +1.54 and +0.21), away from experiment. Its
  error is +6.86 with the first shell and +7.73 with the second.
- The forensics report had placed the aqua deficit in the enthalpy. The second shell leaves that part
  untouched.

Structures of the two pairs:

| | La s1 | Dy s1 | La s2 | Dy s2 |
|---|---|---|---|---|
| Mean M–O, first shell (Å) | 2.578 | 2.448 | 2.592 | 2.462 |
| Range (Å) | 2.542–2.637 | 2.407–2.538 | 2.546–2.673 | 2.400–2.564 |
| Mean M···O, 12 second-shell waters (Å) | 4.591 | 4.503 | 4.496 | 4.412 |
| O···O contacts below 3.2 Å | 38 | 40 | 37 | 42 |
| Modes below 100 cm⁻¹ | 20 | 18 | 20 | 17 |

- Each Dy cluster keeps the arrangement of its La partner: the same nine first-shell waters, 12 second-shell
  waters, and every O···O contact of the La structure (2 and 5 contacts are added as the cluster contracts).
  RMSD between the La and Dy final structures: 0.157 Å (s1) and 0.380 Å (s2), all atoms.
- First-shell-only references: La_aq9 2.561 Å, Dy_aq9 2.429 Å.
- Seed 2 minus seed 1 for the same metal: La +3.64, Dy +1.58 in G. Pairing leaves 2.06.

## Consequences

- Stage B (acetate) stays unreleased, by the stop rule. With the measured aqua shift, hypothesis (a) predicts an
  acetate shift of 0.24 × 2.38 = +0.57, inside the plan's inconclusive band (+0.3 to +0.7).
- In the plan's wording: an explicit second shell at this level does not carry the aqua deficit.
- The static cluster-continuum route keeps a Dy−La hydration error of about +7, with its enthalpic part (7 to 8)
  unchanged by the second shell. The forensics report found the same error insensitive to the CPCM metal
  radius (−0.91 to +0.35).

## Limits

- Two seeds, differing by 2.06 in the shift. The forensics proposal named a spread above 0.5 as a failure mode.
- Both Dy clusters stay nine-coordinate, so the eight-coordinate Dy state with a second shell is unsampled.
  The swap starts from the nine-coordinate La arrangement; the reference rule covered a change that did not
  occur.
- The term-by-term split rests on ORCA's default quasi-harmonic thermochemistry for clusters with 17 to 20
  modes below 100 cm⁻¹. It is a description of this calculation.
- The experimental enthalpy difference is flagged approximate in the forensics report and was not re-derived.

## Unverified leads (nothing started)

The remaining enthalpic error is common to the first-shell and second-shell models. That points at an
ingredient both share: the metal description (4f-in-core pseudopotential and its basis), or the reference
value itself (single-ion hydration values come from thermochemical cycles that include the first three
ionization energies of each metal). Both are hypotheses; neither was examined.

## Files

- `PLAN.md` (frozen), `DEVIATIONS.md`, `CURRENT.md`, `second_shell_analyze.py` (frozen),
  `stage_a_descriptive.py` (post hoc).
- `workspaces/lady_second_shell_20261109/`: `GATE.json`, `analyze_final.log`, `STAGE_A_DESCRIPTIVE.json`,
  `watch.log`, `runs/Dy_ssaq_w21_s{1,2}/`.
