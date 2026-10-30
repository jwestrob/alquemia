# La/Dy ALMO charge-transfer quartet — 2026-10-24

Question: at the archived fixed geometries, does suppressing metal↔ligand-shell
charge transfer materially change the La/Dy water→acetate double difference
C = [X(Dy,mix) − X(La,mix)] − [X(Dy,aq) − X(La,aq)]?

Supersedes the stalled gating in `../lady_almo_implementation_20261023/` (kept
intact). Changes in approach, per Jacob 2026-10-24:

- CT = U − P is computed entirely within one backend (PySCF PBE0, direct J/K,
  unpruned 150×974, same basis/ECP/geometry as the 10-23 baselines). Absolute
  agreement with archived ORCA RIJCOSX energies is **not** a gate for CT; the
  PySCF−ORCA offset is reported as its own double-difference line.
- Mixed-Dy baseline (gradient 1.067e-7 vs 1e-7) is polished by restarting SCF
  from its saved density (`--polish`) instead of being treated as a blocker.
- CP2K validation dropped. Solver check = `scfmi_check.py`: almo.py versus an
  independent finite-difference L-BFGS over occupied→virtual mixings (no shared
  derivative/retraction/normalization code), water dimer and Mg²⁺·H₂O, HF and
  PBE0, plus a 10 Å-separated dimer (CT must vanish).

## Acceptance criteria (declared before any quartet CT number was seen)

Solver: |P_A − P_B| ≤ 1e-6 Eh in every scfmi case; U ≤ P ≤ F; far-dimer |CT| < 1e-3 kcal/mol.
Endpoint: ALMO converged (|g| ≤ 1e-6, |ΔE| ≤ 1e-9 Eh), U ≤ P ≤ F; U replay of saved
baseline within 1e-8 Eh (or polish converged with |g| ≤ 1e-7).
Interpretation: CT is called **material** for the selectivity question if
|C_CT| ≥ 0.5 kcal/mol **and** ≥ 20 % of |C_full|; **negligible** if |C_CT| < 0.2
kcal/mol; in between = minor. The backend-offset double difference is reported
as a scale for numerical noise in the full differential.

## Limits

Conditional gas-phase cluster energies at fixed geometry; not binding free
energies. Effective-singlet ECP model: Dy 4f electrons are in the ECP55 core,
so 4f covalency is not represented (a null CT result does not exclude it).
D4 is the archived scalar, added identically to full and constrained states.
AO-block ALMO with def2-TZVP (no diffuse functions); first-order stationarity only.

## Files

`ct_run.py` endpoint driver · `scfmi_check.py` solver cross-check · `collect.py`
double differences · `job.sh` sbatch wrapper. Outputs:
`workspaces/lady_almo_ct_20261024/`. Jobs: 23037 (scfmi), 23038–23041 (quartet).

## Result (all gates passed; jobs 23037–23041, 2.0 h wall total, longest 47 min)

Solver cross-check: all 5 cases |P_A − P_B| ≤ 6.3e-13 Eh (criterion 1e-6); far
dimer CT = 3.6e-11 kcal/mol. Water dimer CT −1.84 (HF) / −3.56 (PBE0); Mg²⁺·H₂O
CT −8.13 / −13.13 kcal/mol. U replays ≤ 8e-12 Eh; mixed-Dy polish converged
(|g| 8.2e-9, ΔE −3e-12 Eh vs saved). All four ALMO solves converged in 13–24
iterations (|g| ≤ 9.7e-7), U ≤ P ≤ F everywhere.

Balanced double differences (kcal/mol; negative = Dy favoured by water→acetate):

| quantity | C |
|---|---|
| full PBE0 electronic (PySCF) | +0.447 |
| full PBE0 electronic (archived ORCA) | +0.457 |
| backend offset PySCF − ORCA | **−0.010** |
| full PBE0 + D4 | +0.565 |
| CT-suppressed (P) + D4 | +0.979 |
| frozen (F) + D4 | +0.558 |
| ΔFRZ / ΔPOL / ΔCT | +0.440 / +0.420 / **−0.414** |

Per-endpoint CT: La aq −68.10, Dy aq −63.25, La mix −70.52, Dy mix −66.08.
Acetate raises |CT| by 2.41 (La) vs 2.83 (Dy) kcal/mol.

Verdict under the pre-declared criteria: |C_CT| = 0.41 kcal/mol → **minor**
(0.2–0.5 band), although it is 73 % of |C_full|. CT favours Dy on acetate
substitution; removing it nearly doubles the La preference (+0.56 → +0.98) but
does not change its sign. The ORCA/PySCF backend difference cancels to
0.01 kcal/mol in the double difference, so the old 1e-4 Eh per-endpoint parity
gate was irrelevant to this quantity.

## Counterpoise (BSSE) bound — jobs 23067–23074, `cp_run.py`

Each fragment re-converged in the full complex basis (partner atoms as ghosts).
BSSE_A = E_A(full basis) − E_A(own); CT_cp = CT − ΣBSSE (Khaliullin 2007).
Results: `workspaces/lady_almo_ct_20261024/CP_RESULT.json`.

| endpoint | BSSE metal | BSSE shell | CT | CT_cp |
|---|---|---|---|---|
| La aq | −0.009 | −10.08 | −68.10 | −58.02 |
| Dy aq | −0.015 | −8.16 | −63.25 | −55.07 |
| La mix | −0.010 | −10.83 | −70.52 | −59.68 |
| Dy mix | −0.016 | −8.74 | −66.08 | −57.33 |

BSSE is almost entirely the ligand shell borrowing metal functions (the La
basis is more diffuse), and it does not cancel: ΔΔBSSE = +0.18 kcal/mol.
**ΔΔCT_cp = −0.60 kcal/mol** (vs −0.41 uncorrected); CT-suppressed + D4
differential with CP accounting = +1.16 kcal/mol.

## Final verdict

ΔΔCT lies between −0.41 (raw) and −0.60 (counterpoise) kcal/mol, which is right
at the pre-declared 0.5 "material" threshold, and is 73–106 % of |C_full| = 0.56.
Robust statements: CT favours Dy on water→acetate substitution; it is the same
magnitude as the whole differential; suppressing it roughly doubles the La
preference (+0.56 → +0.98–1.16) without changing the sign. Classification:
borderline minor/material, depending on BSSE treatment.
