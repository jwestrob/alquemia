# Experimental La→Dy thermodynamics vs the static cluster model (2026-10-26)

Dy − La, kcal/mol. PDFs in workspaces/lady_relaxed_exchange_20261024/src_untrusted/lit/.

| system | source | ΔΔG | ΔΔH | TΔΔS | model ΔΔG / ΔΔH / TΔΔS |
|---|---|---|---|---|---|
| acetate K1 | Grenthe, Acta Chem Scand 18 (1964) 283, Table 3; 25 °C, 2 M NaClO4 | −0.17 | +0.74 | +0.89 | −1.95 / −2.77 / −0.82 |
| EDTA | NIST SRD 46 v8 (Martell–Smith), 25 °C, I = 0.1: log K 15.36/18.30, ΔH −2.9/−1.2 | −4.01 | +1.70 | +5.71 | −6.35 / −6.87 / −0.52 |
| EDTA (check) | Ots 1973 ΔH, 1 M, Dy interpolated Tb/Ho; Dyatlova 1988 via PMC8582344 | −3.82 | +0.94 to +1.70 | | |
| acetate K1, I = 0.1 | NIST 46: log K 1.80/1.85 (no ΔH at low I) | −0.07 | | | |
| LBT3 ITC | Hatanaka, Sci Rep 10, 19468 (2020), SI Table S1; 10 °C, MES pH 6 | −2.36 | −1.59 | +0.77 | — |
| hydration | ΔG: Marcus 1991 via Li/Song/Merz JPCB 2014; ΔH, ΔS: computed from NBS 1982 tables (+R ln16 for Dy3+(g)) | −66.9 (−67.05) | −69.8 (−68.5 with NIST ASD IEs) | −2.74 | ΔG large-core cluster+CPCM: −59.2 |

Findings:
- Acetate and EDTA selectivity is entropic with enthalpy *opposing* Dy. The model gets ΔΔH with the wrong sign
  (error −3.5 acetate, −8.6 EDTA) and TΔΔS too low (−1.7, −6.2). The ΔG errors are similar (−1.8, −2.5) only because
  the two errors partly offset, so the "ligand-independent offset" hypothesis in README is refuted at the ΔH level.
- The enthalpy error grows with the number of first-shell waters displaced (acetate 1–2, EDTA ~6), and the model
  under-hydrates Dy relative to La by ~8 kcal/mol (ΔΔG_hyd −59 vs −66.9). Working hypothesis: desolvating Dy is
  too cheap in the model, by roughly 1.3 kcal/mol per displaced water.
  Falsifier: a model variant that reproduces ΔΔG_hyd = −66.9 but leaves the exchange ΔΔH errors unchanged.
- LBT3 differs: about two-thirds of its selectivity is enthalpic.

- Differential hydration itself is enthalpic (ΔΔH_hyd ≈ −70, TΔΔS_hyd ≈ −2.7). Complexation entropy is then the
  release of the more strongly bound Dy hydration water. NIST 46 extracts: src_untrusted/lit/nist46_*_ln.txt
  (decoded from the NIST installer data by tools/nist46.py; glycine and EDTA protonation constants reproduced).
