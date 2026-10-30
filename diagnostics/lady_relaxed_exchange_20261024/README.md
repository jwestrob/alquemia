# Relaxed La/Dy exchange ladder — started 2026-10-24 (interim)

Question: once La and Dy complexes are relaxed separately (instead of the archived
fixtures with every M–O at 2.500 Å for both metals), does a cluster-continuum model
recover the known La→Dy selectivity on small ligands?

Benchmarks for ΔG of La·L + Dy(aq) → Dy·L + La(aq), negative = Dy favoured (kcal/mol):
- acetate: flat within about ±0.3 (recalled Smith–Martell log K₁ is about 1.8 for La, peaks mid-series near Sm,
  then falls for the heavy Ln, so Dy is close to La; **not checked against a primary table**)
- EDTA, about −3.9 (log K 15.5 La vs about 18.3 Dy; the Dy value is unverified)
- LBT3 ITC, −2.36 (the target)

## Level of theory and accounting

ORCA 6.1.1 PBE0-D4/def2-TZVP, def2/J RIJCOSX, DefGrid3, CPCM(water), Opt + Freq. La ECP46 and
Dy ECP55 use the archived basis blocks (4f in core, effective singlets). Each released water gets
RT ln 24.46 + RT ln 55.34. `analyze.py` runs every hydration-number combination;
`analyze.py --conformers` replaces each state by its lowest-G structure. That structure is taken
from the primary run plus the CREST-conformer re-optimizations, and must keep the coordination
number and have no imaginary mode below −30 cm⁻¹.

## Result 1 — aqua ions (final for first shell)

| | ⟨M–O⟩ CN8 | ⟨M–O⟩ CN9 | ΔG(8→9) | ΔH(8→9) |
|---|---|---|---|---|
| La | 2.526 | 2.561 | **−3.41** | −9.91 |
| Dy | 2.378 | 2.429 | **+0.77** | −8.02 |

This reproduces the known switch, with La preferring CN 9 and Dy sitting just on the CN 8 side.
The contraction at each metal's preferred CN is 0.18 Å. The hydration-number asymmetry alone
contributes about 4 kcal/mol to the reference state, which the fixed-geometry vertical swap
could not represent.

## Result 2 — acetate (interim: first shell + CPCM, conformer-corrected)

| ligand | ΔG | ΔH | states (La·L, Dy·L, La aq, Dy aq) |
|---|---|---|---|
| acetate mono | −1.44 | −2.36 | La_acmono_w8, Dy_acmono_w7, La_aq9, Dy_aq8 |
| acetate bi | −2.58 | −10.13 | La_acbi_w7, Dy_acbi_w7, La_aq9, Dy_aq8 |

Dy is over-favoured by 1.1–2.9 kcal/mol relative to the roughly flat benchmark. That conclusion
holds for any experimental value within ±0.3. Conformer selection moves single states by up to 1.1 kcal/mol, the
same order as the signal, so a single-structure cluster number is not interpretable at this
precision. La_acmono_w8 and Dy_acmono_w7 are still being sampled with distance-constrained CREST
(`conf_v3/`), because unconstrained GFN2 sampling ejects a first-shell water.

## Result 3 — acetate denticity (mono and bi are one species)

| step (conformer-best, kcal/mol) | La | Dy |
|---|---|---|
| mono_w7 → bi_w7 (isomers) | −1.80 | −0.50 |
| mono_w8 → bi_w7 + H₂O | +0.63 | −1.64 |

At a fixed water count, bidentate wins for both metals. When the second O has to displace a water, La
keeps the water and Dy chelates. That makes the lowest acetate state metal-dependent, so `analyze.py` now
pools mono and bi as ligand `ac`. Merged result: **ΔG = −1.95** (ΔH −2.77; La_acmono_w8, Dy_acbi_w7,
La_aq9, Dy_aq8). La bidentate with 8 waters (CN 10, the CREST-preferred mode) was never computed; it is
running as `runs/La_acbi_w8` and `conf_cn10/`.

## Result 4 — EDTA from the crystal structure (`edtax`)

The xtb-wrapped `edta` starts unwrapped under DFT and are excluded as artifacts. The `edtax` starts come
from COD 7236231, Na[Dy(EDTA)(H₂O)₃]·5H₂O, with La built by metal swap; the CN 8 variant drops the
2.535 Å water. All four stay hexadentate after Opt:

| state | M–N | M–O(carboxylate) range | imag. |
|---|---|---|---|
| La w3 | 2.75, 2.75 | 2.45–2.66 | none |
| La w2 | 2.72, 2.75 | 2.42–2.63 | −79 cm⁻¹ |
| Dy w3 | 2.62, 2.64 | 2.32–2.58 | none |
| Dy w2 | 2.56, 2.61 | 2.30–2.46 | −66 cm⁻¹ |

Best states give ΔG = **−6.35** (La w3, Dy w2, La aq9, Dy aq8). Using the minimum-only Dy w3 instead
gives −5.89. Either way, against the benchmark of about −3.9, Dy is over-favoured by **2.0–2.5 kcal/mol**.
The two w2 states still carry an imaginary mode and have not been conformer-searched.

**Cross-ligand pattern.** Acetate over-favours Dy by about 2 (−1.95 against about 0) and EDTA by about
2–2.5. For a monodentate ligand and a hexadentate one, the error is roughly the same size. That suggests a
ligand-independent error in the common aqua-ion leg, La(H₂O)₉ against Dy(H₂O)₈, rather than an error
in metal–ligand bonding. This is falsifiable: if the error is in the aqua leg, LBT3 should also come out
about 2 too Dy-favoured, and correcting only the aqua leg should fix all three systems. If the EDTA
error grows with denticity once the w2 states are cleaned up, the hypothesis is wrong.

## Ingredient tests (partial)

Fixed-geometry single points: `diagnostics/lady_ingredient_tests_20261026/`.
- **The CPCM contribution to the exchange is +11 to +17 kcal/mol** (the gas-phase exchange is −10 to −28).
  The answer is a small residual of two large terms, so a 10–15% continuum error equals the entire
  discrepancy.
- ORCA used its fallback cavity radius of **2.40 Å for both La and Dy**. The radius sensitivity runs
  (`r_ionic` with Dy 2.14 Å, `r_small` 2.0, `r_large` 2.8; `make_radius.py`) test this directly.
- Bare-ion ΔΔG_hyd(Dy−La) at large core is −59, against about −67 experimental (recalled Marcus,
  unverified). The Dy aqua ion appears under-stabilised relative to La, which is in the direction of the
  cross-ligand offset.
- The small-core 4f and DLPNO-CCSD(T) rows are still running.

## Second hydration shell (running)

`build_ss.py` takes 12 states × 3 seeds = 36 clusters. Each core is the conformer-best
first-shell structure, plus 12 waters at H-bond acceptor sites (61–68 atoms). In the GFN2
pre-opt the core is frozen. An unconstrained pre-opt (build 24894) pulled waters into the first
shell (CN 8→9–10), so it was discarded. Where the frozen-core pre-opt still changed CN, the raw
placement is used. This tests whether the missing explicit second-shell / dielectric
differential is what over-favours Dy. Outputs in `run_ss/`.

## Protocol issues found (for the next person)

- TightOpt on first-shell clusters stalls for 200–300 cycles on the 1e-6 Eh energy criterion,
  through floppy H-bond modes, once the gradients have already converged. `restart_job.sh` resumes
  with normal Opt + Freq.
- CREST 3.0.2 native mode crashes in GFN0 WBO setup (`DSYGVD` illegal value) on Ln clusters, so
  use `--legacy`. Legacy mode then aborts in `trialopt` because M–O sits at the bond cutoff, so
  also use `--noreftopo` and filter conformers by CN (`conf_dft_array_v2.sh`).
- GFN2-xTB over-coordinates Ln³⁺: unconstrained relaxation pulls in extra waters. Constrain or
  freeze M–O whenever xtb touches these clusters.

## Files

`build.py` starts · `make_inputs.py` / `orca_job.sh` primary runs · `restart_job.sh` ·
`crest_job*.sh`, `conf_dft_array*.sh` conformers · `build_ss.py`, `build_ss_job.sh` second shell ·
`analyze.py [--conformers]`. Outputs are in `workspaces/lady_relaxed_exchange_20261024/`
(`EXCHANGE_RESULT*.json`, `JOBS*.txt`).
