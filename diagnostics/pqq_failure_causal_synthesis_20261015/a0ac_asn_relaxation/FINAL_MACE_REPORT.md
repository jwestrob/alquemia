# Completed flipped-basin MACE result: two of ten sources change

Jobs 4831 and 4844 completed all **20 declared searches and 40 final paired endpoint cells**, retaining the full **ten-source denominator** and 80 exact reused baseline cells. All optimizer receipts report success and physical feasibility. The final collection incorporates the actual scout; it does not count that source twice.

Common3 contains origin/adaptive_Ca/adaptive_La; common4 adds the fixed amide flip; common6 adds relaxed_flip_Ca/relaxed_flip_La. Mathematical and operational selections agree for every source. Original-first ties and the strict −0.1 kcal/mol operational rule remain unchanged.

## Actual common4→common6 response

Values below are saved MACE vacuum changes in kcal/mol. The paired sign is ΔE(Ca)−ΔE(La), so a positive change reflects a larger lowering of the selected La energy.

| Source | Common4 winners (Ca / La) | Common6 winners (Ca / La) | Ca change | La change | Paired change |
|---|---|---|---:|---:|---:|
| Ca-conditioned sample 4, target | adaptive_Ca / amide_flip | relaxed_flip_Ca / relaxed_flip_La | −0.2860356038 | −16.8133988373 | **+16.5273632335** |
| La-conditioned sample 4, control | adaptive_Ca / adaptive_La | adaptive_Ca / relaxed_flip_La | 0 | −11.4637900423 | **+11.4637900423** |
| Ca samples 0–3 and La samples 0–3 | adaptive_Ca / adaptive_La | unchanged | 0 | 0 | **0** |

Thus the relaxed alternative changes both sample 4 sources, not uniquely the target. Both exceed the predeclared sign-blind 1 kcal/mol solvent-review gate. This report makes **no composite-energy or class call**; no class bands were applied.

## Qualification and limits

The original four modes/settings, cached flipped q0 energies/forces, rebased projection, covalent checks and full-source relative-contact admission were retained. Four proposals carry boundary flags: both metals for Ca sample 2, La for La sample 0, and Ca for La sample 2. None wins its common6 pool; the new selected sample 4 proposals have no boundary flag.

These remain constrained, local vacuum searches around an admitted flip, not unconstrained or equilibrium minima. Severe-contact admission does not establish clash-free geometry; named-contact interpretation belongs to the separate contact report. Exterior-source atoms affect geometry admission but contribute no exterior energy. No further modes, protonation changes or retries were introduced.

## Measured new-stage cost

| Stage | Actual model forwards | Receipt elapsed seconds | Allocated CPU-seconds |
|---|---:|---:|---:|
| Scout 4831 | 30 | 569.195202 | 3984.366412 |
| Remaining 4844 | 236 | 4196.215674 | 29373.509715 |
| Total | **266** | **4765.410875** | **33357.876127** |

Both workers used seven CPU threads, with zero failed requests; peak worker RSS was 7459376/7569404 KiB respectively. There were zero new SCC, DFT or GPU calls. Totals combine the two application receipts, excluding historical baselines, preparation, queuing and the separate solvent stage; they are not full-campaign latency or scheduler TotalCPU.

## Immutable evidence

Under `workspaces/pqq_failure_causal_synthesis_20261015/`:

| Artifact | SHA256 |
|---|---|
| `asn_relaxation_remaining_v1/execution/FINAL_COLLECTION.json` | `2ed05f5a241f4f1eccbb2a680428a4b032a51449df56ad323e138b83cf32fd82` |
| Remaining `EXECUTION.json` | `0834eba9aa82e8cc33d137dfb973dedcbe0513b88ecdd2c1801b2f0b6f1ad8e4` |
| `asn_relaxation_scout_v1/execution/FINAL_COLLECTION.json` | `e0d0e8e31eddf6954dfcf77ade9a871e871bdf67c8e688fb9f0f67fcc0ce758f` |
| Scout `EXECUTION.json` | `091ea56153952ddd0b1c978594fa2b1235085c94e0c06ac4b8ea2d621f994ce1` |

This report reads existing results and resource receipts only; it performs no molecular or scientific reanalysis.
