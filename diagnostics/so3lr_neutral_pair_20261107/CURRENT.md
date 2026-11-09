# Closed bounded evaluation — 9 November 2026

Read REPORT.md, FINAL_RESULTS.json and DECISION.json. The current Tauri setup is not qualified for intact-ligand La/Dy selectivity. End automatic expansion. All scientific jobs in this diagnostic are terminal; no longer MD or reference retry is released. Preserve M1 and every failed input/output.

The decisive finding is chemical integrity. Trans600 C6–C7/C7–O8/C7–O9 separate during native FIRE, first observed at saved step 20 (step 10 was intact). The longest bond is 13.076 Å before MD and 12.617 Å after 0.25 ps. Atom order, periodic-image mapping and checkpoint/HDF correspondence pass. A separate capped ASE LBFGS diagnostic crosses the original 2 Å range for C6–C7 at step 37 (2.052 Å), with all atom steps within 0.05 Å. It stopped at the declared threshold; do not describe this capped result as completed scission or assign a unique physical cause.

Bridge 66047 failed its fixed 1e-5 eV/Å tolerance (actual 1.964e-5); no optimizer ran there. ASE-only comparison 66081 preserved that failed record without qualifying the native interface. No further optimizer search or MD followed.

Other controls pass: all three original periodic tests (max 0.01560 kcal/mol), the hydrated test (0.02439 versus 0.5), and CPU/GPU agreement. Native MD completed 1,000 steps but began with damaged chemistry. Thermal stability is not affinity evidence. The other five starts and all alchemical calculations were not run.

GPU 52844's 62-second failure was a disabled CPU callback backend. Runtime-only recovery 65869 enabled cuda,cpu, reused two predictions and completed. Hydrated 65896, scout 65917, analysis 65972, integrity audit 66010 and capped comparison 66081 are terminal. Preserve failed audit 66003 and its source-reconstruction note. Exact costs are in ACCOUNTING.json; the GPU phase used 894 seconds including failures.

The experimental common medium and small-core Dy reference remain unqualified. Do not compare endpoint gaps or this damaged trajectory with the conditional experimental contrast. Reserved outcomes, other sessions' jobs, historical calculations and the production baseline remain untouched.

Receipts for 52844/65869/65896/65917 have been consumed and ACKed. Delayed copies are ACK-only: never resubmit from those events. The reproduction command repeats saved analyses only.

Final saved-data replay 66238 completed and exactly reproduced the 254 trajectory ligand records, coordinate/bond audit and final decision. All 38 capped-optimizer frames were independently verified. See FINAL_REPLAY_CHECKS.json and ARTIFACTS.json. Final accounting covers 27 terminal jobs: 23,810 allocated CPU-s, 6,803.388 measured CPU-s and 894 H200 GPU-s. CPU-only allocations were 13,082 CPU-s, within the 14,400-s envelope. TOTAL_TAURI_COST.json includes the earlier phases separately.
