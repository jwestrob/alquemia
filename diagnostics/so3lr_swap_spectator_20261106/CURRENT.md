# Test 1 complete — 6 November 2026

L passes the frozen0.5-kcal/mol energy criterion: maximum distant-water effect
0.0670 atLR1000/0.0677 no cutoff. M sensitivity misses it:0.5821/0.5827.
All62 model calls/checks completed; all numerical tables and decisions reproduce.
Read REPORT.md and RESULTS.json. This is a charged isolated acetate/aquo
prototype, not neutral periodic-box or chemical-accuracy validation.

Force differences still change (largest L component about4.18 kcal/mol/A).
Do not equate energy cancellation at one geometry with force/ensemble validation.
Next recommendation: prepare the intended neutral ligand system and repeat
energy/force checks before trajectories. No Test2, MD or new DFT was started.
The original screen's numerical data/reference limitations remain unchanged.

Jobs52169 and52184 completed. Total640 allocatedCPU-s (10.7CPUmin) within3600;
175.149 measuredCPU-s, zeroGPU/DFT/MD. No live job or duplicate/restart needed.
Raw arrays and executed source copies remain in the new ignored workspace.
reproduce.sh queues saved-only analysis to a fresh directory. Preserve M1,
reserved outcomes, historical files and other-session executors.
