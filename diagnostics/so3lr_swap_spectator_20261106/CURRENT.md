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

7 November continuation: the active goal is to proceed until there is a good
reason to stop. The neutral cis/trans system is prepared and its first controls
are complete; H200 job 52844 is pending with a completion watcher. Read
../so3lr_neutral_pair_20261107/CURRENT.md. This subsequent authorization does
not restart the old quantum references or change the completed Test 1 results.

9 November closure: the subsequent neutral-system evaluation is complete. Its
numerical controls passed, but the intended trans ligand failed covalent
integrity during native preparation and a separate capped ASE comparison.
Read ../so3lr_neutral_pair_20261107/REPORT.md and DECISION.json. Stop that setup;
no longer MD, alchemy or reference retry is released. All jobs in that phase
are terminal. This does not change the completed Test 1 results.
