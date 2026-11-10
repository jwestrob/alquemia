# Deviations from PLAN.md — NIST charge response, 9 November 2026

The frozen design (PLAN.md sha256 `c2c04a1f…9246`, frozen 01:43 UTC 10 November) ran as written.
Everything below lies outside it.

1. **Charge-only baseline Mq.** Requested by the lead agent after the frozen fit had run. The subagent
   fixed Mq and the neutral-donor definition in the docstring of `add_charge_baseline.py` at 01:48 UTC
   before running it.
2. **Ionic-strength cap.** The PLAN's reason for the 0.1 M cap understated the error of mismatched
   strengths: moving Ca alone from 0 to 0.1 M shifts its log K by about 0.43 per unit of ligand
   charge. Six of the seven admitted nearest-strength pairs put La at 0.1 M and Ca at 0 M and sit
   among the ten smallest differences. The same-strength set gives the same reading.
3. **All-pairs charge line.** Written and run by the lead agent after the subagent's results were
   known (`charge_line_all_pairs.py`, sha256 `48c7a3b4…a87bc`). Its design is fixed in its docstring.
   Charge comes from the NIST protonation label; donor classes come from name keywords. An earlier
   exploratory pass by the lead agent over the same table (charge against the difference, no donor
   classes) preceded the script and motivated it.
4. **Protein comparison.** The line's coefficients were on disk before the protein table was read. The
   choice of descriptor (acidic residues in the first shell) and of summary statistics came after the
   lead agent had seen the protein rows. It is descriptive.
5. **Report authorship.** `REPORT.md` and this file were written by the lead agent, because the
   harness returns subagent findings as text.
