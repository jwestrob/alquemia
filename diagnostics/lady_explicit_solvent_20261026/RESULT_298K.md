# Explicit-solvent La->Dy exchange, 12-6-4 / OPC / GAFF2-BCC, 298 K (2026-10-27)

Alchemical La->Dy (linear mixing of the metal 12-6-4 pair terms), 13 lambda x 6 ns, seed 1,
neighbour-BAR (pymbar MBAR failed its normalization check), analytic r^-4/r^-6 tail added.
Data: workspaces/lady_explicit_solvent_20261026/{aq,ac,edta}/fep/T298/s1; FEP_RESULT.json.

| kcal/mol            | explicit 12-6-4 | static PBE0/CPCM | experiment |
|---------------------|-----------------|------------------|------------|
| dG_alch(aq) Dy-La   | -70.2           | -59.2            | -66.9 (12-6-4 fitted to HFE: not a test) |
| acetate exchange    | -2.29           | -1.95            | -0.07 (I=0.1) |
| EDTA exchange       | -10.49          | -6.35            | -4.01 (NIST 46) |

BAR errors are 0.03-0.05 but are not the real uncertainty: EDTA donor count hops 4 <-> 6 across
intermediate windows (slow binding-mode sampling).

End-state structure: La aq CN 9.00, Dy aq CN 8.71 (8: 29%). Acetate La 95% monodentate, Dy 89%
bidentate. EDTA: La tetradentate 71% / pentadentate 29% (crystal is hexadentate), Dy hexadentate 81%.

Conclusion: explicit solvent does not rescue the static model; it over-selects Dy further, and the
excess grows with the number of carboxylates (ac -2.2, EDTA -6.5; static -1.9, -2.3). The
"missing solvent statistics" hypothesis predicted a move toward experiment with the selectivity
shifting into TdS; dG moved the other way. Not decisive (the metal-carboxylate terms are transferred
from water-fitted parameters, and La-EDTA unwraps), but both very different models share the same
sign of error per carboxylate, which points at the metal-carboxylate interaction itself.

Multi-temperature (dH/TdS) windows and seed 2 were cancelled: with dG off by 6.5 and the La
binding mode wrong, the van 't Hoff split would describe this force field, not the experiment.
