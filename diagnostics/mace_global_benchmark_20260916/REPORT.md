# Whole-protein MACE / frozen-solvent development panel

Production baseline unchanged. Five consumed structures, four biological groups.
Alpha structures are one qualified affinity observation; PQQ labels describe functional association.
R = E(Ca) - E(La). Solvent R adds G(Ca) - G(La). All values below are kcal/mol.
No aquo offset, inherited decision band, fitted threshold or absolute class.

| Case | Medium vacuum R | Medium solvent R | Large solvent R | Large−medium solvent R |
|---|---:|---:|---:|---:|
| GGR_1GLG | -405441.555705 | -405476.714009 | -405435.162831 | 41.551179 |
| ALPHA_1F6S | -405365.440185 | -405497.973561 | -405492.318000 | 5.655561 |
| ALPHA_6IP9 | -405376.030182 | -405511.726741 | -405507.889519 | 3.837221 |
| PQQ_1H4I | -405328.594043 | -405432.858906 | -405291.030995 | 141.827911 |
| PQQ_4MAE | -405363.215636 | -405450.487074 | -405353.480553 | 97.006521 |

## Predeclared relative-order checks

```json
{
  "medium": [
    {
      "difference_kcal_mol": -17.628167761780787,
      "expected_order": false,
      "higher_expected": "PQQ_4MAE",
      "lower_expected": "PQQ_1H4I"
    },
    {
      "difference_kcal_mol": -21.25955179228913,
      "expected_order": false,
      "higher_expected": "ALPHA_1F6S",
      "lower_expected": "GGR_1GLG"
    },
    {
      "difference_kcal_mol": -35.01273113646312,
      "expected_order": false,
      "higher_expected": "ALPHA_6IP9",
      "lower_expected": "GGR_1GLG"
    }
  ],
  "large": [
    {
      "difference_kcal_mol": -62.449558199499734,
      "expected_order": false,
      "higher_expected": "PQQ_4MAE",
      "lower_expected": "PQQ_1H4I"
    },
    {
      "difference_kcal_mol": -57.15516990568722,
      "expected_order": false,
      "higher_expected": "ALPHA_1F6S",
      "lower_expected": "GGR_1GLG"
    },
    {
      "difference_kcal_mol": -72.72668889880879,
      "expected_order": false,
      "higher_expected": "ALPHA_6IP9",
      "lower_expected": "GGR_1GLG"
    }
  ]
}
```

Medium primary ordering screen passes: False.
Numerical checks: {'medium': True, 'large': True}.

The ordering screen is a development result, not broad validation or goal completion.
Frozen monopoles omit dipolar/self-consistent solvent response; metal cavity parameters remain unvalidated.
No combined descriptor gradient, relaxation or entropy correction is available.
Exact inputs, exclusions, terminal/H operations and chemical states are pinned per case.
Per-task receipts retain evaluation/initialization time and memory; whole Slurm costs are reported separately.

## Outcome and measured cost

Both medium and large solvent descriptors reverse all three predeclared
relative-order checks. Medium differences are −17.628168 kcal/mol for
4MAE−1H4I and −21.259552/−35.012731 for the two alpha−GGR comparisons.
Large gives −62.449558 and −57.155170/−72.726689 kcal/mol. Do not invert the
sign convention or fit a threshold to turn this development failure into a win.
Vacuum descriptors retained alpha above GGR but also reversed the PQQ pair.

All 24 MACE and 24 GB calls completed, with no failed attempts or new DFT.
Jobs 1200701/1200702/1200711/1200712 total 1,370 allocated GPU seconds
(22.83 GPU-minutes), 21,920 allocated core-seconds, 1,475.483 actual CPU seconds.
Medium primary pair evaluations: alpha17–18s, GGR49.0s, PQQ113–117s.
Large: alpha41–42s, GGR110.1s, PQQ234–242s. Startup/validation overhead is
included in allocation totals. Peak allocated GPU10,208,399,872 bytes medium,
16,012,973,056 large; peak reserved12,041,846,784/22,099,787,776 bytes.
OpenMM GPU memory was not measured. Full details are in costs.json.

Numerical credibility: admitted. Tested medium rigid energy errors <=1.45e-5
kcal/mol, forces <=3.25e-7eV/A; GB follows transformed learned monopoles and
also passes. Large reuses its prior real full-system rotation validation; no
new large rotation call was claimed for this panel. A fresh preparation replay
reproduces all five physical arrays and all ten endpoint XYZ hashes exactly.

Scientific usefulness: this direct whole-protein/frozen-solvent construction
fails the declared development screen. The earlier solvent-induced checkpoint
agreement did not predict biological usefulness. Affordability: minutes per
protein pair on A5000, with a small GB increment; compute cost is no longer the
main obstacle. This is a failed model version, not a failure of the baseline.

Retain baseline. Continue the active MACE goal by separating local descriptor
errors from the whole-protein contribution using actual matched cores and
archived DFT; do not launch a large calibration of this failed direct version.
Combined response/relaxation remains unavailable. No production promotion.
