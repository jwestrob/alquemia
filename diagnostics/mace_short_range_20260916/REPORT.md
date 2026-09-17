# Whole-protein short-range learned component

This is a separate empirical descriptor, not a complete energy or binding free energy.
It uses the original trained `interaction_energy` readout, before charge restoration and field updates.
No new MACE, GB or DFT evaluation. Baseline and failed total-energy descriptors unchanged.

| Case | Full medium Rshort | Full large Rshort | Local corrected-H Rshort | Full−local corrected-H |
|---|---:|---:|---:|---:|
| GGR_1GLG | -17.570592 | -14.992743 | -9.090123 | -8.480469 |
| ALPHA_1F6S | -22.689014 | -21.823952 | -17.165250 | -5.523764 |
| ALPHA_6IP9 | -28.373353 | -26.949454 | -23.313390 | -5.059962 |
| PQQ_1H4I | -9.153793 | -6.267153 | -3.790932 | -5.362860 |
| PQQ_4MAE | 14.442741 | 21.219776 | 21.450287 | -7.007546 |

All values kcal/mol; larger was predeclared as more La-like. No inherited zero/band.

## Predeclared comparisons

```json
{
  "medium": [
    {
      "higher_expected": "PQQ_4MAE",
      "lower_expected": "PQQ_1H4I",
      "difference_kcal_mol": 23.596533612614323,
      "expected_order": true
    },
    {
      "higher_expected": "ALPHA_1F6S",
      "lower_expected": "GGR_1GLG",
      "difference_kcal_mol": -5.118421436893101,
      "expected_order": false
    },
    {
      "higher_expected": "ALPHA_6IP9",
      "lower_expected": "GGR_1GLG",
      "difference_kcal_mol": -10.802760683115597,
      "expected_order": false
    }
  ],
  "large": [
    {
      "higher_expected": "PQQ_4MAE",
      "lower_expected": "PQQ_1H4I",
      "difference_kcal_mol": 27.486929798259155,
      "expected_order": true
    },
    {
      "higher_expected": "ALPHA_1F6S",
      "lower_expected": "GGR_1GLG",
      "difference_kcal_mol": -6.831209433130182,
      "expected_order": false
    },
    {
      "higher_expected": "ALPHA_6IP9",
      "lower_expected": "GGR_1GLG",
      "difference_kcal_mol": -11.95671053889261,
      "expected_order": false
    }
  ],
  "local_archived": [
    {
      "higher_expected": "PQQ_4MAE",
      "lower_expected": "PQQ_1H4I",
      "difference_kcal_mol": 25.11112641752549,
      "expected_order": true
    },
    {
      "higher_expected": "ALPHA_1F6S",
      "lower_expected": "GGR_1GLG",
      "difference_kcal_mol": -11.548778021839649,
      "expected_order": false
    },
    {
      "higher_expected": "ALPHA_6IP9",
      "lower_expected": "GGR_1GLG",
      "difference_kcal_mol": -14.40267696467927,
      "expected_order": false
    }
  ],
  "local_global_H": [
    {
      "higher_expected": "PQQ_4MAE",
      "lower_expected": "PQQ_1H4I",
      "difference_kcal_mol": 25.241219184440887,
      "expected_order": true
    },
    {
      "higher_expected": "ALPHA_1F6S",
      "lower_expected": "GGR_1GLG",
      "difference_kcal_mol": -8.075126535817585,
      "expected_order": false
    },
    {
      "higher_expected": "ALPHA_6IP9",
      "lower_expected": "GGR_1GLG",
      "difference_kcal_mol": -14.223267216419757,
      "expected_order": false
    }
  ]
}
```

Primary development screen passes: False.
No threshold/weight was fitted. Cases are consumed development; alpha is one observation.
The component has finite-neighborhood sensitivity and no explicit total-charge/spin response.
Do not report total-energy forces as its gradients. Independent evaluation remains necessary.
