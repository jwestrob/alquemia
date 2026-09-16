# Stage C: executed analytic-gradient integration

Final review of completed job **1199805**, using the frozen
`workspaces/ggr_mechanism_20260915/stage_c_tasks_v1/manifest.json`
(SHA256 `248b0e30b4a9154e030dac4e7d80eaa9a0c49af131971a752933fd8e7b7db89b`).
All **20 endpoints** are complete: four analytic-gradient centers and sixteen
displaced single points, covering the 58-atom extended-amide and 111-atom
connected models. All **four coordinate blocks / twelve La, Ca and R checks
pass**. The half-step branch was not triggered.

Final collection: `stage_c_tasks_v1/collection_1199805.json`, SHA256
`6706fc8ad34850b1960e8db12def1d6c4a31638ed6dacebf2735aad14e5c5e60`
(workspace prefix as above).

## What actually passed

- All four ORCA 6.1.1 centers terminated normally and produced `.engrad` files.
  The adapter reads the actual scalar gradient/four-column coordinate format.
  Total energies match the output; all true atomic numbers and atom coordinates
  match the pinned XYZ within the declared tolerance. La is atomic number 57,
  with 46 ECP core electrons removed; explicit counts are La 224 / Ca 234 for
  the extended model and La 440 / Ca 450 for the connected model.
- Output reports analytic SCF, CPCM, dispersion, gCP and Cartesian gradients;
  La additionally reports its ECP gradient completed. No numerical-gradient
  fallback appears. Collection now requires these actual component markers.
- Full input/XYZ/output/gradient/runtime-artifact and execution-receipt pins
  pass, including the approved ORCA executable, charge 0/−1 and singlets.
- Physical source mapping and paired `grad(E_Ca−E_La)` extraction work on the
  actual gradient files. Neither synthetic cap coordinates nor a response
  energy are introduced.
- TightSCF-minus-normal center changes all pass the frozen 0.05 kcal/mol bridge
  criterion. This checks numerical compatibility, not an affinity label:

| Representation | La change | Ca change | R change | Bridge |
|---|---:|---:|---:|---|
| Extended | −0.0001270549141 | −0.0000078781711 | +0.0001191767429 | pass |
| Connected | +0.0003799859272 | +0.0004974613641 | +0.0001174754369 | pass |

**Tests:** `python -m unittest discover -s tests -p test_ggr_sensitivity.py -v`
ran **15 tests, all passing, no skips** after the real gradients arrived.
Eleven are software/geometry/archived-energy checks; four use the executed
gradient artifacts, including actual format variants and explicitly corrupted
copies for missing-component, charge and malformed-file rejection. Two earlier
targeted existing `affordable_response` regressions also passed.
After both connected centers arrived, the four relevant real-gradient
integration tests were rerun: **4 passed in 8.822 s**, validating all four
centers without rerunning quantum calculations.
Final integrity review additionally verified all 20 raw endpoint energies and
execution/runtime/gradient pins, then independently recomputed all twelve
odd-energy, projected-gradient and even-energy checks directly from the actual
`.engrad` arrays, frozen QM Jacobians and displaced output energies. Maximum
absolute residual is **0.0005865567598238412 kcal/mol**, connected-model La
metal motion, against the frozen 0.02 kcal/mol tolerance. No extra DFT ran.

## Final directional checks

Both named coordinates pass for each endpoint and for R in both models.
Residual means actual odd energy change minus the center-gradient prediction,
in kcal/mol. All twelve
checks use the frozen 0.02 kcal/mol absolute tolerance floor.

| Model / coordinate | La residual | Ca residual | R residual | Paired projected derivative |
|---|---:|---:|---:|---:|
| Extended metal ±0.02 Å | +0.0004242650 | −0.0000161584 | −0.0004404234 | +7.7812941220 kcal/mol/Å |
| Extended peptide ±1° | +0.0004499546 | +0.0000819113 | −0.0003680433 | +6.6915182697 kcal/mol/radian |
| Connected metal ±0.02 Å | +0.0005865568 | +0.0000293730 | −0.0005571838 | +3.4735357989 kcal/mol/Å |
| Connected peptide ±1° | +0.0004083764 | +0.0001936247 | −0.0002147517 | +5.9382470570 kcal/mol/radian |

The actual R odd changes in table order are +0.1551854590, +0.1164209825,
+0.0689135322 and +0.1034272113 kcal/mol. No block is eligible for the
conditional half-step branch (`half_step_blocks=[]`, all four
`half_step_eligible=false`). No half-step evaluations were needed.

Even terms are `([E(+h)−E(0)]+[E(−h)−E(0)])/2`, in kcal/mol:

| Model / coordinate | La even term | Ca even term | R even term |
|---|---:|---:|---:|
| Extended metal | +0.1129035242 | +0.0498186886 | −0.0630848356 |
| Extended peptide | −0.0107220286 | +0.0060928848 | +0.0168149133 |
| Connected metal | +0.1159409165 | +0.0508639321 | −0.0650769844 |
| Connected peptide | −0.0072994178 | +0.0092563043 | +0.0165557221 |

These terms include ordinary quadratic curvature; they do not by themselves
demonstrate anharmonicity or a stable mechanical basin. Negative values remain
reported without clamping. No stiffness, relaxation energy or entropy is fitted
from these checks.

## Traceable center costs and artifacts

Endpoint directories below are relative to `stage_c_tasks_v1/`. Each contains
`endpoint.out`, `endpoint.engrad` and `endpoint.out.execution.json`.

| Endpoint | Receipt wall-seconds | Assigned rank-seconds | `.engrad` SHA256 |
|---|---:|---:|---|
| `extended_center_La` | 256.308216 | 4100.931456 | `1d16fadef6e1c31fa72d9a2126c3bea73f92ff860e92bc30f7d92cd6cf1ee844` |
| `extended_center_Ca` | 244.401423 | 3910.422768 | `dac641df38ddbf772baca21955e5d7935ef350f8e55d368b80af5b84fe903301` |
| `connected_center_La` | 790.883137 | 12654.130192 | `985972f9033fa6df9967720229da1032f6eb5b7d89effd832bd4c29b0a2937d5` |
| `connected_center_Ca` | 762.249228 | 12195.987648 | `dc9913d15b6f06b52f504fd27f6a6a190ab5da6a487a4cdbc78545e7331d5df5` |

All used 16 MPI ranks, one thread per rank. The extended pair sums to
500.709639 receipt wall-seconds / 8011.354224 assigned rank-seconds; the connected
pair sums to 1553.132365 / 24850.117840. All four centers sum to
2053.842004 / 32861.472064. These are endpoint sums, not elapsed batch time or
the full Stage C cost. ORCA's gradient-module times are 49.295 s (extended La),
36.789 s (extended Ca), 179.919 s (connected La), 151.775 s (connected Ca).

For all 20 Stage C endpoints, receipt sums are **8156.495456 wall-seconds**
and **130503.927296 assigned rank-seconds**. These include the actual centers
and displaced points, not a proposed routine per-site production cost.

## Final interpretation and reproducibility

These checks validate analytic-gradient integration and consistency with the
two prescribed physical displacements in both local models. The paired metal
sensitivity changes materially between representations (7.78 to 3.47
kcal/mol/Å); the peptide sensitivity changes less (6.69 to 5.94
kcal/mol/radian). This is model sensitivity on already inspected GGR geometry,
not broad biological validation, whole-protein convergence or an affinity
correction. Only the two tested directions and amplitudes are validated.
The two approved motions leave link-cap anchors fixed,
so this is not an end-to-end moving-cap force validation.

No acceptance threshold, input, method, amplitude or task count changed during
integration. Technical additions only strengthened actual-output validation and
tests. Collections record the current collector/adapter hashes separately from
the immutable preparation snapshot. No additional scientific runs were launched
by this check. Relaxation and entropy corrections remain null with
`response_model_not_validated`.

Recollect the completed nominal tasks using the pinned Stage A reference manifest:

```bash
python scripts/ggr_sensitivity.py collect \
  --manifest workspaces/ggr_mechanism_20260915/stage_c_tasks_v1/manifest.json \
  --output workspaces/ggr_mechanism_20260915/stage_c_tasks_v1/collection_recheck_1199805.json
```

Use a fresh output filename if that record already exists. This reads existing
artifacts and performs no quantum calculations. All nominal checks passed, so
the already approved conditional half-step branch remains unused.
