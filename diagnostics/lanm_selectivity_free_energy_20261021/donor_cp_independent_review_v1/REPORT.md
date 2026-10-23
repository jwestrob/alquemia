# Independent donor counterpoise review

**The balanced formulas, signs, state inventory and units are correct.** The original full collector has one material admission gap: it verifies output/input hashes and stored success, but does not replay raw energies or bind each result's full cell and executed input to its declared cell. `REVIEW.py` closes that gap separately; the frozen runner, panel, scientific inputs and thresholds remain unchanged.

The strict reviewer requires the exact frozen panel, unique complete row inventory, matching result/cell/source pins, actual `START.json`, job/rank identity, the actual memory-to-MaxCore rule, and byte-exact executed input regenerated for that cell. It replays native admission and compares every energy used in the cycle to the raw output. Saved real-complex PBE0 electronic and D4 components are independently re-extracted. Missing cells never become zero.

## Actual saved evidence

The existing ten DFT rows pass all **294 checks**. Independent 40-digit decimal arithmetic reproduces the saved DFT subset within **2.05×10−11 kcal/mol**. No native programs were invoked; review took 2.84 seconds on one local thread. All twelve CC slots remain explicitly unavailable to this subset review; this does not report their current scheduler state or manufacture a full comparison.

| Donor | PBE0 electronic CP paired exchange | PBE0-D4 CP paired exchange |
| --- | ---: | ---: |
| Acetamide | +4.33338133 | +4.66935492 |
| Bidentate acetate | −6.10343645 | −5.55346267 |

Values are kcal/mol and reproduce already consumed data. The canonical-method comparison awaits the complete 22-cell denominator.

## Interpretation retained before CC outcomes

For each metal, subtract the full-basis ghost-ion and ghost-ligand energies from the real complex, then subtract that metal's water interaction. Take **Dy minus La** last. The equivalent independent arrangement is the donor Dy−La interaction gap minus the water Dy−La interaction gap. Thus, if δL and δW are the corresponding PBE0-minus-CC method discrepancies, the balanced exchange error is **δL−δW**. This tests cancellation rather than comparing unlike metal total energies.

D4 uses only real molecules: complex minus real free ligand, with isolated-ion D4 zero. Free-ligand D4 cancels in each paired metal contrast but remains in the individual interactions. CC receives no D4. The acetate charge balance is (+3)+(−1)→(+2)+0; the amide is (+3)+0→(+3)+0. Fragment bases, frozen-f ECPs and effective singlets remain matched. The 0.25 kcal/mol resolution is unchanged and is applied separately to each donor. These fixed electronic cycles do not establish aqueous affinity, physical spin/SOC accuracy, basis convergence, or force accuracy.

## Full collection use

Pass exactly the 22 actual result paths after root admits their completion:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python diagnostics/lanm_selectivity_free_energy_20261021/donor_cp_independent_review_v1/REVIEW.py \
  --panel workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/donor_cp_v1/prepared_v1/PANEL.json \
  --out WORKSPACE_FOR_NEW_FULL_REVIEW --mode full RESULT_PATH_1 ... RESULT_PATH_22
```

Require `full22_strict_review_passed` before persisting the original collector's result. Full mode also evaluates that frozen collector in memory and compares its balanced cycles and decision predicate to independent decimal arithmetic. It saves signed water-only and final exchange method gaps without selecting a favorable method after inspection. A failed review stops collection; no scientific recovery is implicit.

Source SHA256: `eef93c84d9a469dd3a15fbe2ed2213140b0e4a51f2c6ddacbbfc4299c76cd379`. Actual DFT subset review: `workspaces/lanm_selectivity_free_energy_20261021/donor_cp_independent_review_v1/dft_subset_v1/RESULT.json`, SHA256 `5b8f0047564757cb78419b8760ccf6fc9eeaddfcc02e9a92bdbc7032329d4327`.
