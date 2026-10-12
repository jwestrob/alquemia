# Independent context × accommodation review

12 October 2026. **The integrated arithmetic is correct. Both context and the
candidate pool affect the result; enlarging the context does not remove the
large differential response of the old q9 La proposal.** The common five-candidate
pool retains positive two-anchor ordering only in the expanded context, with a
smaller margin than the expanded-context-only proposals suggest.

## Recomputed ordering

The table gives `R(c5axv8) − R(q9z4j7)`, where `R = E_Ca − E_La` and
`E = MACE_vac + GFN2_water − GFN2_vac`, in kcal/mol. Positive values preserve the
expected ordering of these two consumed La/Ca anchors. These are **not new
25-reference calibration gaps, classifications, affinities or validation cases**.

| Shared candidate pool | Original context | Expanded context |
|---|---:|---:|
| Origin only | +8.867782 | +12.175484 |
| Origin + original-context Ca/La proposals | −4.556602 | −1.212842 |
| Origin + expanded-context Ca/La proposals | +0.437421 | +3.472171 |
| Origin + all four saved proposals | −3.714195 | +1.265264 |

Both metals receive the same pool within each context. Selection retains the
origin unless the best candidate improves its own metal energy by more than
0.1 kcal/mol, with the recorded deterministic candidate order. Mathematical and
operational selections coincide for these completed pools. No cross-context
absolute energy minimization is used.

## What this establishes

The old three-candidate pool still reverses these anchors when evaluated in the
expanded context. Therefore the earlier favorable expanded result cannot be
attributed solely to pocket enlargement: its different donor proposals matter.
Conversely, with a fixed common five-candidate pool, changing context improves
the two-anchor ordering by **4.979459 kcal/mol**, so proposal selection alone
does not explain the difference either.

For q9, La selects the **original La proposal** from the common pool in both
contexts; Ca selects the **expanded La proposal**. In expanded context these
lower their respective energies by 14.091344 and 3.695339 kcal/mol. Giving La
access to its old proposal reduces the expanded-only pool's apparent ordering
from +3.472171 to +1.265264. Both c5 metals select its expanded La proposal.
Excluding the unfavorable old q9 proposal would hide this real search sensitivity.

At the matched old q9 La geometry, motion changes the paired contrast by
**+12.086147 kcal/mol** in original context and **+12.111762** in expanded context.
The context–motion interaction is only **+0.025615**. This large response is
therefore preserved under enlargement. It is not solely a solvent artifact:
the original response contains +8.517684 from MACE and +3.568463 from solvent.
At the origin, q9's context effects are −11.693344 MACE and +8.716626 solvent,
giving −2.976718 net; c5's net context effect is only +0.330984.

Context expansion changes real fragment membership, caps and cavity together.
These results do not uniquely identify an electrostatic, mechanical or learned-
potential cause. The selected fourth torsion also differed in the original
searches; these were cross-scores of their real saved proposals, not new searches.

**Recommendation:** retain the common-pool comparison and investigate why the
old q9 La movement obtains its large relative stabilization. The two-anchor
positive result does not yet support promoting enlarged-context accommodation
or restoring the failed corrected accommodated calibration.

## Independent checks and actual work

[REVIEW.py](REVIEW.py) imports no scientific executor or integrator. It checks all
40 component pairs, all 160 solver execution receipts and corresponding raw
stdout energies, all atom populations and charge closure, model receipt energies,
input/artifact hashes, and exact Ca/La paired coordinates. The largest numeric
difference between archived solvent XYZ and replayed MACE XYZ is
5.551115123125783e-17 Å, below the previously declared 1e-12 Å mapping tolerance.
Both original file identities remain recorded; different byte hashes are not
misreported as identical files.

All 40 solvent differences pass the separate v2 component gate. The maximum
independent discrepancy is 1.245193743670825e-7 Eh; maximum all-atom population
error is 9.19400207077814e-8 e. All 80 original absolute-energy comparisons still
fail. The v2 component check does not retroactively admit those absolute totals.
Actual output samples show converged iterative solutions and energies matching
the collected full endpoints.

The independent calculation reproduces every integrated component contrast,
selected geometry, accommodation work, context/motion interaction and two-anchor
ordering within 1e-7 kcal/mol. [RESULT.v2.json](RESULT.v2.json) pins the reviewed
integrator result, SHA256
`1845da30da4fa9d828eb8becf748749d437b07b5317ece9b4862157173deb1c7`.
The first independent recount is retained in `RESULT.json`.

Measured factorial process CPU: **486.153929 s** for 80 primary xTB calls and
**1960.283493 s** for 80 independent tblite checks. MACE records 20 new calls
(16 cross-context plus four replay checks), **590.695933374 CPU-s** and
113.883849568 s wall including model load. Its 24 archived endpoint reuses are
separate from those calls. The separate four-call dispersion diagnostic and
preparation costs remain in their source reports. These are process measurements,
not Slurm allocations or a whole-pipeline latency benchmark. **This review ran
zero molecular calculations.**

## Reproduce without molecular execution

From the repository root, with a fresh output destination:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python diagnostics/orca_parameter_audit_20261012/context_comparison_review/REVIEW.py \
  --mace diagnostics/orca_parameter_audit_20261012/context_MACE/RESULT.json \
  --solvent workspaces/orca_parameter_audit_20261012/context_gfn2_component_v2/FINAL_COLLECTION.json \
  --comparison diagnostics/orca_parameter_audit_20261012/context_comparison/RESULT.json \
  --output /tmp/nikasha-context-factorial-review-20261012.json
```
