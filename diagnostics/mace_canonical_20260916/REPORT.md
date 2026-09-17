# Direct canonical PQQ MACE candidate: completed, rejected

**Both checkpoints fail calibration separation on the exact 25-member panel.**
The baseline remains unchanged. The apparent success of the earlier two-crystal
ordering did not establish a usable classifier across the canonical set.

| Method | Calibration gap (min La minus max Ca), kcal/mol | Result |
|---|---:|---|
| Released DFT baseline | +8.60293200268643 | Separated |
| MACE-POLAR medium + frozen OBC-II, primary | -5.816631373832934 | Overlap |
| MACE-POLAR large + frozen OBC-II, sensitivity | -1.6092710972880013 | Overlap |

The MACE extrema in both checkpoints are La-class A0A3F2YLY8 and Ca-class P38539.
No bands are released under the frozen max-Ca/min-La rule. All three transfer
scores were computed but their classifications remain **unavailable because
calibration failed**. This is not three wrong transfer classifications. No
threshold was shifted, no model selected after scoring, and no sign reversed.
Unrounded endpoint energies, components, groups, source pins and baseline
comparators are in [result.json](result.json).

## What ran

Exact original canonical cores, caps, hydrogen coordinates, PQQ(3−), singlets
and zero waters: 25 calibration structures plus consumed 1H4I, 4MAE and 1KB0.
MACE jobs 1200776/1200779 and solvent jobs 1200781/1200782 all completed.
There were **108 new MACE and 108 new GB evaluations**, four compatible cached
medium endpoints reused per method, no scientific failures and no new DFT.
The original audit/preflight failures happened before scientific execution and
remain preserved. Preparation v2 fixes packaging only; it changes no input atoms.

The frozen model is direct vacuum MACE plus its endpoint-monopole OBC-II reaction
energy. It is not a correction to CPCM. Both Ca and La use the declared common
1.8 Å radius; no radius or dielectric was tuned. Existing solver and analytic
kernel checks are prerequisites, not newly repeated claims of chemical accuracy.

Seven real-fixture canonical parser/algebra/runner tests passed in 5.416 s.
The existing runner regression passed 10 tests in 8.173 s; two checks requiring
the isolated MACE installation were explicitly skipped in that Python environment.
Actual GPU calculations above are reported separately from software tests.

## Cost and feasibility

Recorded incremental allocation: **1,212 GPU-seconds, 19,392 allocated
core-seconds, 1,236.464 reported actual CPU-seconds**. One A5000, 16 CPUs and
64,474 MiB host allocation per job. Median model inference per La/Ca pair was
2.179 s (medium) and 3.563 s (large); launch, load, validation and collection
make the complete jobs longer. Newly computed solvent pairs took about 0.009 s
inside the solver. Cached medium values retain their original timings.
Peak MACE allocated GPU memory was 473,664,512 / 904,410,624 bytes for
medium / large; maximum worker host RSS was 1,673,660 KiB.

[costs.json](costs.json) includes actual job and endpoint receipts. Preparation,
failed preflight, software setup and local checks were not fully profiled and
are explicitly unavailable, not zero. Earlier reused calculations retain their
original development costs. No matched-hardware production DFT speedup is claimed.

## Interpretation and next action

- Numerically executable: yes, with the previously validated numerical backend.
  The solvent radius and frozen-charge approximation remain physical limitations.
- Predictively useful as this classifier: **no**; even the required calibration
  separation fails. No broad-affinity or incremental-information claim follows.
- Affordable: small-core inference is seconds and memory is modest. Affordability
  does not compensate for the failed prediction criterion.

Both 1H4I/P16027 and 4MAE/I0JWN7 share calibration accessions: they are structural
transfer observations, not independent proteins. All cases are retrospective;
PQQ functional class differs from direct affinity, and motif/core charge already
separate the canonical panel. These limitations are recorded before interpretation.

**Recommendation: retain baseline and abandon this frozen direct MACE+GB
classifier.** The broader MACE goal remains active. The transferable result is
that a successful two-structure ordering was insufficient; full-system learned
charge response and small-core vertical scoring have now both failed declared
predictive tests. Validated local curvature remains a separate useful component.
No production default, reference, benchmark label or immutable experiment changed.

Runnable read-only audit, preparation, collection and report commands are in
[COMMANDS.md](COMMANDS.md). The completed report lives at
workspaces/mace_canonical_20260916/report_v2/; report_v1 is preserved and differs
only in how unavailable calibration is displayed. Next read-only command:

```bash
python -m unittest discover -s tests -p test_mace_canonical.py -v
```
