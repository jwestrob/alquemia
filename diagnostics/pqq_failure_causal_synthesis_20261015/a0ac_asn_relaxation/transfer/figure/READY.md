# Individual-context comparison figure — source ready, rendering unrun

The thin `PLOT.py` CLI consumes the **actual completed saved reporter outputs**, not a synthetic dataset or newly computed score. It produces editable SVG and PDF plus a 220-dpi PNG. No numerical/render execution, imports, tests, binding or submissions occurred during preparation. Gauss's independent static schema review found no blocker; runtime and visual checks remain unrun.

## Exact future command

After the actual `asn_transfer_report_v1/execution` files exist and its reporting `EXECUTION.json` is complete, root may run from the repository root:

```bash
sbatch diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/transfer/figure/RUN.sbatch \
  "$PWD/diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/transfer/figure/PLOT.py" \
  e3b04f953f2f0c6fd22fa2fb929171b32d774f5e571bd062af943a33d886215a \
  "$PWD/workspaces/pqq_failure_causal_synthesis_20261015/asn_transfer_report_v1/execution" \
  "$PWD/workspaces/pqq_failure_causal_synthesis_20261015/asn_transfer_figure_v1/execution" \
  operational
```

One CPU, standard queue, no memory/time/priority override; root owns submission and completion wake. The explicit script path avoids Slurm's spooled-script location. Output must be a new directory. No watcher, binding engine, molecular code or source/scoring mutation is included. The optional `mathematical` policy uses the same recorded report fields and fixed bands in a separately named output directory; it does not select whichever policy looks better.

Direct CLI inside that allocation is:

```text
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python3.11 ABS_PLOT.py --script-sha SHA --report-dir ABS_REPORT_DIR --output ABS_NEW_OUTPUT_DIR --policy operational
```

## Figure design and population

- One row per **257 distinct source/context pairs**, grouped by recorded biological protein identity. All **260 memberships** are read; the three exact shared contexts are drawn once only after their identity, scores, selections, outcomes and regression flags agree. Both membership IDs remain in the figure-data files.
- Every individual common3 and common6 score is drawn: hollow versus filled conditioning-specific markers, with a connector. Ca, La, two-Ca and unspecified conditioning have explicit legend entries and row labels. Duplicate short row labels get deterministic suffixes; the full source/context IDs remain in editable data and SVG element IDs.
- Recorded regressions or lost coverage have red connectors/labels and `!`. No group averages replace individuals. Missing scores get an `x` in the separate **NA 3 / 6** gutter, never a fabricated zero. Partial scientific coverage is permitted if the saved reporter completed with the full denominator.
- The only score transformation is subtracting the unchanged canonical **Ca_max** from every score and both boundaries. All protein panels share one x scale. The two saved boundaries, shaded support/inconclusive regions and raw offset are retained in `FIGURE_DATA.json`. Different common3/common6 bands fail explicitly; no thresholds, probabilities or decisions are fitted or recomputed.
- Titles/caption state consumed development evidence and correlated structures. The current plot includes the crystal/probe memberships; the 18 superseded historical primary records remain in the report's `SUPERSEDED_PRIMARY.tsv`, outside the declared current 260-membership population. This is not independent biological validation or production promotion.

The figure uses four columns of protein panels with enough height for every source row; it never samples points or clips extreme scores. SVG text remains text (`svg.fonttype=none`); PDF uses embedded TrueType text (`pdf.fonttype=42`). Main artifacts are `COMMON3_COMMON6_operational.{svg,pdf,png}`, exact input/output pins and all individual rows in `FIGURE_DATA.json`/`.tsv`, and the Slurm `EXECUTION.json`. Rendering failure retains an execution receipt and does not synthesize missing data. Root should visually inspect the actual first render before delivery.

## Exact code/environment provenance

| Artifact | SHA256 |
|---|---|
| PLOT.py | `e3b04f953f2f0c6fd22fa2fb929171b32d774f5e571bd062af943a33d886215a` |
| RUN.sbatch | `a30c1df5e465cbac1d33ee8648bf258f8204c9fb7b91db57164c74ff6fa4f718` |
| Required actual REPORT.py producer | `5c526e22a182ae24ef6c0f6cfefdd506963765b4eafbcefab8fc990df0004f99` |
| `/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python3.11` | `53244dfc692ccbfd6a065f053dfa2026e2f3e4755beef84d600c8f7dbc9b5b27` |
| Existing Matplotlib 3.10.9 package METADATA | `a37f68bc5fec598f62d7641cdde0b10e137514400edd912cca88836177798559` |

The plotting environment is the one that actually rendered the A0AC saved-summary figure in job5030. Its executable and package metadata pins are checked before Matplotlib import, which occurs only inside Slurm. Actual Matplotlib version/source, font file, input table hashes, display limits and outputs are recorded at runtime. This does not claim a complete lockfile for every transitive environment dependency.

Style conventions reuse `diagnostics/nikasha_manuscript_20260922/render_structural_comparison.py` (SHA `589ed7e616aa3fadbdc4e6f021b3c2b526f3062533b6e1a4415d8a157a62683e`): DejaVu Sans, editable text, white background, restrained spines and the existing teal palette. Its count-bar logic is not reused because this figure must expose every structure and regression. The complete reporter schema and actual canonical unchanged-band report were read before implementation; final transfer outcomes were not invented or assumed.
