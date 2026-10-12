# Independent native-GFN2 runtime review — 12 October 2026

**The reported 60-control result is supported by the actual artifacts.** An
independent parser, without importing the execution collector, re-read every
stdout/property pair and verified 720 input/artifact hashes. It independently
reproduced 52 admitted endpoints and the exact eight retained density failures.
No molecular program or new upstream solver was run for this review.

## Results checked

| Claim | Independent result |
|---|---|
| Original eight Ca/La controls | All eight reproduce their pinned archived native energies and full atomic populations |
| Corrected full La file imported into Ca | Both media give bit-identical parsed energy and all 222 charges to the corresponding current default |
| Corrected transition-metal aquo controls | All 14 pass against both archived xTB and tblite: 28 accepted comparisons |
| Unchanged aquo export/import | All 14 pairs have exactly equal parsed energies/populations; only ten pairs have both endpoints admitted |
| Closed-shell aquo controls | All eight native endpoints admitted; seven cases pass both upstream comparisons; Zn water fails both population comparisons |

Corrected aquo maximum absolute errors are **5.278201342662214e-8 Eh** and
**7.279609413046728e-7 e**, under the declared 1e-6 Eh / 1e-5 e gates.
Zinc water's largest charge error is **8.742845572728974e-5 e**, despite its
energy difference of only 2.0978340842248144e-8 Eh against tblite.

The eight failures are the default and unchanged-roundtrip versions of Cu2
vacuum, Cu2 water, Fe2 water and Fe3 water. Each fails both printed density
criteria despite normal termination. The largest is Fe2 water: maximum density
change 1.1752e-4 against 1e-5, RMS density change 7.7515e-6 against 1e-6.
Agreement between two inadmissible outputs does not rescue those endpoints.

## Provenance and scope

[RESULT.v2.json](RESULT.v2.json) retains all individual residuals, receipt hashes,
archived comparisons and failures. Its initial smaller recount is retained in
`RESULT.json`. Both describe the same completed calculations. The source report
reviewed has SHA256 `2f139e12a174d61ccfa171633fcabb1eb989971ad23e9c8f11ce0433c57964a3`;
the archived aquo source table has SHA256
`296404712010b6072522fdf30c8aabd93a73c93590587b9b4e02f8832c1735c7`.

The independent comparisons use the supplied accepted **archived macOS xTB and
tblite records**, with matching geometry hashes, charges and multiplicities.
They do not imply fresh Linux executions of either upstream solver. The eight
native replication controls are separately compared against archived native
stdout and full-precision property files.

Two metadata conventions need care when extending admission: the explicit Ca
import receipts retain copied null `parameter_source` metadata; corrected aquo
receipts use that field for their original correction base. Actual imported
identity is established by task-local file hashes and runtime import lines.
The current Ca/La guard handles the first case. A future aquo registry needs an
explicit base/import adapter, rather than treating the ambiguous field as the
imported file or rewriting completed receipts.

This qualifies the stated parameter/runtime controls, not protein prediction,
physical open-shell energetics, mixed-metal parameter unions or new builds.
The explicit Ca/La guard remains restricted to its actually covered conventions.

## Reproduce the independent check

From the repository root, with a fresh output destination:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python diagnostics/orca_parameter_audit_20261012/runtime_review/REVIEW.py \
  --workspace workspaces/orca_parameter_audit_20261012/runtime_qualification_v1 \
  --reported diagnostics/orca_parameter_audit_20261012/runtime/RESULTS.json \
  --output /tmp/nikasha-native-runtime-review-20261012.json
```

The script refuses an existing output path. It performs only file reads,
hashing and scalar parser/comparison operations.
