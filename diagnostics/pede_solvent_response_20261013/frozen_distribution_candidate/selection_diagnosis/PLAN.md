# Saved-matrix selection diagnosis

## Question and fixed scope

Separate the completed frozen-distribution candidate's raw score changes into
(1) changing the energy expression at the baseline's selected geometries and
(2) selecting geometries again under the candidate. Separately account for the
movement of the independently frozen decision bands. This explains transitions;
it does not choose a new classifier, tune a threshold or estimate an affinity.

Input is the completed transfer collection with SHA256
`82981de62dfa2f2ff5eaedd273e38d8b20bf36034dc39110c104908599470d25`, beneath
`workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/`.
Include every declared record and all three policies: operational, mathematical
and origin. Preserve unavailable rows, cohort membership, source identities and
the existing 125-source recovery mapping. Overlapping cohorts are not independent
observations and are not pooled into one accuracy claim.

## Exact decomposition

All saved scores and works are already kcal/mol; convert no units again.
For each source/policy let old and new selected geometries be `o_M` and `n_M`.
Using the candidate's saved per-metal works relative to its origin:

```
R_fixed = R_new + [W_Ca(o_Ca) - W_Ca(n_Ca)]
                    - [W_La(o_La) - W_La(n_La)]
delta_model = R_fixed - R_old
delta_selection = R_new - R_fixed
```

`R_fixed` is a counterfactual using the new energy expression at the old
selected geometries; it is not an optimized state or a proposed hybrid rule.
An independent reconstruction uses the new origin contrast plus the two old-
selection works. Preserve each metal's reselection work separately.

For expected class define `s=-1` for Ca and `s=+1` for La, and use its own
correct-class boundary `b` (`Ca_max` or `La_min`). Signed margin is
`m=s*(R-b)`: nonnegative means the score reaches its expected-class band;
negative can mean inconclusive or wrong, distinguished by the saved decision.

```
delta_margin = s*delta_model + s*delta_selection - s*(b_new-b_old)
```

Report raw and counterfactual scores, old/new selections, saved and counterfactual
decisions under their stated frozen bands, signed margins, all three margin
contributions and closure residuals. Do not infer an error mechanism from an
arithmetic contribution alone. The candidate still changes both solvation
response and its calibrated scale.

Require absolute closure at most **1e-8 kcal/mol**, declared before running this
analysis, for raw-change and margin accounting, selected-work reconstruction,
and the independent origin reconstruction. Missing matrix works/bands remain
explicitly unavailable; no zero substitutes. Inconsistent available data fail
the analysis visibly. Verify the collection and its canonical/inventory pins;
require exact equality to the already frozen candidate and baseline bands.

## Execution and outputs

Only standard-library saved-result analysis; zero molecular calls, optimizations,
recalibrations or tests. Execution requires a real Slurm job. Root submits one
CPU on the standard partition, with no new memory/time limits. Code/metadata
inspection and AST checks on login do not execute the scientific analysis.

Explicit command (root chooses a fresh output directory):

```bash
python diagnostics/pede_solvent_response_20261013/frozen_distribution_candidate/selection_diagnosis/DIAGNOSE.py \
  --collection workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/FINAL_COLLECTION.json \
  --sha 82981de62dfa2f2ff5eaedd273e38d8b20bf36034dc39110c104908599470d25 \
  --output workspaces/pede_solvent_response_20261013/selection_diagnosis_v1/execution
```

Outputs: `RESULT.json`, flat `ROWS.tsv`, compact `REPORT.md`, and a success/failure
`EXECUTION.json`. These are diagnostic tables only. Production and the original
collection, classification decisions and six unavailable triples stay unchanged.
