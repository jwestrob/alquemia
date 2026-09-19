# First task-specific classifier experiment

**The classifier and scoring interface work. This frozen first feature set does
not improve accuracy.** Structural features preserve PQQ performance; adding
the two MACE summaries introduces one held-out error. No model was retuned.

## PQQ: hold related proteins out together

| Fitted classifier | Correct / 25 held out | Ca correct / 14 | La correct / 11 |
|---|---:|---:|---:|
| DFT contrast | 25 | 14 | 11 |
| DFT + structural features | 25 | 14 | 11 |
| DFT + structure + MACE summaries | 24 | 13 | 11 |

Four sequence-homology groups were held out in turn, with 18, 5, 1 and 1
canonical cases. Accession duplicates and transitive >=50% identity groups
remain together. Scaling and coefficient fitting use training rows only.
The three classifiers share exactly the same folds and rows. This is development
testing within the PQQ fold family, not independent cross-family validation;
all source cases had already been inspected during earlier development.

All three final fits correctly replay 1H4I and 4MAE. Both crystals share a group
with training sequences and remain consumed transfer/replicate checks. 1KB0 is
unsupported, retained in the three-transfer denominator. The released DFT bands
and default scorer remain unchanged; these fits have their own declared logit
decision rule and do not modify any historical reference.

## What the added MACE information did

Q9Z4J7, labelled Ca, becomes La under the held-out MACE-feature model. Its
logit is +0.02342588740702891. The 6–12 A response contributes +0.33862635668857816;
the intercept and remaining terms together contribute -0.31520046928154925.
This explains the fitted decision algebraically, not a causal physical mechanism.
The DFT-only and DFT+structure fits retain the correct held-out call.

The local water-probe accessibility feature is zero throughout the canonical
PQQ panel and supplies no information there. Its inactive coefficient is fixed
to zero using training variance, not selected from classification results. No
probe radius, shell boundary, regularization or threshold was adjusted afterward.
This tests these particular readout summaries, not every possible MACE representation.

## Direct affinity: insufficient independent training groups

The fully matched site-level data have three GGR structures (Ca direction) and
two alpha-lactalbumin structures (condition-qualified La direction). They are
**two biological groups**. Holding either out leaves one training class, so all
five held-out predictions are unavailable. This is not a measured 0% accuracy.

Each model fits all five training structures correctly, including DFT alone.
That is not a correction of the known failures: it cannot be separated from
learning the two training groups. Exported affinity models explicitly state
`not_evaluable_missing_labelled_groups`. Protein-level or supporting labels
were not promoted to site-level affinity truth to increase sample count.

## Delivered and tested

- `scripts/site_classifier.py`: inventory, feature extraction, grouped evaluation,
  export and score operations. Eight feature coefficients at most, fixed L2
  logistic regression; no checkpoint training or hyperparameter search.
- 32 supported feature rows, exact source/readout checks, six exported research
  models with evaluation status, and all held-out predictions/coefficient records.
- Six real-fixture tests pass in 3.465 s: rigid transforms, corrupted coordinates,
  objective gradients, replicate weights, homology grouping, saved-fit replay,
  protocol rejection and unlabelled prediction handling. These are software checks.
- The actual exported DFT+structure model scores 1H4I as Ca with logit
  -0.9276845019065618. This is a replay, not a new independent success.

Final feature extraction: **13.237680656835437 s wall / 12.908226509 CPU s**.
Grouping and all fits: **1.039798479527235 s wall / 0.911060191 CPU s**.
No new DFT evaluation, MACE forward, GPU allocation or Slurm job. Earlier
development extraction/tests are separately retained; these timings are not
claimed as the total agent-development cost or a production speed benchmark.

## Next accuracy requirement

Do not spend another campaign retuning these features on the same PQQ cases.
For the affinity classifier, the immediate missing test is another independently
labelled La-favoring family and another Ca-favoring family with matched inputs.
The existing evidence ledger's PqqT and aqualysin candidates are worth resolving;
their preparation/assay mapping gates must be settled before using them as labels.
No expansion, new calculation or altered classifier was launched here.

[Commands](COMMANDS.md), [exact results](RESULT.json),
[all features and predictions](FEATURES_AND_PREDICTIONS.tsv),
[frozen scope](AGREEMENT.md). The completed experiment is not an accuracy gain
and does not authorize production promotion.
