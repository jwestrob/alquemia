# Disconnected sodium consistency check

Execution: complete; 20/20 endpoints.
Numerical gate: True. Spectator consistency: False.
Modified development directions: False.

| Case | Original R_coord | With spectator | Change, kcal/mol | Consistent within0.1 |
|---|---:|---:|---:|---|
| GGR_1GLG | 62.71484481271285 | 58.633528901355916 | -4.081315911356931 | False |
| ALPHA_1F6S | 99.86840743209439 | 85.10310583241724 | -14.765301599677144 | False |
| ALPHA_6IP9 | 90.80313269254084 | 76.59966790583779 | -14.203464786703051 | False |
| PQQ_1H4I | 101.61412567754543 | 101.61412572049906 | 4.295362998618657e-08 | True |
| PQQ_4MAE | 110.64796003887037 | 101.18908554304907 | -9.458874495821291 | False |

| Contrast | Original | With spectator | Change, kcal/mol | Direction retained |
|---|---:|---:|---:|---|
| PQQ_4MAE minus PQQ_1H4I | 9.033834361324935 | -0.42504017744998634 | -9.458874538774921 | False |
| ALPHA_1F6S minus GGR_1GLG | 37.15356261938154 | 26.469576931061326 | -10.683985688320213 | True |
| ALPHA_6IP9 minus GGR_1GLG | 28.088287879827995 | 17.966139004481875 | -10.12214887534612 | True |

The constructed Na+ is at least10000A from every original atom and has zero model graph edges. Original coordinates/protonation/waters remain identical. Total charge increases by one; electron parity remains valid.
The checkpoint receives total charge and cannot constrain independent fragment charges. Interpret sensitivity as a representation issue for the intended separated ionic state, not as measured sodium effects or an exact ground-state error.
All cases are consumed development. No threshold adjustment, new biological label, absolute affinity claim or production promotion.

## Cost and conclusion

Job1200823 completed20new forwards with no failures:667GPU allocation-seconds,
10672allocatedcore-seconds,612.926reportedactualCPU-seconds. Native inference
total301.920377seconds; peakGPUallocation11,537,681,408bytes. Preparation
514.63wall/309.90CPU seconds and reporting9.35wall/11.44CPU seconds are separate.
Most preparation overhead revalidates archived qualification artifacts.

Numerical implementation: passes. Intended disconnected-fragment consistency:
fails. Predictive usefulness: original development ordering remains observed,
but is fragile to global charge conditioning and cannot justify promotion.
Affordability: one A5000 handles whole chains in tens of seconds per endpoint;
this does not rescue the representation problem. No newDFT,solver,training or
force evaluations.

Read CHARGE_EMBEDDING_AUDIT_REPORT.md for the exact checkpoint degeneracy
that explains MxaF's invariance. The raw immutable collection's inherited
contrast `pass` refers to the original direction; this report names original
and modified directions explicitly and separately from the consistency gate.
Three real preparation/corruption tests and six existing OMOL regressions pass.
The first new test run had a test-field typo (electrons versus all_electron_count);
correcting that assertion changed no scientific input or output.

Retain the production baseline. Continue the MACE goal with an explicit charge
representation strategy; do not repair these findings by adjusting thresholds.
