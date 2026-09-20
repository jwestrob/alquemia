# Independent numerical review

The water-basins agent independently recomputed the six native-tight corrections
from raw final energies, finding exact agreement with the reported maximum score
change0.025544881546637945kcal/mol and endpoint change0.025577379644665527.
The original0.10/0.20 criteria and actual printed TolE values were verified.

All48 task coordinates are byte-identical to source, with unchanged charge/spin
and exported GFN2 parameters. The25 successful property files independently
confirm electron/atom counts; max|charge|0.99426e and maximum charge-closure error
1.947e-8e. Earlier xtbw and simple-MORead controls actually use SAD/ignore MOInp.
All eight explicit-block GBW checks use actual MOREAD, but only one converges;
no complete solvent or Ca/La pair exists. Reviewer found no material discrepancy
or overclaim in the native-precision-only qualification. No files or calculations
were changed by the review.
