# Existing-score geometry diagnosis, before extraction

Question: does the PLM score spread reflect changes in the actual presented donor
geometry within the same fixed PQQ chemistry, and are those arrangements covered
by the established reference panel? This is model diagnosis, not an accuracy test
on PLM predictions.

Inputs: all 176 rows of the completed PLM XoxF cohort (137 scored), and all 25
canonical PQQ calibration cases. Retain every unavailable row. Read existing
source/carve manifests and endpoint coordinates, with hashes; no molecular
calculation, fold, optimization or source edit.

Extract all actual metal distances for each fixed-core donor role; PQQ N6/O5/
O7A/O7B separately; nearest and second-nearest carboxylate O for each acidic
role; cationic-partner distance; existing typed CN at 3.1 A; hydrogen-preparation
residuals where present; endpoint composition and charge. Preserve distances
beyond coordination cutoff. Compare reference-class and PLM predicted-band
distributions and report descriptive rank correlations with the existing score.
No fitted classifier or new threshold. Missing roles stay missing, not zero.

All PLM labels remain unknown. References were already used for calibration.
Any representative follow-up chosen from these results is explicitly development,
never blind validation. Freeze its identities and intervention before scoring.
Use the observed contrasts to distinguish a concrete geometry/state hypothesis
from merely asking another method to reproduce existing predictions.
