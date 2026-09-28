# Existing La response constrains interpretation of the small new signal

Matched source, La A/B XYZ and point-charge hashes are exactly identical between archived native r2SCAN-3c and the new PBE0-D4/lcecp1 calculation. No molecular evaluation was repeated.

Native-3c La A→B work: -0.1766543900 kcal/mol. New PBE0 recipe: -0.0435310966 kcal/mol. Recipe difference: 0.1331232934 kcal/mol. Both favor the same motion weakly, but their magnitude differs by more than the new +0.10778728 kcal/mol Dy-minus-La response.

This is **not** a measured uncertainty on the La/Dy contrast: native explicit-f Dy never converged, and method effects might cancel or reinforce between metals. It does show that the small single-source signal should not be promoted as method-independent specificity. Functional, basis/composite treatment and numerical settings differ; this comparison cannot identify a unique cause. No absolute cross-method totals are compared.

Use the already available physical force projections to choose a meaningful structural test. Do not launch a functional sweep or resurrect failed explicit-f jobs solely to resolve this small effect. The source reversal remains the more consequential challenge. Raw collection pins and exact numbers are in LA_RESPONSE_METHOD_COMPARISON.json.
