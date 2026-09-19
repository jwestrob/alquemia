# Fixed extension, before inspecting transferred logits

Parent approved 2026-09-19: reuse the already fitted all-25 PQQ DFT-only and
DFT+radial coefficients, without refitting, and report alpha-minus-GGR logit
differences on all four combinations of two alpha structures and two GGR
representations. Calculate coefficient times feature difference divided by the
saved training scale; intercept and training mean cancel. Do not produce an
absolute class or probability across the different targets/preparation protocols.
Positive difference orders experimentally La-associated alpha above
Ca-associated GGR. Labels and biological groups come directly from the earlier
`site_classifier` feature/evidence ledger. This is exploratory cross-target
transfer on consumed cases, not independent validation or new calibration.

Reporting correction: initial feature v1 manually reversed the direct labels;
the numerical features and all PQQ fits were unaffected. Preserve v1 with an
explicit superseded-metadata notice, and write corrected v2 metadata from the
authoritative ledger. Do not change coordinates, numerical features or models.
