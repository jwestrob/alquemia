# Archived collection schema correction

Analysis job3338 failed in one second before writing a scientific result. The
file named CANONICAL_COLLECTION.json contains135 rows, including transfer and
separate-probe roles; the declared analysis uses its25 `canonical_calibration`
rows. The initial reader did not filter by role and encountered an unavailable
policy record. Version2 selects the exact predeclared canonical role and handles
explicit null policy pools. No reference band, score, outcome or scientific
denominator was changed. The failed receipt and exact executed code remain in
calibration_v1; version2 uses a new output directory and frozen implementation.
