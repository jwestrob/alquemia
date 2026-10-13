# Actual partial-output test passed

Slurm3376 completed0:0 on one CPU. All six declared real-fixture checks passed: exactly one accepted native endpoint and seven unresolved; every incomplete DFT transfer/contrast remains unavailable; raw final-SP energy is preserved; saved corrected GFN2 response algebra is reproduced. The accepted scout also passed raw population, source-map, state/electron and native execution/hash checks.

The test used the pinned3370 corrected collection and no invented molecular values. Analysis0.170148s wall,0.124041s process CPU,24800KiB peak RSS; scheduler elapsed1s on one CPU. It made zero molecular calls. Exact evidence is `PARTIAL_VALIDATION.json`.

The complete eight-cell DFT/CPCM comparison remains unavailable until the dispatched jobs finish and the final collector runs. This test validates the saved-output join and missing-data handling, not a successful solvent model, a physical density response or improved classification. `COMPARE.py`, its snapshot and CONFIG remain unchanged after this test. Root owns the final afterany collector and wake.
