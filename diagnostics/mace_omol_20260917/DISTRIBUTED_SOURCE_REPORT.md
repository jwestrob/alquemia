# Charge/dipole source: substantial field improvement, incomplete qualification

All8fits and8native spatial-validation utilities completed as1201022. Separate
probe field errors fell from29–41% for the projected monopoles to8–14% for the
charge/dipole representation. **Six of eight endpoint field screens still fail**
the frozen10%/1e−4au criterion. No full-model energy or new classification ran.

| Preparation | Separate-probe field error Ca / La | Paired U0 error, kcal/mol |
|---|---:|---:|
| GGR extended |7.930% /9.824%|−0.399517|
| GGR connected |11.072% /13.036%|−0.087109|
| Alpha1F6S |10.327% /11.385%|−0.534054|
| Alpha6IP9 |12.126% /13.922%|−0.430531|

The paired Ca−La vector-field screens, all potential screens, numerical
resolution and all four1kcal paired bare-U0 flags pass. U0 remains an undamped
diagonal diagnostic, not an AMOEBA energy or score correction. Training field
errors are1.50–2.38%, substantially below the separate-probe errors: sampling
at atom centers alone does not establish local interpolation accuracy.

All fitted charges close to the formal QM charge within1e−9e; moment rotation
and rigid spatial replay pass1e−8au. Five tests pass4.700seconds on the real
outputs, including analytic field sign/tensor order versus spatial potential
derivatives. No surrogate native output was used. The integration tests were
explicitly skipped until the actual calculation completed.

Fitted moment magnitudes need scrutiny: maximum absolute charges range1.988–
4.378e and maximum local dipoles1.828–3.577eÅ. These are fitted expansion
coefficients, not uniquely measurable atomic populations. They cannot be
inserted into arbitrary atom-wise damping or covalent exclusions while
assuming the validated total source field is preserved.

## Measured cost and provenance

Job244seconds×8CPUs=1952allocated core-seconds;819.184actualCPU-seconds;zeroGPU.
Eight fits took22.247078992summed wall-seconds and9.506616981summed threadCPU-s.
Eight utilities took798.954364952summed wall-seconds and795.69CPU-s. Maximum
individual utility RSS318228KiB; sampled batch MaxRSS2498716KiB. Preparation
4.033089373wall/3.791939549CPU-seconds. No failed attempts, newSCF, MACE,
force-field energies, geometry changes or population analyses.

Protocol `vacuum_density_physical_charge_dipole_source_v1`. Workspace
`workspaces/mace_omol_20260917/distributed_source_v1/`, report_job_1201022;
manifest SHA256 `98820483e24705f632a993400f01aaabbd55ffe6551005c8a7db998040e67282`.
All source/validation densities, IDs, positions and executable hashes were
matched; fits were saved before their native validation calls. Full singular
spectra, coefficients and arrays are retained. Compact result/cost pins:
DISTRIBUTED_SOURCE_RESULT.json.

Baseline/default remain unchanged. This improves a physical representation
test, not demonstrated La/Ca predictive accuracy. Recommendation: retain the
baseline and make one declared denser-sampling follow-up, keeping the same
charge/dipole form, regularization rule and criteria. The first spatial
validation observations become development data in that new version; they
must not be reported as fresh validation again. No full-hybrid promotion.
