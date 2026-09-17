# Denser source sampling: follow-up complete, full source screen still fails

Job1201026 completed8fits and8native validations. All native numerical checks,
potential screens and Ca−La vector-field screens pass. **All four alpha endpoint
field screens fail**, and alpha6IP9 exceeds the1kcal paired-U0 flag. The declared
single sampling follow-up is closed; no regularization sweep or added multipole
orders follow automatically. Baseline/default and previous results are unchanged.

| Preparation | New-probe field error Ca / La | Paired U0 error, kcal/mol |
|---|---:|---:|
| GGR extended |7.383% /8.973%|−0.150486|
| GGR connected |8.443% /9.631%|−0.291933|
| Alpha1F6S |13.250% /13.903%|−0.737790|
| Alpha6IP9 |12.022% /13.517%|−1.231046|

These probes are on the opposite side of each original environmental atom.
The previous positive-offset observations were explicitly used as development
training data. V1andV2 thus have different spatial validation observations;
do not turn their aggregate errors into a controlled same-probe improvement
claim. Both versions preserve the fixed molecular states and biological strata.

Maximum fitted local dipoles reach6.383eÅ; maximum absolute charges reach4.110e.
The coefficients approximate a whole quantum distribution and are not validated
atomic populations. A good summed field cannot automatically survive applying
atom-wise native damping or covalent exclusions to these coefficients.
This is another reason not to pass the source representation into a scorer.

The calculation gives no new La/Ca classification. U0 is the previously declared
bare diagonal response diagnostic, not a native induced energy. Source-field
errors do not uniquely diagnose earlier discriminator failures. The useful
engineering result is that explicit saved-density fields are affordable,
whereas these compact fits have not earned their use for induction.

## Real execution, software checks and cost

Eleven tests pass in11.623seconds across both versions, including actual native
receipts, charge closure, derivative signs/tensor ordering and use of actual
native observations as V2training data. Integration was explicitly skipped
before outputs existed. No simulated native scientific output was substituted.

Job244wall-seconds×8CPUs=1952allocated core-seconds,804.084actualCPU-seconds,
zeroGPU. No newSCF/DFT, MACE or force-field energies. Exact fit/utility/preparation
timings and memory are in DISTRIBUTED_SOURCE_SAMPLING_RESULT.json. The initial
preparation failed before any fit or native call because a local list shadowed
a function. Its artifacts are preserved; its elapsed/CPU times were not measured.
Recovery renamed that local variable with unchanged scientific inputs/model.
There were no failed fit/native attempts or extra scientific retries.

Protocol `vacuum_density_physical_charge_dipole_two_center_source_v2`.
Workspace `workspaces/mace_omol_20260917/distributed_source_v2_recovery_v1/`,
report_job_1201026; manifest SHA256
`f07f44bd2c41d3887a8ebeb6894a1eed9d438dc216b64cf73e55a35552a1108a`.
V1andV2together used16actualfits/16nativeutilities,3904allocatedcore-seconds,
1623.268actualjobCPU-seconds, noGPU. These are development costs, not a measured
production score. Fitted parameters/singular spectra/source hashes remain available.

Recommendation: retain baseline and stop this charge/dipole fitting trial.
Investigate direct quantum-field input to the existing native polarization
solver, with its energy expression verified first. No new scored model is
declared by that investigation; native d/p scaling and GK accounting must be
respected. No fitted correction is filled into an unavailable score.
