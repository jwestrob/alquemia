# Quantum source fields: distant agreement hides a near-field weakness

All8saved-density potential calculations completed as1201017. Numerical
resolution passes; **both the original CHELPG and projected monopole models
fail the declared all-site endpoint field screen**. No new discriminator or
correction is reported. The production baseline and prior scores are unchanged.

## Actual result

Every environmental atom outside the fixed QM source support was tested:
4640GGRextended,4587GGRconnected,1880alpha1F6S,1843alpha6IP9 atoms per endpoint.
Each received12observation-point offsets for two central-difference spacings.
No nucleus moved and no SCF, new charge fit, MACE or force-field energy ran.

All-site polarizability-weighted relative field errors are21.4–33.5% for the
original fitted charges and27.9–40.7% after their physical cap projection.
All8endpoint screens fail. In the predeclared ≥3Å stratum the respective
ranges are1.98–5.72% and2.61–7.90%; all8pass. These are separate strata:
the close sites remain part of the failed all-site result.

The Ca-minus-La field differences are more accurate:3.28–4.81% fitted and
4.54–5.83% projected over all sites; all4paired vector screens pass. However,
that cancellation alone does not qualify an induction model:

| Consumed preparation | Projected Ca−La vector error | Paired projected−exact U0 error, kcal/mol | ≥3Å U0 error |
|---|---:|---:|---:|
| GGR extended |5.63%|−1.792884|−0.000243|
| GGR connected |4.54%|−0.852195|+0.087725|
| Alpha1F6S |5.81%|−3.141362|+0.045973|
| Alpha6IP9 |5.83%|−2.728425|+0.058780|

`U0=-0.5 Σ alpha |E|²` is a **bare diagonal diagnostic**, not the native
AMOEBA induction energy or an environmental correction. It omits protein
permanent fields, mutual polarization, solvent, damping and covalent scales.
Three of four paired projected U0 errors exceed the frozen1kcal flag. Original
unprojected fits already yield errors−0.885/−0.567/−2.526/−2.231kcal in the
table's order, so cap projection is not the only representation issue.

The transferable point follows directly from the quadratic expression:

```
U0_Ca − U0_La = −0.5 Σ alpha (E_Ca + E_La) · (E_Ca − E_La).
```

A common endpoint field error can couple to the physical Ca−La difference;
accurate differential fields do not guarantee accurate differential quadratic
energies. The earlier exterior-potential/direct-coupling pass did not test
this. Conversely, this diagnostic does not identify the cause or magnitude
of error in a complete damped/covalently scaled model that has not run.
Close covalent boundary sites need explicit treatment; do not simply remove
them because the far stratum looks better.

## Numerical checks, cost and limits

Maximum vector change on halving the observation-point spacing was
1.454029979e−8atomic field units, versus the1e−6tolerance. Largest change in
paired U0 was1.132347037e−6kcal/mol, versus0.01. These errors are far below
the representation effects. Four tests pass in5.788seconds, including native
output/receipt/vector replay; integration was explicitly skipped before output
availability. The point-charge derivative test is algebra on real coordinates,
not fabricated quantum evidence.

The job took203wall-seconds with8CPUs:1624allocated core-seconds,
711.484actualCPU-seconds, zeroGPU. Native utility wall times sum to
702.264582902seconds;699.04utilityCPU-seconds. Largest individual utility RSS
318076KiB; sampled batch MaxRSS1511192KiB. Input-preparation wall time is in
the manifest; its CPU use was not measured. All outputs/wavefunctions/receipts
are retained. There were8utility calls and no failed attempts.

Protocol `normalized_vacuum_QM_electric_field_diagnostic_v1`; workspace
`workspaces/mace_omol_20260917/qm_electric_field_v1/`, report_job_1201017.
Manifest SHA256 `0b639b8502fdade3b8dd7727563908cbdf1d2237ae4008b8cfe329ce872d0673`.
Compact pins/costs/results: QM_ELECTRIC_FIELD_RESULT.json. Vector arrays,
individual-site coordinates and exact native potentials permit further audit.

Judgment: native field evaluation is numerically credible and inexpensive.
The current monopoles are not qualified as the complete induction source.
No biological improvement or broad generalization was tested; these remain
the same2consumed biological groups with qualified alpha evidence.

Recommendation: retain baseline and continue the source-representation work.
Test one inexpensive distributed charge/dipole representation against separate
spatial probes before any new polarization score. Do not rerun a large panel
or replace this failure with a favorable subset. Existing MBIS monopoles were
already examined historically and were slower/worse for some potential/coupling
tests; do not present MBIS as an untried automatic fix. Native ORCA also supports
MBIS atomic multipoles, but that capability alone does not demonstrate improved
fields or affordable extraction here. See the [native population manual](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/population.html#mbis-charges)
and historical diagnostics/density_embedding_20260916/CHARGE_FIT_RESULT.md.
