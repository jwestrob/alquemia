# Prioritize predictive usefulness alongside numerical qualification

Jacob agreed on 2026-09-18 to test whether complete methods improve benchmark
discrimination in parallel with qualification. He explicitly asked to reduce
provenance/consistency work that does not answer that question, and to email
actual discrimination gains, including better PQQ separation where possible.
This supersedes using the full solver qualification as a prerequisite for
computing an exploratory complete score. Preserve all failed gates and label
exploratory results clearly; no production/default change or threshold fitting.

## First complete comparison

Use the seven already consumed, prepared states: GGR 1GLG extended source,
2FW0 and 2FVY; alpha-lactalbumin 1F6S and 6IP9; parvalbumin 4CPV CD and EF.
Reuse their actual responsive-density DFT, direct coupling, MACE short-context,
AMOEBA environment and induced dipoles. No new high-level endpoints or new
biological labels. Replicate crystals and sites remain grouped by protein.

For each endpoint evaluate the complete frozen-response transfer already
derived and checked against native GK:

    C_model = E_model(P+(mu_d+mu_p)/2) - E_model((mu_d-mu_p)/2)
    R_new = R_old + (C_CPCM,Ca-C_GK,Ca) - (C_CPCM,La-C_GK,La).

All intrinsic QM, direct-density, vacuum-response and MACE terms remain fixed.
The environment-only transfer cancels within each Ca/La pair only after
confirming its identical native preparation/reference in the archived score.
The result is a frozen-response descriptor, not self-consistent polarization
or a complete binding free energy. No old aquo reference or calibrated bands.

Compute two resolutions alongside the biological comparison: primary18/974
and refined24/2030, with the existing conductor model, radii, epsilon78.3,
eta.1/shift0 and exact native multipole source potential. Eight continuum
solves per case (Ca/La times average/difference times two grids): 56 total,
seven independently allocated 64-CPU groups. Native GK cross components are
reused when exact inputs match, otherwise evaluated through the checked native
wrapper. No response iteration, DFT, MACE or forces. No project time/compute cap.

Use the isolated source-loop OpenMP build after a short full-array comparison
with a real stock source fixture. Source representation is independently
checked against Cartesian potentials during execution. The separate proposed
12-call thread-scaling experiment is deferred; measure actual pilot timings.
Source errors fail the affected calculation. Existing endpoint/grid failures
do not hide otherwise computed exploratory contrasts or become passing gates.

## Decision-relevant outputs

Report all seven old/new raw contrasts and all twelve alpha/parvalbumin-minus-
GGR margins, retaining the earlier 0.02 kcal directional screen. Report both
grids and each margin's numerical change; a positive but unresolved margin is
not a robust improvement. Alpha evidence remains qualified cross-study and
parvalbumin supporting cross-study; these are not twelve independent gold tests.
GGR's same-assay evidence remains Ca-favoring. Preserve every invalid result.

The main question is whether this improves the failed cross-family orderings
consistently across structures. Complete results may reject this route; do not
extend engineering merely to make it succeed. Keep numerical qualification,
predictive usefulness and measured cost as separate judgments.

PQQ separation is the next requested comparison if this complete correction
earns further work. First reuse compatible existing PQQ states; do not mix its
production CPCM baseline with this vacuum/trial-density hybrid algebra. Report
panel spread and transfer behavior, not calibration accuracy alone. No new
PQQ preparation or high-level run is included in this first seven-case manifest.
