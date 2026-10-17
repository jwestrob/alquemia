# LBT Cartesian implementation specification

Declared before new evaluations. This implements Modules 1–2 of PLAN.md under
the unchanged whole256 MACE-vacuum + fixed-FF-charge ALPB objective. It does not
resume the stopped clamped-release or projection diagnosis.

Four starts are exactly the admitted own-metal terminals of 12762/12765,
replayed through completed pool12767. All256 atoms move. No source-displacement
window, fixed donor, torsion subset, source spring or frozen bond length is used.
Physical state is charge−2, La multiplicity1/Dy6; the original force/energy
recipe and native units remain unchanged. No GFN, MM or DFT is added.

## Physical admission and state labels

The immutable original source defines graph, stereochemistry and inherited
contacts. Every covalent bond may change length inside0.65–1.35 times its pinned
native force-field equilibrium length. An inherited outside-range length is
retained as an allowed endpoint, without allowing further worsening. These are
broad connectivity guards, not a force-field energy or relaxed-bond restraint.
Reuse source chirality signs and peptide cis/trans basins; do not enforce
conjugated-group rigidity. All physical nonbonded contacts use the existing
1Å heavy /0.55Å H severe-overlap rule, excluding actual1–2/1–3 pairs and
preserving inherited shorter contacts at the existing1e−7Å roundoff rule.

For each target, separation from every peptide N/O by more than4.5Å labels
`departed_protein`. Separation from every declared original site oxygen by
more than4.5Å, with another peptide N/O closer, labels `departed_declared_site`.
Coordination-number changes alone are not rejected. These distance rules are
the same for both metals and sources, used only as state labels. LBT has no
retained waters or spectators. A decreasing state-changing endpoint keeps its
energy/gradient, terminates as a state change and is excluded from the bound
accommodation pool. Other evaluated departures remain in the trial records.

## Search

Limited-memory BFGS uses the existing direction helper, history8 and Armijo1e−4.
Remove only global translations/rotations from search directions, using the
Euclidean Cartesian metric and normalized-Gram relative cutoff1e−10. Report
both raw and rigid-motion-removed force norms; convergence uses the RAW full
gradient: maximum atom-vector norm≤0.10 and atom-vector RMS≤0.03kcal/mol/Å.

Initial maximum atom step0.025Å; adaptive trust radius between1e−6 and0.10Å.
Each line search has at most24 trials with multipliers1,1/2,… . An accepted
full trial with actual/predicted decrease>0.75 doubles the radius; a ratio<0.25
halves it. No cumulative displacement bound or energy-dependent basis change.
Backend failure stops the affected search. No automatic chemistry retry.

The7CCO/La scout stops at12 accepted steps or40 distinct new composite requests.
The pilot stops at80 accepted steps or200 requests per state INCLUDING the
scout's counts. Exact checkpoint coordinates, energy, gradient, L-BFGS history,
trust radius, visited cache and counters continue into its pilot run. Rejected
but evaluated points count; exact cache hits do not. Every accepted point has
an immutable checkpoint. Stops with residual force are finite outcomes, not
stationary minima. A geometrically blocked unchanged start remains available.

## Qualification and executable stages

`preflight` uses saved files only and requires no GPU/model initialization. It
replays all four terminal gradients from raw MACE and fixed-charge ALPB, original
source mappings, graph, native-state recipes and immutable source guards. Cache
key tests alter coordinate/state/protocol metadata only and call no energy.

`qualify` uses two predetermined directions at each of four exact starts:
(1) a deterministic coupled all256 Cartesian direction given by
`sin((3*i+a+1)*sqrt(2))`, with only global rigid motion removed; (2) opposing
motions along actual backbone A/8//C–A/9//N. Normalize each by its maximum atom
displacement. Evaluate±1e−4 and±5e−5Å, at most32 new composites total. Reuse all
central analytical gradients. Finite points must pass physical guards; rejected
slots remain unavailable, without choosing replacement directions.

Each MACE, ALPB and total central derivative must agree with its analytical
projection at BOTH steps to0.005+1e−5*abs(analytic)kcal/mol/Å. Step-to-step
agreement uses the same bound. No force or energy clipping, fitted tolerance,
retry with another step, or automatic extension. The fixed absolute resolution
is below the0.03 atom-vector RMS convergence target. This checks numerical
consistency of these directions, not physical accuracy or all possible motions.

`scout` requires the complete admitted qualification. `run --state SOURCE_METAL`
requires admitted scout execution/physical state, irrespective of preference
sign. For7CCO/La it continues the checkpoint toward the pilot's cumulative cap;
other states start their declared original terminals. Pool integration is owned
by root and must retain the original44 cells and incompatible-state outcomes.

All numerical work, including saved preflight, runs on Slurm. OneCPU suffices
for preflight; molecular stages use the pinned actual MACE-hybrid environment
and one allocated GPU. Native MACE uses the CPU grant; ALPB runs afterwards,
single-threaded. Root owns configurations, launches and shared status.
