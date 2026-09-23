# Activated native restart removes the Q88 numerical discontinuity

**A confirmed continuation of the native electronic state removes the large
Q88JH5 energy jump.** At both nearly identical geometries, self and cross starts
reach the same lower vacuum energy within 0.000093 kcal/mol. The original
4.804994 kcal/mol difference between La solvent-transfer descriptors shrinks
below 0.000090 kcal/mol. This is a useful numerical remedy on one consumed case;
no new La/Ca classification, reference, or production change has been made.

## What ran and why the activation recovery was necessary

Jacob authorized old/new Q88JH5 La geometry ×vacuum/ALPB ×self/opposite native
seed: eight logical tests. The initial attempt copied only the saved xtbw and
removed NoAutostart. All eight calls completed as **SAD initial guesses**, with
zero confirmed restarts. Their observed energies exactly reproduce the original
values; their qualified restart fields remain null. See
[the preserved first-attempt result](SEED_ONLY_RESULT.md).

The documented AutoStart trigger is an existing same-basename GBW; native xtbw
takes precedence when restarting. Root approved adding each seed's matching,
already archived GBW as a technical activation correction. The recovery retained
the original Hamiltonian, source coordinates, charge−2/singlet, 438 native
electrons, 300K, native mixer, MaxIter500, eight ranks and exact solvent setting.
No extra input keyword or ordinary SCF was introduced. Both files came from the
same pinned source calculation. [Recovery agreement](RECOVERY_PLAN.md).

Recovery job1210185 **confirms all eight native restarts**: every output reports
`Checking for AutoStart` and `INITIAL GUESS: XTBRESTART`, followed by the special
xTB mixer. Source-matched parameters, electronic state, charge closure, normal
termination and executor receipts all pass. Immutable original GBW/xtbw, their
pre-run hashes, rewritten files, retained GES, full iterations, occupations and
charges are saved. This demonstrates actual native activation rather than an
assumption based on the input keyword or the presence of seed files.

## Full energy matrix

Unrounded parsed energies are in the immutable collection. Values below are Eh;
“opposite” always means the other geometry's seed in the same medium.

| Destination | Medium | Original SAD | Restart: self | Restart: opposite |
|---|---|---:|---:|---:|
| old adaptive | vacuum | −252.116197615034 | −252.116198367217 | −252.116198514926 |
| new template | vacuum | −252.108540025437 | −252.116198513619 | −252.116198392367 |
| old adaptive | ALPB | −252.570227444230 | −252.570236658914 | −252.570236665730 |
| new template | ALPB | −252.570227099694 | −252.570236692564 | −252.570236685757 |

The new-geometry vacuum self continuation lowers the energy by **4.805774 kcal/mol**;
the old-geometry self continuation changes it by −0.000472 kcal/mol. ALPB
continuations change by about −0.0058 to−0.0060 kcal/mol. The full four-initialization
ranges are 0.000092689 kcal/mol in vacuum and 0.000021116 kcal/mol in ALPB.

| Destination | Original ALPB−vacuum | Restart: self | Restart: opposite |
|---|---:|---:|---:|
| old adaptive | −284.9080192991 | −284.9133295986 | −284.9132411869 |
| new template | −289.7130131196 | −284.9132588457 | −284.9133306610 |

Transfers are kcal/mol, converted once with the established factor627.509474.
The new-minus-old transfer differences are +0.000070753 for self starts and
−0.000089474 for cross starts. These comfortably pass the predeclared0.1 endpoint
and0.2 differential-transfer diagnostic tolerances. No preferred branch was
selected or used to overwrite an earlier score.

## Interpretation and remaining limits

The original high-energy result is not stable to its own native continuation.
Self continuation alone reaches the lower solution; cross initialization agrees.
This supports premature stopping or insufficient settling of the original SAD
iteration, rather than a persistent geometry-driven energetic split. It does
not uniquely identify the internal stopping defect or prove the ground state.

All restarted printed maximum/RMS density residuals still exceed their generic
displayed tolerances: maximum density0.000363–0.001263 versus1e−5, RMS density
1.48e−5–2.57e−5 versus1e−6. These are recorded diagnostics, not proof of the
native internal charge/multipole convergence criterion. The observed agreement
is an energy-continuity result for this case, not universal solver qualification.

Only La was recomputed. No full Ca−La score, biological call, new reference,
DFT-agreement claim or broad accuracy gain follows from this diagnostic. No MACE,
DFT, standalone xTB, geometry search or additional state was evaluated. Production
and historical outputs are unchanged.

## Actual cost and tests

| Job | Scientific status | Elapsed | CPU allocation | Core-seconds |
|---|---|---:|---:|---:|
|1210179|8 completed SAD calls;0 native restarts|26s|64|1664|
|1210185|8 completed, confirmed native restarts|27s|64|1728|
|Total|16 actual attempts /8 logical tests|—|—|**3392**|

Both used128GiB host memory on node-224-2t-8gpu-1 with eight8-rank workers.
No GPU was requested. Inner worker allocation totals2588.326568core-seconds;
unmetered local preparation/reporting is additional. Scheduler peak memory was
unavailable, not zero. Seven focused real-artifact tests pass with no skips,
covering exact recipes/seed pairing, corrupted-source rejection, partial-status
handling, actual activation markers and transfer algebra. Molecular execution
and parser/fixture tests are distinct checks.

Actual result:
`workspaces/native_xtb_restart_20260922/recovery_v1/collection_1210185.json`.
`SUMMARY.json` retains raw changes and continuity; both run directories retain
their own `COSTS.json`, scheduler and endpoint receipts. Exact pins are in
[ARTIFACTS.json](ARTIFACTS.json); report-only commands in [COMMANDS.md](COMMANDS.md).

## Minimal predictive follow-up — proposed only

First test whether the same self-continuation changes useful discrimination,
without another geometry search: freeze the existing five-candidate pools for
1H4I,4MAE,Q88JH5 and the still-wrong A0A3F2YLY8 Ca-sample1. Restart both metals
in both media at every existing candidate, at most80 calls, conditional on
retained compatible seed pairs. Reuse exact native MACE energies, rescore each
complete fixed pool, retain all old components and compare raw contrasts plus
explicit frozen-band transfer. These are consumed development cases. Do not
recalibrate four cases or inherit validated decisions for the changed numerical
policy. A declared repeat-continuation check would still be needed to qualify
the stopping policy before any wider reference or production change.

This is a proposed next scope, not a submitted job or an automatic retry.
