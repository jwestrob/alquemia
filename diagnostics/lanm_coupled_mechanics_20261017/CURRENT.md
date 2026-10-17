# Coupled La/Dy mechanics is active

17 October 2026 UTC. Jacob activated “Implement the saved La/Dy mechanics plan. Work autonomously, validate the physics, and document the results.” PLAN.md is the finite protocol. Root owns all submissions and integration; PQQ production remains unchanged.

## Actual results

- Mex qualification12952 passes all34 composite evaluations and all component/total derivative checks, including backbone, boundary, water and spectator movement. Search preflight12955 passes.
- LBT partition12953 passes all four preparations, with native MM XML unchanged from failed12950. Baseline172/expanded210 electronic atoms each include two mapped caps.
- LBT original preflight12954 failed only old exact-transform replay. Saved diagnostic12956 and reviewed V2 preserve the archived transform and all physical/scientific gates. Preflight12957 and all32 derivative cells12958 pass.
- Old LBT44-cell pool replay12959 passes with zero calls and an explicit required-missing-cell test.
- LBT scout12960 passes:12 accepted steps/12 new composites; energy−63.9031 kcal/mol and max force83.2803→12.7867 kcal/mol/Å. Physical/state checks pass, but it is not converged. Read LBT_QUALIFICATION_AND_SCOUT.md.

## Owned execution

Root submitted the six planned searches after reviewing actual scout admission and cost. LBT7CCO/La12961 completed/admitted80 steps, maxforce2.14091/RMS0.921311 and energy−98.5724 kcal/mol from the prior restricted terminal; still nonconverged. The other five initial jobs12962–12966 failed before molecular calls because their node lacked a pinned historical compiler path. Configuration-identical host-only recoveries use the qualified H200 node:

| System/state | Job |
| --- | --- |
| LBT7CCO La, exact scout continuation | 12961 |
| LBT7CCO Dy | 12967 (pre-call failure12962 retained) |
| LBT7CCN La | 12968 (pre-call failure12963 retained) |
| LBT7CCN Dy | 12969 (pre-call failure12964 retained) |
| Mex La | 12970 (pre-call failure12965 retained) |
| Mex Dy | 12971 (pre-call failure12966 retained) |

All four LBT searches 12961/12967/12968/12969 are complete and admitted at 80 steps, all nonconverged. Their model energy reductions are 98.57–137.65 kcal/mol and residual maximum forces are 2.14–10.17 kcal/mol/Å. Residuals and component work are pinned in BROAD_RESPONSE_CHECKPOINT.json. These are finite search endpoints, not minima or preference results. LBT paired scoring 12974 is queued, with at most four new cross cells.

Mex La 12970 is also complete and physically/state admitted: 80 steps, energy reduction 2,945.11 kcal/mol, residual maximum/RMS forces 11.05/1.42 kcal/mol/Å. Most lowering is in the MM ledger (−3,170.81 kcal/mol), partly opposed by solvent (+552.93); this is model-specific relaxation, not affinity. It took 716 scheduler seconds on 12 CPUs and one H200. Mex Dy 12971 is running. The full paired comparison must precede interpretation of metal-specific work. No automatic continuation beyond the frozen budget.

Zero-call hybrid preflight 12972 passes all eight saved component replays and freezes actual A2CA/A1N directions; nine-composite hybrid scout 12975 is queued. Mex common5 baseline replay 12973 passes all ten old cells and its missing-cell test. All jobs have completion wakes.

Embedded DFT preparation 12976 passes both reused origin-gradient checks and constructs the eight declared displaced tasks with zero DFT calls. The first external-charge La ± pair is queued as 12977/12978, 12 MPI slots each on the qualified CPU host. Completion wakes and manifest-aware, alert-only SCF health supervisors are recorded in each SUBMISSION.json. Read mex_dft/README.md for exact collector commands. Actual scout derivative admission and measured cost precede the remaining six endpoints.

All searches use 12 CPUs/one GPU and the frozen 80-step/200-evaluation limits. Scout work counts toward its continuation. Actual submission commands, pins and detached completion watchers are under each workspaces/lanm_coupled_mechanics_20261017/{lbt,mex}_run_*/SUBMISSION.json. Inspect actual terminal results before pool evaluation; never duplicate submissions.

## Ownership and next commands

- Root: source/config pins, scientific submissions, result integration, shared status.
- Leibniz: paired pools and remaining LBT derivative qualification are implemented and reviewed; no submissions.
- Harvey: final saved comparator and DFT qualification independently reviewed; no remaining static blocker.
- mechanics_mex_adapter: saved comparison staging is complete; preparing a saved-only Mex motion/MM diagnostic to select a defensible scaffold-stiffness segment from both actual metal responses. No new energies or shared-status edits.

Next: collect Mex Dy, then bind the two actual searches into the reviewed Mex paired pool (at most two new opposite-metal composites). Collect queued LBT paired scoring and hybrid scout. After actual hybrid scout admission, run its eight-origin matrix with exact reuse before the prepared expanded-region derivative check. Saved Module2 comparison source 062f1167 is independently reviewed and staged with five completed searches, all five host-startup failures and the original scout counted once; missing Dy/pools remain explicit until actual files arrive.

The saved comparison retains the prior conditional exchange for continuity, but labels whole-MACE LBT versus local-electronic/MM Mex as different scaffold Hamiltonians. It cannot qualify a matched-Hamiltonian preference. Module3 tests that distinction separately.

Current observed standard policy is12CPU/job120CPU/user3days, with standard_expand24CPU/job500CPU/user5days available. All numerical work uses Slurm, allocation-derived resources, unchanged parameter files and actual completion wakes. No priority changes, memory/time overrides, source-position springs or old clamped-core diagnosis. Goal remains active and unmet.

H200 partition separately limits one concurrent job/GPU per user; the account 120CPU/10jobs does not override this. CPU checks ran on qualified standard node-64-768g-1 after a recorded scheduler-only update. Keep jobs queued normally. Independent DFT pair 12977/12978 has been submitted; no derivative result exists yet.
