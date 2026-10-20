# Independent native endpoint review

**15545 meets the unchanged native admission rules. Its admitted reference forces do not show stationarity at the model-converged geometry.** Native SCF convergence is not an all-residual, force-accuracy or physical-state qualification.

## Exact evidence and admission

| Artifact | SHA256 |
| --- | --- |
| Run RESULT | `a53125b98e597ea0473e06d58316ac2b32c914e9647da595f163cf0aa76a8c19` |
| Executed REFERENCE.py | `f251ca1e9514fb0c76659faf603b03ce2daf3684204fe4e98ef44e5d448bb240` |
| Run configuration | `e995bda312a070725daa866f9964c0df42d2bb4cd865f4ba787a44908ffe5f92` |
| Native execution receipt | `ed7f078a89b45d581729343350415384f5f1be763ec4a7ed1220e0e2725ae5f1` |
| Native endpoint.out | `6cddbc036f4a10ff99ff27e42ef3c289f4ba2f252ce42c36102d6f55495fbf29` |
| Native endpoint.engrad | `67e125366f884f659ff208fe1b4060e963813c248b718e658c7606eec72fcbfb` |
| New force report | `e6c0d039500c9deb87a53f31f974dfdc741c1cd801e481b51dc3316d100b7a7c` |

The workspace is `workspaces/lanm_coupled_mechanics_20261017/lbt_7ccn_converged_reference_v1/run/execution/`. I checked 36 selected artifact pins plus the original POINT, search, final checkpoint and source pins. Scientific binding, input evidence and both solvent admissions exactly match actual preflight 15543. The old native row is an exact reuse of 15321, including its original `new_La` task identity. The new endpoint binds POINT `6605f43c`, checkpoint `0461e30d` and the unchanged source; no rebasing occurred.

Actual output and receipt agree on ORCA 6.1.1 normal termination, return code 0, native SCF convergence and complete analytic engrad. The energy is `-6765.682831403961 Eh` in both output and engrad. Both admitted rows and both force tables contain all 256 atoms in the declared physical-ID order. One new DFT call and one runner invocation occurred; every other molecular, geometry and optimizer count is zero.

The input template hash `b94c22fc` is identical to 15321: PBE0-D4, def2-TZVP/def2-J, RIJCOSX, DefGrid3, VeryTightSCF, EnGrad, NoAutostart and PModel, with the same La basis/ECP46. Molecular charge is −2, physical/effective multiplicity 1, with 1030 explicit electrons, no external charges or caps. The printed atomic-guess multiplicity 2 does not replace the molecular singlet. Execution used 12 MPI ranks, one thread/rank and allocation-derived MaxCore 3100 MB from 48000 MiB with 25% reserve.

## Actual SCF stopping evidence

The output explicitly reports **Gradient check signals convergence**, after 18 cycles. MaxIter 125, ConvCheckMode Total+1el-Energy and ConvForced 0 are unchanged.

| Printed residual | Value | Threshold | Within threshold |
| --- | ---: | ---: | --- |
| Energy change | 7.1741e-9 | 1e-9 | No |
| Maximum density change | 5.1354e-5 | 1e-8 | No |
| RMS density change | 3.0492e-7 | 1e-9 | No |
| DIIS error | 1.2471e-3 | 1e-8 | No |
| Orbital gradient | 2.1020e-7 | 2e-6 | Yes |
| Orbital rotation | 5.8796e-7 | 2e-6 | Yes |

`all_printed_residuals_pass` remains false. The reused old endpoint stopped with a passing energy-change criterion but a failing orbital-gradient criterion; the two stopping routes must not be conflated. The early small-gap warning −0.061 is preserved, without treating it as a measured final gap or stability proof.

Final gridx3 recomputation completed and printed an exchange correction of −0.000002320 Eh. This scalar is not a nuclear-force error estimate. The earlier failed Mex external-point-charge derivative is a distinct case and cannot be transferred mechanically to this nuclear-gradient calculation.

## Reported force stationarity

The saved step 725 force report gives PBE0+ALPB maximum/RMS forces **7.868978/3.587757 kcal/mol/Å**, compared with model **0.075055/0.029611**. Model-minus-reference force-error RMS is **3.589135**. At step 92 the reference RMS was 3.639505: the admitted reference force magnitude remains substantial despite model convergence. The final normalized step 92→725 chord projection is −1.754370 for reference+ALPB and −0.030674 for model+ALPB; that chord is not a new optimizer direction.

This establishes disagreement under the admitted analytic recipe, not a calibrated uncertainty bound or proof of which model is physically correct. Independent derivative qualification remains absent. Conditional free termini, the vacuum dianion, post-added solvent without density polarization, ECP/no-SOC assumptions and electronic-state uncertainty remain. Both `physical_state_qualified` and `force_consistency_qualified` stay false.

The inspected scheduler event records COMPLETED, 25985 seconds/12 CPUs, exit 0:0. Native elapsed time is 25975.357699 seconds; application timing is not CPU utilization. Root owns final accounting and acknowledgement. Native evidence supports the declared saved collector, not a repeated endpoint or automatic optimization.

Review used source/output text, hashes, metadata equality and existing reported force metrics only. No gradients or energy differences were recomputed, no scientific package imported, and no jobs or source changes made.
