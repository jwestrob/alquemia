# Finite Stage2 saved collector

The question is numerical energy/gradient/virial consistency of the frozen old e41f La response model. Sources, 32 prepared displacements, eight fixed directions, three bypass decks and the original numerical contract are pinned in AUTHORITIES.json. This package performs no native calls, FFT/operator replay, assembly, factorization, optimization or parameter selection.

## Exact modes

`central_scout` requires exactly four records: active central plus old-library/absent, clone/absent and clone/zero bypasses. It checks the actual Stage1 central state/certificate, raw solver termination and counters, both common h additions, stationary unequal-dual polarization energy and exact complete E/G/virial numeric tokens across bypasses. Success is `stage2_central_bypass_scout_passed`; zero FD checks are claimed.

`full` requires exactly 36 gradient records, all 32 separate displaced certificates/initializers and the inherited central certificate. Every active gradient must bind its own expected and actual initialized-state hash. All 24 derivative checks are retained: two FD/analytic comparisons and one two-step comparison per fixed direction. Missing, failed or unadmitted records remain explicit and cannot produce a partial pass. Success is `stage2_finite_derivative_qualification_passed`. This is finite numerical consistency, not physical accuracy, full Cartesian qualification or trajectory stability.

## Arithmetic and readbacks

Native output is parsed once with exact ordered records and finite tokens. Decimal arithmetic with 70 digits computes literal finite differences, gradient or virial projections, print allowances and central dual contraction. Tangent binary64 values are used exactly via Decimal.from_float. FD gates stay 0.001 kcal/mol/angstrom plus literal token printing bounds. The two-step check retains the inherited sum of both slope allowances, including each analytic print term; that conservatively includes the common analytic uncertainty twice.

The central energy check is -electric/(2*dielec) sum u dot corrected_Ep; u and v are kept distinct. The original gate is 1e-8 kcal/mol plus 128*N*machine-epsilon times the maximum scale, N=35874. No physical threshold or extra induction application is introduced. Both reported and independently reconstructed shared-h errors use the inherited 128*N*epsilon bound. PCG's original pre-peek iteration/residual controls remain unchanged; post-peek u/v are reported, with no new residual solve. Old archive internal counters are unavailable and never filled with zero.

Three bypass E/G/virial numeric token sequences must match bytewise after whitespace separation; roles, counters and polarization-only records are outside this exact parity comparison. Every runtime key must equal its exact prepared key plus the declared read-only suffix `\nPME-NORM-EXPORT 1\n`; xyz and parameters remain exact. Source_inputs and actual runtime inputs are both bound.

## Execution interface

Planck's NATIVE.py produces `stage2_gradient_complete` records with implementation, state_id, role, gradient_output, input_binding, configuration, release, build, source_manifest, source_inputs, inputs, certificate, and commands.native.returncode. Active records additionally carry expected_state and actual_state pins. Counts must equal native_programs=1, gradient_calls=1, energy_only_calls=0, extra_operator_applications=0. Saved runtime and cost records are retained.

CERTIFY.py produces `stage2_state_nominal_certificate_passed` for each changed state; PREPARE.py produces the legacy central certificate metadata. Both bind their implementation, candidate, manifest, source_manifest, state_id, exact inputs, expected_state, pme_export and sufficient_numerical_test_passed. Changed certificates bind build plus initializer/assembly/factor receipts; initializer carries native_programs=1. Central metadata binds actual Stage1 expected state/PME and job20053 root factor admission, with no new initializer. Nominal U_B/U_R/delta_odd/U_nominal/tau are retained; e_app remains null. Assemblies/factors are not repeated.

Populate a fresh copy of the appropriate CONFIG_*_UNBOUND.json after terminal evidence exists; bind actual producers/build/execution manifest/root release and each result/certificate. Change status to `bound_stage2_saved_collection`. Run COLLECT.py CONFIG_RUN.json fresh-output-directory in one thread. It writes immutable RESULT.json and REPORT.txt and exits nonzero on incomplete/failing qualification. No automatic retries follow. Scheduler accounting may be supplied as an exact saved pin; task costs and scheduler totals remain separate to avoid double counting.

## Source controls and pending evidence

Twenty-five manufactured controls cover strict record parsing, unequal duals, changed h/energy, missing or duplicate atom rows, original convergence, nonfinite/trailing records, derivative sign and strain virial, both steps, exact36/24 denominator, missing-cell failure and bypass parity. They use one synthetic site and synthetic panels only. No actual Stage2 output exists or was read. Runtime producer/build/release and the 32 certificates remain unbound. Root owns source review and all execution.
