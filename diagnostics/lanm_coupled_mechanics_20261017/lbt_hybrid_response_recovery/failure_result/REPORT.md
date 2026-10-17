# Certified hybrid-response recovery: closed after a second native failure

The 7CCN/La search remains **unavailable**. Job13048 correctly resumed the original checkpoint and accepted the certified third trial, but embedded GFN2 failed again at cumulative evaluation11. The fixed second-failure stop rule closes this finite recovery. No paired four-state pool is available, and no partial endpoint is substituted.

## Actual outcome

- 10 accepted steps /11 distinct charged evaluations; original80/200 limits retained.
- Failed new point008: embedded GFN2 reached500 SCC iterations with last printed energy change−15.0112Eh and RMS charge residual0.742e. This is large unresolved oscillation, not a near-threshold miss.
- The same trial's isolated GFN2 converged in84 cycles; MACE/MM completed. ALPB was not called for the failed trial.
- Original default solver, parameters, geometry/state,300K electronic temperature and numerical gates remained unchanged after the one explicitly certified replacement.
- The last complete point has work−15.295612kcal/mol and maximum/RMS forces5.702562/1.574863kcal/mol/Å. These are diagnostic partial values, not an admitted search endpoint or a minimum.

This demonstrates inadequate numerical reliability of this embedded-GFN response route for the declared four-state search. It does not establish an erroneous experimental structure, the true electronic ground state, or which mechanical Hamiltonian is accurate. The original successful four-call qualification remains valid only for its tested coordinates.

## Cost and disposition

Actual151 scheduler seconds on12CPUs and oneH200:1,812 allocated CPU-seconds and151 allocated GPU-seconds. Reported TotalCPU177.249s and MaxRSS1,386,964KiB are retained separately; utilization is not inferred. Application149.124s. New requests:8MACE,16GFN2,8MM configurations/32 group queries and8ALPB (one certified-point completion plus seven successful new evaluations). Earlier failure, qualification and preflight costs remain separate.

Both original and recovery failures are immutable. Unexecuted paired-pool/comparison bridge drafts are gated closed. The independent whole256-atom PBE0 force/cost scout13051 and Mex derivative scouts12977/12978 continue under their own gates. No additional optimizer or SCC repair is launched.

Exact raw paths/hashes and actual scheduler accounting are in ROOT_REVIEW.json and SACCT.txt. Copied RESULT/SEARCH/failed-component receipts preserve the complete unsuccessful result; this report is a proactive inspection, not a delivered-event acknowledgement.
