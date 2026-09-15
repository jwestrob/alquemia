# Quantum endpoints complete; parser recovery — 2026-09-15

This updates the earlier queued status in AUDIT.md; old receipts remain intact.

- All six ORCA endpoints terminated normally. Job 1198934: 1491 seconds on
  64 allocated CPUs, **95,424 allocated core-seconds**, batch MaxRSS 10,318,208 KiB.
- Follow-on 1198939 exited without physical calculations: the original parser
  also matched the MBIS LARGEPRINT atomic dipole rows. A scheduler COMPLETED
  status did not imply a successful scientific result; all six failures were
  retained explicitly in the solver result.
- The corrected parser selects only the named ATOM/CHARGE/POPULATION/SPIN
  table ending at TOTAL. All six real outputs now pass atom order, completeness
  and the unchanged 1e-4 e charge-sum gate. No ECP shift or renormalization.
- Tests after the fix: **21 tests; 20 pass, one APBS identity integration
  skipped because that result is still absent**. Real MBIS integration now runs
  successfully. Additional corrupted-header and actual accounting tests pass.
- Original wavefunctions/densities are retained. Recovery job **1198958**
  performs only the already planned ESP/APBS stage, in `solver_recovery/`.
  No new high-level endpoint evaluations. Its admission budget subtracts the
  full 688 allocated core-seconds used by 1198939, leaving **52,288**. A campaign
  lock and admission ledger require terminal accounting for subsequent recovery.
- Scientific settings, inputs, tolerances and physical checks are unchanged.
  There is still no validated environmental correction or production-cost claim.

Current commands:

```bash
squeue -j 1198958
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/affordable_compare.py \
  --root /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs \
  --solver-result workspaces/affordable_challenger_20260915/solver_recovery/solver_result.json \
  --output diagnostics/affordable_challenger_20260915/comparison_after_recovery
```

Run the comparison after recovery has produced its result. Its own job writes
`solver_recovery/REPORT.md` automatically. The task-owned watcher records
`solver_recovery_terminal_receipt.json`; submission metadata records its PID.
Do not rerun successful quantum endpoints to address a parser problem.
