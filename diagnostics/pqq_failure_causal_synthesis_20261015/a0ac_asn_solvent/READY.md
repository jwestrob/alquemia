# Asn290 solvent stage submitted

Root submitted job **4623** on15October2026 at02:07:29UTC, using the existing
standard-partition seven-CPU wrapper. Completion watcher PID3466505 is armed;
the exact submission and wake commands are preserved in
`workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_solvent_v1/SUBMISSION.json`.
Do not resubmit this active stage.

## Actual admission and current status

MACE4617 completed all20 new cells and retained all ten sources. The declared
sign-blind1kcal/mol mathematical common3→4 gate passed only for Ca-conditioned
sample4, with recorded change+1.460981253638387kcal/mol. Its La row selects the
flip; the other nine common4 outcomes remain unchanged at this vacuum stage.
This is not a calibrated class result.

Actual4623 `execution/PREFLIGHT.json` is **passed**. It revalidated all60 existing
solvent cells and reproduced390 historical pool comparisons with maximum
residual0.0kcal/mol against the existing1e−8 tolerance. The actual MACE matrix
and stage gate replay also passed. No terminal execution receipt was present at
this readiness update; solvent results remain pending.

Bound configuration SHA256:
`0c22b4056a64b4d0dabb3922eff20b72e3a3bd6455dbd740c7501a484415010c`.

Reviewed source hashes:

- `SOLVENT.py`: `ada49e46ebe58530c47bd28489332a60b0d78ea1806b46ac902995c347792365`.
- `RUN.sbatch`: `6c954d32079d40d87d9caf7e6e13c0a1b85e8e9a05f202cfb86886452912a5d4`.
- `PLAN.md`: `262a19de6f9d60678a8fafc595087004abcaa1b9b9e877fd842026aaf74c72e3`.

Root and independent review cleared the frozen adapter. Local checks were only
AST/Bash syntax and metadata identity checks. The actual scientific preflight
ran in Slurm4623.

## Exact stage commands, recorded for reproducibility

Root ran this metadata binding after actual MACE completion and gate review:

```bash
python diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_solvent/SOLVENT.py prepare \
  --score-config workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_scoring_v1/CONFIG.json \
  --score-collection workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_scoring_v1/execution/FINAL_COLLECTION.json \
  --output workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_solvent_v1
```

The following is the submitted command, expressed from the repository root;
**it is a record, not a retry instruction**:

```bash
asn_solvent_dir="$PWD/workspaces/pqq_failure_causal_synthesis_20261015/asn_rotamer_solvent_v1"
sbatch --parsable \
  --output="$asn_solvent_dir/slurm_%j.out" \
  --error="$asn_solvent_dir/slurm_%j.err" \
  "$asn_solvent_dir/RUN.sbatch" "$asn_solvent_dir/CONFIG.json" \
  0c22b4056a64b4d0dabb3922eff20b72e3a3bd6455dbd740c7501a484415010c execute
```

## Preserved scope and outputs

The stage evaluates all20 fixed flips with the unchanged standalone GFN2 and
native vacuum-charge ALPB recipe, reusing60 original common3 cells. Maximum80
SCC calls and40 native replays; zero new MACE, DFT or optimizer calls. The focused
dispatcher copy only corrects the protein metadata label. Actual native and
component checks, precision limits, parameter lineage and all four solvent terms
remain unchanged. Historical absolute-v1 failures remain visible.

Final outputs will be `execution/FINAL_COLLECTION.json`, `PARTIAL_RESULTS.json`
and `EXECUTION.json`, with original3 and complete common4 shown separately under
operational, mathematical and origin policies. All ten sources remain visible;
failed required cells cannot produce a partial-pool score. The fixed short H
contacts are represented in the204-atom inputs and have not been repaired or
relaxed. This experiment does not test every Asn conformation or introduce class
bands, recalibration, production changes or automatic follow-on work.
