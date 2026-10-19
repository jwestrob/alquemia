# Independent static review

19 October 2026. **PASS for the declared zero-call Slurm preflight.** No material
source/configuration blocker found. Molecular execution still requires the
actual passed preflight and its exact root review binding. This review performed
source, JSON-identity and hash inspection only; no scientific imports, numerical
replay, tests, submissions or edits to implementation files.

## Scope and admission

- All 13 direct CONFIG_PREFLIGHT pins verify. The state is exactly 7CCN/La,
  with 12 additional accepted steps and 40 new distinct evaluations: cumulative
  limits 92/120 from the unchanged 80/80 consumed counters. Future preflight and
  run-decision pins are null.
- The original optimizer remains byte-identical (`d93821da`). The isolated copy
  changes only the hash import and entry/budget admission. It rejects a missing
  resume state and accepts only the complete original checkpoint's canonical
  digest, its original 80/200 limits and the new 92/120 limits. I independently
  verified that digest against checkpoint `4e3a6123`. The consumed-count check
  handles the lower new evaluation cap without resetting or generally unlocking
  budgets. No later checkpoint is admitted by this copy.
- All numerical optimizer bodies and settings remain unchanged: full-system
  maximum force 0.10 and atom-vector RMS 0.03 kcal/mol/Å, rigid projection,
  L-BFGS history/fallback, trust updates, Armijo/backtracking, native-component
  energy differences and source-relative physical guards.

## Replay, calls and outputs

The adapter retains the original initial point, terminal resume point, eight
history pairs and trust 0.1 Å. It invokes the original loader and Evaluator,
seeds all 99 archived point records, and checks original request/checkpoint
joins. The actual 99-record/80-request schema, original source and original
scout-admission flag match. Cached initial/terminal gradients and components,
unchanged state/recipe/versions and original-source guards must replay on Slurm.

Preflight never calls `Evaluator.evaluate`; native initialization remains lazy.
It requires no model object, no requests and zero component costs. Run admission
requires exact scientific binding and preflight-evidence equality, plus the
root decision for this one state/checkpoint and 92/120 caps, before search.

The preserved evaluator counts each new distinct attempt before native work;
the preserved search checks the cumulative limit before evaluation and stops on
native failure. It rejects retrying a failed exact state. Cached states do not
spend new-evaluation budget. Final counters separately retain new MACE/ALPB
attempts and cumulative evaluations; failures and the one-state denominator
remain visible. There is no opposite-metal, DFT, GFN, MM or additional-search
execution path. A state-changing outcome cannot support a bound-state conclusion.

Both staged launchers match the declared host and original environment:
preflight one CPU/no GPU; run 12 CPUs/one H200; default scheduler memory/time,
no requeue. They supply the interpreter themselves, so READY argv is passed
starting with the script path.

## Limits

Actual zero-call replay remains to be demonstrated. Static approval does not
claim that numerical or native execution has passed. Audit 15176 establishes
nominal local descent only; neither it nor this bounded model continuation
establishes finite-step reference descent, improved accuracy or convergence.
Budget exhaustion remains incomplete relaxation. No additional continuation,
reference or tolerance change follows automatically.

| Artifact | SHA256 |
| --- | --- |
| CONTINUE.py | `0c14a03e0c3d13884739d8be33a6c2f2c39f734c5374a36de0f4d511d7573b8e` |
| CARTESIAN_OPTIMIZER_SCOUT.py | `3badac29106412f8b29e79e874678c330af620fe5e6d4cc548dd0c3b3f0abfbc` |
| SPEC.md | `b8effde9bb57e0e9e4bf3809024a3f30e127ec289b97961b8b2d268e68e6f541` |
| CONFIG_PREFLIGHT.json | `ad0d09964c043e89840a84a7fd7d8fd8b72b20e0805b8c35eca691eebbb184f6` |
| READY.json | `95a03d4afd77a854be06dffedc0b51b15d4f3abd3bbf75106e8be87e21b1dd38` |
| Canonical original checkpoint | `2be9f7cd5ca1da78aebec332268a30a26a350b659af4f7d567ad00cb22847436` |
