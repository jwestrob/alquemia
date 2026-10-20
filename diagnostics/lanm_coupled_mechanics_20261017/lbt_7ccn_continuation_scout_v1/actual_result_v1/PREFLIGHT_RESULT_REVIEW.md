# Independent actual preflight review

19 October 2026. **PASS.** Actual Slurm preflight 15192 satisfies the frozen
saved-input gates. The separately bound run decision and configuration preserve
the declared single 7CCN-La continuation. No admission blocker found.

Review was limited to saved fields, literal JSON identities and file hashes.
No scientific imports, numerical replay, source edits or submissions were made.

## Evidence

- RESULT is `preflight_passed` for the reviewed implementation and configuration,
  with the exact expected scientific binding. It ran on the qualified
  `node-224-2t-8gpu-1`, one CPU, no visible GPU. Recorded application elapsed
  time is 10.559697 seconds; this is separate from the zero native-call costs.
- Original initial point `50a597c8`, terminal point `cf59ef7a`, checkpoint
  `4e3a6123` and prepared source `4d064ff4` are preserved. Counters remain
  80 accepted steps/80 distinct evaluations, with eight history pairs and trust
  0.1 Å. Settings match the original checkpoint, including force convergence
  0.10 maximum and 0.03 atom-vector RMS. New cumulative limits are 92/120.
- All 99 cache rows match the original point-record list, with 97 unique keys.
  Each record and its selected cache authority hash-verify, share the saved key,
  and retain complete 7CCN-La identity. Initial/terminal replay and the original
  derivative/scout admissions passed in the executed frozen adapter.
- Both saved physical guards pass with compatible state, all 256 atoms movable,
  and no bond, chirality, peptide-basin or severe-overlap failures. The terminal
  guard exactly equals the original search's final guard. The source was not
  rebased; the physical recipe and native package versions match the terminal.
- Requests are empty. Every native cost counter and new-component counter is
  zero; `model_initialized` is false and `no_model_initialized` is true.
  No other metal or search ran. Physical-state and force-accuracy qualification
  remain false. In total, 117 referenced file/hash pairs were independently
  verified during this review.

## Run binding and limits

CONFIG_RUN differs from CONFIG_PREFLIGHT only in `preflight_result` and
`run_decision`. Both pins verify. The decision binds this actual preflight,
the identical scientific binding, original checkpoint, state `[7CCN, La]`
and cumulative 92/120 limits. This admits only the declared additional
12 accepted steps/40 new distinct evaluations, with unchanged native settings
and stop-on-failure behavior. The explicit authorized QoS override is a
scheduling change; it does not alter these scientific caps.

This pass establishes saved-input readiness, not the outcome of molecular
continuation, finite-step reference descent, force accuracy or convergence.
The actual run must preserve partial/failure evidence and undergo its own review.

| Artifact | SHA256 |
| --- | --- |
| Actual preflight RESULT | `2dd9db660a7907077872b35a50b71e571fa2782ab5a12fa50200bdc3deadd245` |
| CONTINUE.py | `0c14a03e0c3d13884739d8be33a6c2f2c39f734c5374a36de0f4d511d7573b8e` |
| CONFIG_PREFLIGHT.json | `ad0d09964c043e89840a84a7fd7d8fd8b72b20e0805b8c35eca691eebbb184f6` |
| CONFIG_RUN.json | `7127a7fcba909e71a0ec3675a89159eecf627e23ce61def59e8aebc6935a14f7` |
| ROOT_RUN_DECISION.json | `a21d2c53f7dc2f5c6598f12462ce0b468370d3a9a1b1e9e1d6edfecf6a6fae32` |

Original result:
`workspaces/lanm_coupled_mechanics_20261017/lbt_7ccn_continuation_scout_v1/preflight/execution/RESULT.json`.
