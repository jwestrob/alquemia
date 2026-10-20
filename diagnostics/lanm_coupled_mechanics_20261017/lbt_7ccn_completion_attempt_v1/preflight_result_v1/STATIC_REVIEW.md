# Independent static review — finite 7CCN–La completion attempt

**Verdict: pass for the declared zero-call Slurm preflight.** No material source, state, cache or scope blocker was found. Static review does not replace the actual preflight or authorize a run without its separately bound root decision.

## Frozen bundle

| Artifact | SHA-256 |
| --- | --- |
| `CONTINUE.py` | `551d964badc41141cfc53a85d704135ec74c21f185e12b6cc50d578dc47910d2` |
| `CARTESIAN_OPTIMIZER_COMPLETION.py` | `c5092f2480fbf393d6ca0d67733ebf2131580c922a62a67925dea8945eb5acb9` |
| `SPEC.md` | `443c02213579a9e76a8536a011ac5170fd10e2c3d08d51e9067b167820f7b96d` |
| `CONFIG_PREFLIGHT.json` | `664ffb95e45d297904364e888d12dac63f0895b848f79e9d64de39afef52c710` |
| `READY.json` | `5bb2f87a3669ea710a93d871399fd5eff9a0b6cfd89e76df48216023aac1ba39` |
| `SCIENCE_DECISION.json` | `d0138c99f7ded2bd6f943d9542295e0a55e4f32eb35b3d4c022da5c97ac4ebea` |

Review covered source text/AST, literal saved JSON joins, record counts and 132 selected file/hash pairs. No scientific modules were imported, numerical checks replayed, molecular calls made or jobs submitted. Only this review file was written.

## Optimizer and checkpoint

The optimizer is exactly the prior reviewed scout source except for entry messages, the required old92/120→new172/292 limit transition and the permitted checkpoint digest. All other source bytes are unchanged. This preserves the direction calculation, rigid projection, fallback, Armijo backtracking, trust/history updates, physical guards and full-force convergence criteria: maximum atom-vector force0.10 and RMS0.03 kcal/mol/Å.

Admission requires the complete exact step92 checkpoint `317c2b24`, canonical JSON SHA `981ec084b57fa86599125619b9ecce7613c883597b23e38d3521be1eaa5c0da2`. Its actual92/92 counts, eight history pairs, trust0.1Å, original initial point and original prepared source are preserved. The adapter passes the original initial point and this checkpoint to the unchanged search. No other checkpoint, fresh start, later retry or generic budget reset is admitted.

## Cache, calls and state

- Pinned predecessor `0c14a03e.prepare` receives the exact previous run configuration `7127a7fc`. It performs the original saved source, component-gradient, version, derivative-qualification and guard checks without model initialization. The adapter reuses that same evaluator; it does not construct a second evaluator or recreate its directory.
- The actual old and scout records form113 archived records with109 unique cache keys. All92 historical requests exactly equal original80 plus scout12, and every requested/evaluated point is represented. Each record is re-admitted by the original native-component helper. The merged cache retains every prior key and replaces only fresh saved-loader aliases with archived authorities.
- Historical requests stay in the exact resume checkpoint. The unchanged checkpoint writer appends new requests; the new evaluator's costs/requests start at zero. Model load is lazy. Exact cached geometries reuse admitted components. Native failures retain their evidence and stop the search under the existing behavior.
- Whole256 La state, charge−2/multiplicity1, model weights/runtime, fixed-charge solvent, original source guard and initial geometry remain unchanged. Saved initial, step80 and step92 component/gradient identities plus the original-source terminal guard must replay in the actual preflight. This review did not perform those numerical operations.
- The new allowance is at most80 accepted steps and200 distinct composite evaluations, with cumulative172/292 ceilings. The search checks the remaining evaluation budget before a new request. There is at most one MACE and one ALPB call per attempted distinct composite; no DFT, GFN, MM, alternate-metal/source, intermediate reference or extra search path is introduced.

## Preflight, run and limits

Future success bindings are null in the preflight configuration. Run requires the actual passed preflight, identical scientific binding and preflight evidence, and a separately pinned root decision for this exact checkpoint/state/limit set. Fresh output directories preserve prior artifacts. No stage submits another job.

The inspected launchers retain the original Python environment and qualified H200 host: one CPU/no GPU for preflight,12CPUs/one GPU for run, no memory/time override and no requeue. Literal QoS headers remain `standard`; root must supply the declared `--qos=standard_expand` override. The source enforces the qualified host/runtime and12CPU/GPU run allocation.

The specification accurately preserves original physical-trial rejection/backtracking and accepted state-change handling. Budget exhaustion is incomplete relaxation; a state-changing outcome remains excluded from bound-state conclusions. There is no automatic restart, expansion or reference follow-up.

Actual collector15456 is pinned rationale, not a newly fitted accuracy gate. Its support for the prior incremental work does not establish force accuracy or predict the next segment. This attempt tests convergence of the unchanged model only. Reference SCF/gradient, vacuum-dianion, frozen-f, fixed-solvent and conditional terminal-chemistry limitations remain; no affinity, preference or stable-state claim follows.
