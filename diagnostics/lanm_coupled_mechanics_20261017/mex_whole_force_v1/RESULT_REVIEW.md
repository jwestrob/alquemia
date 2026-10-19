# Independent actual-result review

19 October 2026. **PASS for the declared BEFORE directional consistency check.**
No admission or interpretation blocker was found. This review inspected saved
fields, source, native receipts, exact identities and file hashes only; it did
not recalculate derivatives, import scientific code or run molecular work.

## Actual evidence

Preflight 14810 is `preflight_passed` with zero molecular calls. Run 14811 is
`passed`, with `all_declared_passed: true`. Its binding and manifest exactly
match preflight. The run configuration differs only in the byte-identical
implementation copy's path and the actual preflight-result binding.

All 32 declared state/direction/step entries are present and passed, including
all 96 separate MACE, ALPB and total flags. Both original steps, 0.0002 and
0.0001 Å, retain the gate
`abs(FD - analytic) <= 0.005 + 0.0001 * max(abs(FD), abs(analytic))`
in kcal/mol/Å. There are no failed, omitted or unavailable slots.

All 64 raw MACE receipts match their requests: 2049 atoms, charge +2, one graph,
target index 2045, the three Nd spectators, and the declared La/Dy multiplicities
10/15 and 4/9. Each contains the full force array and declares the negative
derivative of the same vacuum energy. All 16 new Dy ALPB executions returned
zero; their commands, raw-input pins and physical mappings match the admitted
binary and receipts. The run retains four saved MACE centers, two saved ALPB
centers and 16 saved La solvent displacements. Across the inspected result,
configuration and 64 requests, 498 distinct referenced file/hash pairs verified.

Recorded calls are exactly 64 MACE forwards and 16 ALPB evaluations, with one
model load and zero central, MM, GFN, DFT or optimization calls. The unchanged
source admits completed native slots only after the 600-second forward and
90%-device-memory gates. Runtime reports NVIDIA H200, 12 CPUs, torch 2.8.0,
mace-torch 0.3.16 and ASE 3.29.0; elapsed adapter time is 78.995696 seconds.
The archived scheduler receipt records 81 seconds for 14811 and 8 seconds for
preflight: total allocated cost 980 CPU-seconds and 81 GPU-seconds.

## Limits and decision

Close this finite check as passed. It supports energy/gradient consistency only
for the four sampled directions at the saved BEFORE geometry under each stated
aggregate input. It does not qualify the AFTER geometry, other directions,
physical accuracy, localized spins, an electronic ground state, optimization or
natural-source transfer. Both aggregate hypotheses and their prior energy
sensitivity remain unresolved. No accuracy inference or new molecular run
follows; the independent 7CCN reference remains a separate question.

## Exact pins

| Artifact | SHA256 |
| --- | --- |
| Actual run RESULT | `c8c6166452996b862da4dc9fb35a758bc163ff44a60d60ceeba79051f27c4c07` |
| Run CONFIG | `94718e39aedd697e067764a3f0353137d65abb1130f59ac8998c4754b1988d00` |
| QUALIFY.py | `bf32bc2322bbfe20deebeba1a6c7de5783483d934864db9decd365e7184438b0` |
| Shared MANIFEST | `27c91c6c883a3835af1ff2d210d972d4bec65a8435f8e96e67ac4ff51e7171cc` |
| Actual preflight RESULT | `9a66cf943451d227f86e6ced46eedfc14f41aae602c4495530104fe0eeab3c88` |
| SPEC | `5a025950c6527c062984d6fc7c56c63618736566fc010fdfbf407ae630148b61` |
| Original independent source review | `67adefadec32a31f66a10c3396fc2ec519c02c1babd8efa860c8384dd0058a7a` |
| Archived COST | `86e5eed6322ce58bb02b491be358dbd9edef654300cc55b34dab04faf747fb26` |

Readable native results and requests remain under
`workspaces/lanm_coupled_mechanics_20261017/mex_whole_force_run_v1/execution/`.
The [result pin](actual_run_v1/RESULT_PIN.json) and
[accounting](actual_run_v1/COST.json) identify the archived evidence.
