# Opt-in preparation commands

These operations use existing source protonation and reconstruct the fixed4.3Å
three-source pocket. They launch no molecular calculation. The older
`pqq_three_source.py` and released scoring commands are unchanged.

Run from the repository root. The following real07ab request creates a separate
reproduction directory; repeat invocation must use a new output directory because
scientific artifacts are immutable.

```bash
py=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
src=workspaces/accommodation_goal_20260920/plm_source_requests_v1/PQQSEQ_07ab500e3df76b30d71c_sources.json
ref=workspaces/motion_envelope_20260923/REFERENCE_PINNED_v1.json
out=workspaces/pqq_three_source_envelope_20260923/reproduction_v1
mkdir -p "$out"
"$py" scripts/pqq_three_source_envelope.py request \
  --sources "$src" --config "$src" --reference "$ref" \
  --protein-id PQQSEQ_07ab500e3df76b30d71c \
  --agreement diagnostics/pqq_three_source_envelope_20260923/PLAN.md \
  --output "$out/REQUEST.json"
"$py" scripts/pqq_three_source_envelope.py prepare \
  --request "$out/REQUEST.json" \
  --source-preparation workspaces/plm_fold_sensitivity_20260920/PQQSEQ_07ab500e3df76b30d71c/source_preparation/preparation.json \
  --output "$out/prepared"
"$py" scripts/pqq_three_source_envelope.py dry-run \
  --preparation "$out/prepared/PREPARATION.json"
"$py" scripts/pqq_three_source_envelope.py report \
  --preparation "$out/prepared/PREPARATION.json" --output "$out/REPORT.json"
```

Immediate inspection of the completed fixture, without reconstruction:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/pqq_three_source_envelope.py dry-run \
  --preparation workspaces/pqq_three_source_envelope_20260923/plm_v1/PQQSEQ_07ab500e3df76b30d71c/PREPARATION.json
```

The corresponding8344 source request, archived preparation and completed fixture
use the exact ID `PQQSEQ_83440678cbbd658047c9` in the same directories. Both
three-member source lists were declared previously; no favorable source choice
or biological classification is implied.

Focused tests:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  -m unittest discover -s tests -p test_pqq_three_source_envelope.py -v
```

## Execution handoff

`plm_v1/FINITE_TASKS.json` enumerates the six known paired origins,24 actual
strict scalar origin inputs, and12 pending source-connected proposals. It is a
preparation manifest, not an executor-compatible manifest. No execute command is
advertised or launched here. The fixed34 envelope frontends remain intentionally
restricted to their study; a later small execution adapter must call the reusable
`motion_envelope_run.build_tasks` after fresh paired origin forces, existing warm
MACE/optimizer kernels, and the qualified strict scalar runner. It must preserve
the exact three-source group identity, reference, original/accommodated contrasts,
complete common three-candidate matrix and strict three-member median/spread.

Prospective total:12 origin-forceMACE calls,12 bounded searches, at most12
cross-MACE calls and72 nativeGFN2 scalar calls. Search evaluation counts are data
dependent. WarmMACE uses the established oneH200/32CPU/200000MiB allocation;
scalar work uses rank1 with up to32 workers on32CPU/64GiB. Molecular cost is
unmeasured for these inputs. Root holds molecular integration until its separate
scientific decision; there is no automatic launch or default promotion.
