# Three explicit source requests: prepare and dry-run

This interface requires no historical canonical member or experimental label.
It currently supports the frozen `protenix_generic_PQQ` normalization contract
(the existing contract also handles the supplied AF3 CIFs), with exactly three
La-conditioned structures and explicit selectors. It performs no scoring.

## Inspect a completed real PLM preparation

```bash
THREE_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
"$THREE_PY" scripts/pqq_three_source.py dry-run \
  --preparation workspaces/pqq_three_source_20260923/plm_v1/PQQSEQ_07ab500e3df76b30d71c/PREPARATION.json
```

The same operation supports the8344 preparation in the adjacent directory. Both
remain unlabeled predictions, not validated affinity calls. No endpoint energy
was calculated for this new preparation.

## Explicit request and exact source-preparation reuse

The actual existing source file below declares all three paths/selectors and
contains the source configuration. This example creates new metadata/prepared
artifacts only; it invokes no protonation or molecular calculation.

```bash
THREE_BASE="$PWD/workspaces/pqq_three_source_example_v1"
THREE_SOURCE="$PWD/workspaces/accommodation_goal_20260920/plm_source_requests_v1/PQQSEQ_07ab500e3df76b30d71c_sources.json"
"$THREE_PY" scripts/pqq_three_source.py request \
  --sources "$THREE_SOURCE" --config "$THREE_SOURCE" \
  --protein-id PQQSEQ_07ab500e3df76b30d71c \
  --reference workspaces/union_triple_pilot_20260923/REFERENCE_v1.json \
  --agreement diagnostics/pqq_three_source_20260923/PLM_ADDITION.md \
  --output "$THREE_BASE/REQUEST.json"
"$THREE_PY" scripts/pqq_three_source.py prepare \
  --request "$THREE_BASE/REQUEST.json" --source-mode reuse-exact \
  --source-preparation workspaces/plm_fold_sensitivity_20260920/PQQSEQ_07ab500e3df76b30d71c/source_preparation/preparation.json \
  --output "$THREE_BASE/prepared"
"$THREE_PY" "$THREE_BASE/prepared/implementation/pqq_three_source.py" dry-run \
  --preparation "$THREE_BASE/prepared/PREPARATION.json"
```

For another explicit request, provide a JSON object with `sources` (or the
existing `cases` list), exactly three descriptors. Each needs `case_id`,
`source_structure` (absolute path or SHA256 pin), `source_conditioning_metal`,
`raw_source_metal`, normalized `metal`, `pqq`, `roles`, `assembly` and
`normalization`. The real file above is a complete example. Optional historical
labels/provenance remain metadata; a canonical flag is neither required nor used.
The config file must contain the qualified configuration itself or its `config`
member. No source edit or hidden path default is needed.

`prepare --source-mode fresh` omits `--source-preparation` and uses the existing
source protonator/graph policy in a new output directory. That option was not
executed in this delivery; use the existing allocation policy for preparation.
It does not run molecular endpoint scoring. A failure remains explicit, with no
substitution of an archived state or a reduced-member union.

PREPARATION.json retains all three statuses, common state anchor, exact union,
source/cap/water mappings, paired coordinates/charges, frozen reference, immutable
implementation and deferred executor handoff. A complete triple exports six
q0 tasks; an incomplete triple exports none. Four-mode selection waits for real
paired origin forces. No generic scientific execution command is enabled yet.

Focused actual-fixture tests:

```bash
"$THREE_PY" -m unittest discover -s tests -p test_pqq_three_source.py -v
```
