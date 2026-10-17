# Downstream bridge drafts: unexecuted, gate closed

Recovery job **13048** is unavailable after **10 accepted steps and 11 distinct
evaluations**. Its saved SEARCH reports `physical_state_admitted: false` and
`Composite unavailable: ValueError('Native GFN component unavailable')`.
Under the frozen stop-on-another-backend-failure rule, this response branch is
closed. The original failed search and the failed recovery remain separate.

Saved evidence:

- [Recovery RESULT](../../../workspaces/lanm_coupled_mechanics_20261017/lbt_hybrid_response_recovery_v1/execution/RESULT.json), SHA256 `f1c4f9c2f359528f1f934fac7d0c29d0e355d695405349fe6eb51efa59fe3b73`.
- [Recovery SEARCH](../../../workspaces/lanm_coupled_mechanics_20261017/lbt_hybrid_response_recovery_v1/execution/SEARCH.json), SHA256 `ffbb32508d7e39ab5033847c6aab0150ca5d30eec98771a840b485a73f19748c`.

The following source-only drafts are preserved unchanged and **were not
executed**:

- `PAIRED_POOL.py`, SHA256 `91487e1390e81f12aa6651d0269a874adbe555dd05981f81a3b37a2ba6ae8d6c`.
- `COMPARE.py`, SHA256 `cebc2af1fd4bd3404425f42dab0cb3a11f55ebb3ba794def2dd4db6403c96180`.
- `POOL_SPEC.md`, SHA256 `d8c3b61b6b2cf6e6e954568a7d729c98460443ec1c3cc102a89915c927834c16`.
- `POOL_CONFIG.template.json`, SHA256 `412e8e22d6dc26826f037e0173255c87b86c1f58bffad65df3da1a4276bba2eb`.
- `COMPARE_CONFIG.template.json`, SHA256 `ba54a22e1545b4dd207dfa503553e8d05de1e5a1f83433a5c424667ed0e34212`.

Their complete four-state admission requirement is unmet; future result fields
remain unbound. AST and source-pin checks passed, but independent bridge review
was deferred when the branch closed. These drafts are not evidence of a paired
pool or completed comparison. No partial-survivor pool, further SCC repair,
optimizer restart, renderer or additional implementation follows this result.
Root owns the raw-failure diagnosis and shared status.
