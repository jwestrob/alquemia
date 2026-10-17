# Matched-response monomer comparison: complete

**The corrected search objective does not establish a robust preference.** Saved
comparison 12769 gives a global balanced exchange score of **+0.982619 kcal/mol**,
within the predeclared 1 kcal/mol interpretive resolution, while the two LBT
sources still differ by **24.120991 kcal/mol**.

The score is R=(Dy−La)_Mex−(Dy−La)_LBT. Positive R favors conditional
Mex-La plus LBT-Dy. The resolution is neither a fitted classifier threshold nor
an uncertainty interval.

| LBT source/menu | Previous R | Matched-response R | Change, kcal/mol |
|---|---:|---:|---:|
| 7CCO source-local union | +1.019378 | **+0.982619** | −0.036759 |
| 7CCN source-local union | +29.281571 | **+25.103610** | −4.177962 |
| Global union | +1.019378 | **+0.982619** | −0.036759 |

The same completed Mex common5 is used throughout. All 15 earlier comparison
columns remain unchanged in the result. Each new LBT source-local menu includes
its two matched-objective terminals; the global menu selects the new 7CCO La and
Dy terminals (`geometry_18` and `geometry_19`). Selecting one source for both
metals does not remove the source-sensitivity limitation.

Matching LBT's search objective to its reported fixed-FF-solvent energy produces
a material 7CCN response, while barely changing the global result. The source gap
falls from 28.262193 to 24.120991 kcal/mol but remains large. This corrects the
objective used for the finite proposals; it does not demonstrate improved
biological accuracy, protein affinity or source-robust classification. Fixed
chemistry/occupancy, approximate electronics/solvent, consumed sources and
nonstationary constrained endpoints remain explicit limitations.

**This finite search-objective branch is closed. No further search or bound
extension follows from it.** The completed result records
`finite_objective_test_closed=true` and `same_search_automatic_extension=false`.

## Evidence and actual cost

Job 12769 completed with exit 0:0 in **3 seconds on 1 CPU: 3 allocated CPU-seconds**,
no GPU, no molecular calls and no new geometry. Slurm reports
`TotalCPU=00:01.216`, batch `MaxRSS=80304K`. This reporting cost is separate from
the [completed search/pool costs](../lbt_matched_response/FINAL_ACCOUNTING.json).
[ACCOUNTING.json](ACCOUNTING.json) pins the actual result, source, configuration,
plan, parents and receipts; [SACCT.psv](SACCT.psv) is the raw scheduler record.

| Artifact | SHA256 |
|---|---|
| [Result](../../../workspaces/lanm_coupled_response_20261016/monomer_matched_response_comparison_v1/execution/RESULT.json) | `71399f8fb20396a8aa6f52a6dc998033d6895fdd0dce82a6945c7c764adfe919` |
| [Executed source](../../../workspaces/lanm_coupled_response_20261016/monomer_matched_response_comparison_v1/MONOMER_MATCHED_RESPONSE_COMPARISON.py) | `eab60c068fd4418786ae8c93bc88925344c35d0bc729ed7b554e56bb6da6b778` |
| [Configuration](../../../workspaces/lanm_coupled_response_20261016/monomer_matched_response_comparison_v1/CONFIG.json) | `a6b7c450bf15590f7bc0a448ef5a33a6af28982825f3876e94717b6bde4af7dc` |

Values are copied from the actual saved result. No scientific computation,
submission or source change was performed for this report.
