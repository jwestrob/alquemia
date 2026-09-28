# Repaired full-region exchange: all six declared origins running

This comparison asks whether the larger electronic region plus actual protein field changes the Hans/Mex source reversal. It is a direct conditional model test, not merely a force-parser qualification. Both Hans crystal sources remain separate; no source averaging or calibrated classification is applied.

The matrix has six cells. Hans8DQ2 La/Dy run as1220332/1220323; Hans8FNR pair as1220334; Mex8FNS pair as1220336. Their automatic collectors are1220333/1220324/1220335/1220337. No new displaced geometries run. Each112rank endpoint uses existing native infrastructure with full available node memory (mem0), normal scheduling and terminal/live-health wakes. Exact commands/PIDs are in each workspace SUBMISSION.json.

The first table has0/6 completed cells; missing energies, contrasts and exchange differences are null. This is a preparation/collection result, not a molecular result. Three actual archived tests verify extraction, sign/units, missing values and duplicate rejection. No fabricated energies were used.

## Reproducible collection

Once the existing collectors finish, run from the repository root:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/metal_environment_response_20260926/repaired_exchange_20260928/COMPARE.py --config diagnostics/metal_environment_response_20260926/repaired_exchange_20260928/COMPARISON.json --output workspaces/metal_environment_response_20260926/repaired_exchange_final_v1.json
```

This refuses overwrite. It requires each accepted collection to pin its declared manifest/source/state and actual output/gradient/receipt artifacts. Both states of a source must use the same preparation; all sources must have matching electronic method, executable and orbital/ECP/aux hashes. Numerical guesses/SCF solver choices remain explicit. Raw C=E(Dy)−E(La) is reported separately from D=C(Hans)−C(Mex). Positive D means conditional greater relative La preference in that specified Hans state, not an absolute affinity.

Full preparation, scientific question and limits are frozen in ../electronic_diagnosis_20260928/REPAIRED_EMBEDDED_EXCHANGE_PLAN.md. Source-dependent water inventories, spectator occupancies, different QM regions, absent bulk-solvent terms and protein-level versus site-level evidence prevent a causal attribution or broad accuracy claim from signs alone. No classical component is silently included, no old compact energy fills a missing new cell, and no missing correction is set to zero. Preserve the completed compact reversal and all original-H results.

## Matching classical terms completed

All six source/metal classical origins now have actual source-mapped component energies and Cartesian gradients. Four Hans8FNR/Mex evaluations cost5allocatedCPU-s; two Hans8DQ2 evaluations cost3CPU-s. Hans8DQ2 archived repaired component energies reproduce within4.1e−12kcal/mol; no old-H arrays were reused.

Classical Dy-minus-La contributions are+0.0580692412 (Hans8DQ2),+0.0733472559 (Hans8FNR),+0.0667013621 (Mex), allkcal/mol. Thus classical exchange additions are−0.00863212085 and+0.00664589381 for the two Hans sources. Only target realQM–MM LJ differs; other classical terms cancel exactly. These small terms cannot on their own account for the earlier tens-of-kcal compact reversal. They do not predict the pending electronic contributions.

Use the separate COMPARISON_WITH_CLASSICAL_v2.json configuration with COMPARE.py and a new output path, for example workspaces/metal_environment_response_20260926/repaired_exchange_finite_final_v1.json. The script validates actual source/state/input pins, physical IDs/coordinates, component sums and electrostatic exclusions, then reports electronic-only and finite sums separately. No missing quantum or classical term is replaced with zero. Four real-artifact tests now include all six actual classical arrays with deliberately unavailable quantum collections; finite contrasts correctly remain null. Earlier partial artifacts retain their original hashes and are not regenerated.

The first classical configuration had an inherited electronic-only scope sentence; v2 corrects that explanatory text without changing cells or calculations. Original config/hash and first partial output remain preserved. Its exact intermediate analysis source is archived under workspaces/metal_environment_response_20260926/repaired_exchange_analysis_versions/COMPARE_classical_initial.py (SHA658d319b9cb6067dcda86f3de937909164dee82a962439c6aa70497bbbb7ce4f).
