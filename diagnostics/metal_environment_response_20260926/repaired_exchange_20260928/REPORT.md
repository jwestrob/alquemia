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
