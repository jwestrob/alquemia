# Technical reader recovery v2

Job3353 failed before producing a result because the archived, unsupported La sample1 stores unavailable metadata in its `local` field. The v1 reader treated those strings as local geometry records. No chemical failure or scored result follows from that exception.

`ANALYZE_v2.py` changes only this reader branch: a source without admitted endpoints retains its raw PQQ metrics, an empty local-metric collection, and explicit `unavailable_original_preparation` status. The same ten raw sources, nine admitted origins, three crystal references, chemical graph and numerical tolerances remain fixed. Original driver, manifest and failed receipt are unchanged.

The fresh manifest and executable command are pinned in `READY_v2.json`; output is under `workspaces/pede_solvent_response_20261013/cofactor_review_v2/execution/`. Input hash and Python AST checks ran on login; scientific analysis/tests did not. Root owns the one-CPU Slurm execution and wake. This is a parser recovery, not a new chemical model or a successful scientific result.
