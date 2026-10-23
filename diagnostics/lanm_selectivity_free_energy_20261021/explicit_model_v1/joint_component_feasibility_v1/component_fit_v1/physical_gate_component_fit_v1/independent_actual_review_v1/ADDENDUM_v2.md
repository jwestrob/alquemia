# V2: separate logical proposals from solver calls

Planck identified an internal accounting mismatch: `check_proposal` can recognize an explicit zero-call pre-solve failure, while the final v1 count equated all logical proposal records with real calls plus imports. V2 counts each reviewed proposal record independently, preserving separate counts for real solves, imports and explicit pre-solve failures. Their sum must match the logical count. The original eight logical, seven new-solve and one imported-proposal caps remain unchanged.

This is a defensive reviewer accounting correction. The frozen current DRIVER also rejects a later zero-call proposal before recording it, so no claim is made that a successfully completed current campaign can reach this edge case. No DRIVER, CORE, campaign policy, mathematical certificate, physical gate or native admission changed. Failed campaign results retain the separate non-admitting failure-audit path.

`REVIEW.py`, its READY, manifest and all prior controls remain immutable. The diff to `REVIEW_v2.py` changes only `run` counters/reporting and adds `check_proposal_counts`. Ten manufactured controls check first import, the full cap, explicit zero-call failure, mixed counts, corrupted totals, duplicate imports and excess solves/records; all other pre-existing function ASTs are identical. The earlier fifteen artificial controls are inherited and were not repeated. No campaign result/vector, real solver or native program was read or executed.

The terminal pass status and downstream result/configuration/DRIVER-implementation/selected/native-admission fields are unchanged. The previous CONTRACT and consumed-trimer binding guidance continue to apply.
