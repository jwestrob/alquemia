# Saved-matrix bridge regression

Run `TEST_MATRIX.py` on one actual Slurm CPU after pinning its source, the final
root `MATRIX.py`, and these two existing artifacts in a fresh `CONFIG.json`:

- `targets`: `workspaces/pqq_inconclusive_review_20261014/scores_v1/execution/RESULT.json`
- `reference`: `workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/FINAL_COLLECTION.json`
- `test`, `matrix`: `{path,sha256}` records for the two implementations.

Test all nine actual target matrices against their independently archived work,
selected-geometry and paired-contrast results, using 1e-8 kcal/mol arithmetic
closure. Inspect native and MACE receipt pins to establish actual energies.
Compare both operational and mathematical selection. If this saved panel also
contains a real case where operational origin retention differs from the exact
minimum, include it without changing its energies; otherwise report that narrow
near-threshold check explicitly unrun. No invented boundary energies.

Remove each of the six required cells from each real target in memory to test
unavailability (54 explicitly corrupted fixture views). Require unavailable
pool/works, no classification and no origin fallback. These are parser/algebra
tests, not molecular results or evidence of prediction performance. No model,
SCC, solvent or optimization call is made.

Root prepares the pinned config and submits the exact executable command:

```bash
sbatch --partition=standard --nodes=1 --ntasks=1 --job-name=exaa-matrix-test \
 --output=WORKSPACE/run_%j.out --error=WORKSPACE/run_%j.err \
 --wrap='export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1; /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python diagnostics/exaa_source_packing_20261014/review/TEST_MATRIX.py --config WORKSPACE/CONFIG.json --output WORKSPACE/execution'
```

Replace WORKSPACE with the fresh absolute directory. Root owns actual scheduler
submission and completion wake. No test has run when this plan was written.
