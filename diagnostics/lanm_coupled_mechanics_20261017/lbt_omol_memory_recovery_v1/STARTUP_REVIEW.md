# Independent startup observation:16680

Observed21October2026 at04:02:44UTC. **Actual allocation and runtime input match the resource-only recovery. No parsing error is visible.** This bounded review used saved metadata, text and hashes only and did not wait for initial-guess or SCF completion.

## Exact evidence

| Artifact | SHA256 |
| --- | --- |
| CONFIG_RUN.json | `908678ef7712be8280f35c57829d77556c3cbea6c8e7f5fc60e69b66e7a62c75` |
| ROOT_EXECUTION_DECISION.json | `c996c2e2573ff24bb0c57796d54a7f16896942e81087823a6488b3ee9ac952ed` |
| execution/START.json | `bc875ab7ab9ddd1ad788376948c8034d8d295c944819f1d3509c04207c417451` |
| execution/MANIFEST.json | `764da123fcc1dd3d7af4578fab0cedd68589078b305b8f10b8c15efa21df0654` |
| execution/SCONTROL_JOB.txt | `0c1be81a7382096c3f6eb3fde58a7c2751aea3daf5e6f2ed6c584c328b7d9d78` |
| central_La/endpoint.runtime.inp | `073e4a8d5c11d2de2119bc8e03a6879635030d24e9c9f5cf669d608214880d1a` |
| central_La/endpoint.runtime.json | `57f57556c2cfc3fbf7b560166dc71baa2008a03ea7177faa65c1f03fbe56fe43` |

Selected source/configuration/release/manifest/evidence/input/XYZ/basis pins match. The task signature exactly matches actual preparation16667. Source remains `5cad2623…`, evidence `67bc9bca…`, scientific input `4359537c…`, XYZ `55e07d81…` and external basis `7389c56f…`.

The actual runtime input equals16611's input after only `%maxcore3100→6200` and `nprocs24→12` substitutions. All other bytes, including PModel/NoAutostart, UKS, beta rotation `{514,515,20,1,1}`, EnGrad and scientific controls, remain identical.

## Allocation and observed native stages

Saved scheduler metadata reports RUNNING on one standard/standard_expand node, `node-64-768g-20`: **12 tasks,24 CPUs,two CPUs per task,96000MiB**. The runtime sidecar records one worker,12 MPI ranks/local slots,24 allocated CPUs, unchanged25% reserve and MaxCore6200MB per rank. The inspected launcher sets OMP/MKL/OpenBLAS/NumExpr/BLIS threads to1, and the passed startup allocation guard requires those exact environment values. This records12 single-thread numerical processes with24 charged CPUs, not24 active MPI ranks or measured utilization.

The live output was98784bytes/1838lines, modified04:02:19.461080UTC; no final immutable output hash is asserted. Native line54 identifies ORCA6.1.1;205/208 identify external orbital basis `def2-tzvpd.bas` and def2/J;224 again reports switching to SCNL. Lines1835–1837 identify external-file La ECP46,lmax3 at atom255. No explicit native parsing failure is present in this snapshot.

PModel and beta rotation are confirmed in the runtime input, but their executed stages are not yet printed. Molecular state/SCF residuals, actual peak memory, finite energy and gradient remain pending. RESULT.json is absent. Completion and health receipts at04:02:39UTC are RUNNING/watching for16680; health process91377 is armed. These receipts do not establish useful numerical progress or convergence.

No intervention is indicated by this startup observation. The prior OOM is not resolved merely by observing the intended memory policy: actual native completion, method/state/gradient admission and costs still require review. Preserve both earlier failures and the single-attempt scope; no automatic escalation or additional endpoint follows.
