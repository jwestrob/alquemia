# Independent terminal review of job 15525

**Admitted model convergence; no blocker found.** Actual RESULT and SEARCH report `full_force_converged` at accepted step 725 / distinct evaluation 739. The raw full-gradient maximum atom force is **0.0750547504** and atom-vector RMS is **0.0296109468 kcal/mol/Å**, satisfying the unchanged maximum 0.10 / RMS 0.03 criteria. Source inspection confirms that convergence uses the full gradient, before the rigid-projected gradient is used for direction construction.

## Scope, lineage and guards

- The run made 567 complete MACE and 567 complete fixed-FF ALPB calls, yielding 553 accepted steps and 14 rejected evaluated trials. No DFT, GFN, MM, opposite-metal or other-search calls are reported. All 567 native ALPB execution receipts return zero; all MACE input-state records retain 256 atoms, charge −2, multiplicity 1 and target index 255.
- The exact bound source/config and saved audit remain intact. The original initial point, source and recipe are preserved; the new final checkpoint retains the exact prior request/evaluated-point prefixes plus the new records. It contains 739 requests, eight history pairs and trust 0.1 Å. The 193 pre-existing point records are retained within the 760-record final inventory.
- All accepted trace states retain the original physical/state guards and all 256 movable atoms. The final state has no recorded bond, chirality, peptide-basin or severe-overlap violations and no state-departure labels. Maximum displacement from the original source is 4.567726 Å; the unchanged protocol explicitly records source displacement as diagnostic, not a gate.
- Terminal positions and full composite gradient are both float64 arrays with metadata shape `(256, 3)`. All request recipes equal the saved admitted recipe, with no extra field or MM term. The terminal trace has no further trial.

Inspection verified 9,280 selected file pins covering the bound source/config, point/request lineage, native component receipts/artifacts and terminal arrays. Only source, JSON, hashes and array headers were inspected; no gradient, energy, optimizer or molecular calculation was rerun.

## Interpretation

The additional step172-to725 work is MACE **+11.053362**, ALPB **−24.547593**, composite **−13.494231 kcal/mol**, as stored by the producer. The run reaches the declared model force tolerance despite this opposing component behavior. It does not prove a physical minimum, positive Hessian, electronic-state stability or force accuracy; those qualification flags remain false. All prior SCF/reference limitations remain.

The application reports 403.841720 seconds, including 69.992559 seconds of MACE calls, 28.555505 seconds of ALPB calls and 6.410947 seconds of model loading. Runtime metadata retains 12 CPUs and one visible GPU on the qualified host. These application timings are distinct from root's scheduler accounting.

## Proposed terminal reference

**One exact whole256 La PBE0 energy/gradient at the terminal POINT is a useful next mechanics test.** Pair it with the admitted step92 reference 15321 and reuse both saved fixed-FF ALPB endpoints. This tests the later model lowering and whether the reference-plus-ALPB residual force is also small at the model-converged geometry.

Bind the exact two endpoint identities before execution; preserve the original method, state, ECP, grid/SCF settings, 12-MPI contract, native admission and 1 kcal/mol reporting resolution. Admit the cached/native chain before the single new call; preserve failures without retries or settings changes. No further optimization, other metal or other geometry is implied.

Report vacuum and composite work separately, along with the full terminal reference-plus-ALPB force and descriptive model/reference force comparison. Shared ALPB means that the work discrepancy compares the electronic models; it cannot independently validate solvent physics. The step92 reference retains its EnergyCheck termination and failed printed orbital-gradient/density flags. Native admission must not be relabeled as all-criteria SCF convergence or a qualified force oracle. Even a small terminal reference gradient would not by itself establish a minimum.

## Exact result pins

- RESULT: `95d870ebbfdedede5166c7499e4b8181dd8c24289503bdb2e35eaba02c4a5068`.
- SEARCH: `a9495fd5f42c6fc115a26e00cd33dd9d2369ce4b7ffcdcb4fe5aa6ee0115c25d`.
- Final checkpoint: `0461e30d7b285cee7b659de4865d27688bfb2f936edd55f3b57a4628b8df88c1`.
- Terminal POINT: `6605f43cc3289509cee65c9be97a195abaad43d2d247ab488f5de889aa6f5419`.
- Terminal MACE: `bc57d82692d923aa2321e1780961e4e355fbb608ff8d528ffc4f3286d79b9a2a`.
- Terminal ALPB: `9a3b7c6e4ae290ad00c4f81ca19ceeef3306df20b21b420d8e9ced85fabb98cc`.
- Terminal full gradient: `cee38bd2763c48b9b9023b1f330a3601432b7bc38f28c650d533d71888a27cff`.

Raw inputs are under `workspaces/lanm_coupled_mechanics_20261017/lbt_7ccn_generous_completion_v1/run/execution/`. The frozen implementation/config are `1fc6b168` / `1bfaebc0`, reviewed in `STATIC_REVIEW.md`.
