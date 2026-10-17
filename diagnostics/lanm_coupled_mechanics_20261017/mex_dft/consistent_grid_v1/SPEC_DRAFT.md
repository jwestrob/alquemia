# Draft: Mex La O97 consistent-COSX-grid qualification

**Draft for root review; no execution or launch authorization.** This is a separate three-cell proposal. The [original eight-point protocol](../SPEC.md) and failed collector13271 remain immutable. Its remaining six points and eight response endpoints remain stopped regardless of this draft.

## Question and evidence

Does converging the density on the original final COSX grid restore the original La external-O97 directional derivative agreement at the original step? Saved diagnosis13277 localizes most of the −0.414201118901417 kcal/mol/Å mismatch to the post-SCF shift (−0.40307711778703115). Diagnostic subtraction still leaves −0.011124001114385607, outside the original gate. The printed COSX correction does not fully close that shift. This proposal tests a numerical consistency hypothesis; it does not establish a native implementation bug. See [result review](../gridx3_diagnosis_v1/RESULT_REVIEW.json).

The manual assigns grid2 to late SCF and grid3 to final energy/gradient evaluation; it describes modern COSX derivative terms but supplies no explicit guarantee about external `.pcgrad` response to the post-SCF correction. [Official RI/COSX documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/RI.html#cosx-numerical-grids). Density nonstationarity for the final grid is an inference; SCF/root and finite-step uncertainty remain.

## Exact source lineage

Paths below are relative to the repository root. Bind all dependencies and exact new decks before any execution.

| Source | Path or reference | SHA256 |
|---|---|---|
| Original specification | `diagnostics/lanm_coupled_mechanics_20261017/mex_dft/SPEC.md` | `6ca3f276638365a58b026e60d82365941e4e6be8ef4de1c329ae87bb5dd1524c` |
| Saved diagnosis | `workspaces/lanm_coupled_mechanics_20261017/mex_gridx3_saved_diagnosis_v1/execution/RESULT.json` | `08e4d217dbe5a6391afadc09ddc8aa6a84c51e2eb7d316b17e6b9362b2a508a4` |
| Central native template | `workspaces/metal_environment_response_20260926/lady_repaired_origins_Mex8FNS_v1/La_A/endpoint.inp` | `c8d5db2987858632c74baa1cadf48427c5cb37642bcfc18d594cbdfaa2ca8d9f` |
| Original plus/minus SP template | each original scout's `execution/endpoint.inp` | `e8e542380900fd0e7e6a694fc85504a5b679821649feb8e6697da54db64059d8` |
| Original preparation | `workspaces/lanm_coupled_mechanics_20261017/mex_dft_preparation_v1/execution/RESULT.json` | `fc97c1feb145855df7c9ae0a3e93a179389ae07df445dfab8571cc4cb70d8d5d` |
| Frozen directions | same preparation directory, `DIRECTIONS.json` | `166b8ea4e19ec29c6d64b5bbd285178b52685c598a3367457e95f2bda387f4a6` |
| ORCA6.1.1 binary | original endpoint CONFIG `orca` entry | `38b5f057452fef275c0a1b98d270ad03d0411dab70453820d1c678bd732f6c83` |

The diagnosis pins all three original outputs, fields, mappings and runtime receipts. Reuse exact central and prepared ± fields/core/physical geometry; replay their joins and byte identities. Core SHA256 is `7d3790074fb01647c0a892bfe9ae0eca0376575799edc456cf75bb7aac66d9b6` in all three. Only external row832, A/97//O, physical1019, is displaced. Its frozen direction is `[-0.610565917407258,-0.40030480677694336,0.6833486095485287]`; h=0.001 Å. Preserve all charge values and original geometry/cap/charge-shift rules.

## Sole numerical candidate: grid2 equals original grid3

Actual original central and ± outputs report:

| Integration | Original angular / maximum Lebedev | Original printed IntAcc | Original points | Proposed effective setting |
|---|---:|---:|---:|---|
| XC | 6 / 590 | 4.959 | 2341093 | unchanged |
| COSX grid1 | 2 / 110 | 4.020 | 256114 | unchanged |
| COSX grid2 | 3 / 194 | 4.338 | 549243 | original grid3, angular4 and its default radial construction |
| COSX grid3 | 4 / 302 | 4.871 | 1026059 | unchanged |

Preserve OptM3/GC(2021), adaptive pruning4, mBecke2022 weights, BFCut1e-12, WCut1e-14, H/He angular reduction, rotationally invariant construction off, COSX partial contraction and UseSFitting on. These are observed native settings. Keep XC and final integration enabled.

Proposed replacement for the existing `%method` block only:

```text
%method
 DoEQ false
 AngularGridX 2,4,4
 IntAccX -1,-1,-1
end
```

`AngularGridX` controls the SCF COSX angular grids; `IntAccX=-1` requests the angular grid's default radial accuracy. The official table associates angular4 with COSX IntAcc4.871. Using its default preserves the original radial construction without assuming that a three-decimal printed IntAcc specifies every internal parameter. [Official grid controls and radial defaults](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/numericalintegration.html#scf-grid-keyword-list).

**Installed parser and default-reset behavior have not been tested.** The native output must demonstrate grid2 and grid3 angular, radial, pruning, cutoff, weighting, fitting and point-count equality, with grid1, grid3 and XC matching the original fingerprints above. A matching IntAcc label alone is insufficient. Reject an incompatible execution; do not silently substitute syntax or settings. This draft contains no prepared input or parser run.

## Three cells and preserved admissions

1. Fresh central La `EnGrad`, including native QM gradient and external `.pcgrad`.
2. Fresh plus La singlepoint at the exact original +h field/geometry.
3. Fresh minus La singlepoint at the exact original −h field/geometry.

Copy the native templates, introducing only the grid block above; remove only `EnGrad` for SP. Preserve PBE0-D4, def2-TZVP/def2/J, all custom La basis/auxiliary/ECP46 bytes, RIJCOSX, NoAutostart, PModel, DefGrid3, VeryTightSCF and DoEQ=false. Preserve native SCF tolerances, maximum125 iterations, Total+1el-Energy/ConvCheckMode2, ConvForced0 and solver defaults. Do not enable restart transfer, change step, tighten SCF, replace final energies, subtract a correction, disable final grid or select a new state by outcome.

Retain Q209 with1840 external charges, electronic charge0, multiplicity1,832 explicit electrons/416α416β, physical La f0, spin-free frozen-f representation and SOC=false. Preserve full assembly2049/native2045, three Nd spectators, original water/state assumptions and fixed-FF solvent bookkeeping. `physical_state_qualified=false` remains explicit.

All three use ORCA6.1.1 and the existing pinned runner/renderer/parser on one node with12 actual MPI ranks,1 CPU/rank and the same allocation policy. MaxCore derives from the actual allocation using the existing reserve; require matching CPU/memory-derived settings across cells. Prior12-rank scouts used48000 MiB and MaxCore3100 MB. No scheduler memory/time override is proposed. Each cell uses a fresh directory. The historical central used112MPI/MaxCore10800; the new matching central is essential, and a pass would not isolate grid effects from every historical runtime/root difference.

Admission retains native normal termination, native SCF convergence and exact source/state/runtime checks. Record every printed residual, tolerance, stopping route, gap/warning and all-residuals flag; false flags remain observations under the original admission policy. Native convergence is not electronic stability. [Official SCF convergence semantics](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html#convergence-tolerances).

## Frozen comparison, finite stop and cost review

Use the **new central gradient and both new final native energies**, with exact printed-string subtraction in native units before conversion. Retain the original gate for PBE0, MM, fixed-FF ALPB and total:

`abs(FD−analytic) <= 0.005 + 0.0001*max(abs(FD),abs(analytic))` kcal/mol/Å.

Reuse compatible archived MM/ALPB energies and gradients after exact geometry/field/state joins; no new MM/ALPB calls. Report components separately and retain printed quantization, binary64 scales, labeled SCF sensitivity proxies and final-grid closure diagnostics. SCF proxies are not error bounds. Step-halving/truncation remains unqualified.

Proposed staging is central first, then the two SPs only after native/grid admission and root's actual-cost review. Maximum is three molecular calls, no retries or added cells. Any native/source/grid/state failure or scheduler termination stops dependent cells; final derivative failure closes this candidate. Output silence triggers review only. There is no automatic SCF fallback, grid sweep, optimizer, step change or remaining-six/response continuation on either outcome.

Measured original costs were26666s×12CPU for plus,22132s×12CPU for minus, and6s×1CPU for collector:585582 allocatedCPU-s. Saved diagnosis added2 allocatedCPU-s. These are sunk costs, not predictions for the denser SCF grid or new central EnGrad. Before finalizing an executable specification, root must review this cost context, finite staging/stops and actual scheduler limits; record the approved resource ceiling/policy and bind the final manifest. Keep actual allocation cost distinct from measured process CPU and count each allocation once.

A pass would establish only this new recipe's one La external-charge directional check at h=0.001 Å. It would not rehabilitate the original failed protocol, qualify Dy/boundary/global forces or stable electronic states, establish physical accuracy/affinity, or promote a production default.
