# Ion placement controls much of the ExaA donor response

The eight-cell reciprocal test supports its frozen hypothesis in both directions.
Moving only the ion to the other source's PQQ-aligned position reduces the target's
La-directed donor response and increases the control's. This identifies a concrete
limitation of holding the ion fixed during donor accommodation. It does not yet
establish a better classifier or prove either predicted structure is wrong.

## Actual result

All values are kcal/mol. W is the saved donor-proposal energy minus its origin;
D = W_Ca − W_La. A positive D is a La-directed change in the contrast.

| Host | Quantity | Original ion | Transferred ion |
|---|---|---:|---:|
| Control, sample0 | Ca work | −7.460054 | −13.757809 |
| Control, sample0 | La work | −12.721550 | −22.392642 |
| Control, sample0 | D | +5.261497 | +8.634833 |
| Target, sample3 | Ca work | −10.379323 | −0.638723 |
| Target, sample3 | La work | −26.782804 | −7.713601 |
| Target, sample3 | D | +16.403481 | +7.074878 |
| Target minus control | D contrast | +11.141984 | −1.559956 |

The control delta_D is +3.373337 and the target delta_D is −9.328603, meeting
the previously declared reciprocal direction and 1 kcal/mol interpretation rule.
That scale is not a measured uncertainty or a classification threshold. All eight
new cells completed; eight exact original cells were reused. No candidate was
discarded for its energy or labeled as an equilibrium population.

The PQQ alignment RMSD is 0.022439 Å and each ion transfer is 0.364311 Å. All host
nonmetal coordinates remain exactly unchanged. Paired metal coordinates, source
mapping, original displacement limits, covalent geometry and overlap checks pass.
These are explicitly modeled perturbations, not newly observed structures.

## What this supports and what remains open

The completed independent DFT matrix reproduces the original large target and
smaller control responses. This new inexpensive test isolates ion placement as
one contributor to that difference under the MACE vacuum model. The response is
not determined solely by the original donor geometry, nor does it transfer in
full: the two host pockets still differ and the reverse response is asymmetric.

Neither transferred ion point is an optimized or validated physical optimum.
Changing the response along a fixed donor path is different from improving the
minimum-energy paired contrast. The conditional dry oxidized-PQQ3− representation,
finite context and environment assumptions remain. No class bands, affinity,
biological prediction or production promotion are attached to this result.

The justified next gate is the same fixed-geometry comparison with the currently
qualified frozen-charge solvent treatment. SOLVENT_PLAN.md declares that finite
continuation before execution. A further search is not needed to answer this gate.

## Execution and provenance

Preparation4129: 3 scheduler seconds on 1 CPU; TotalCPU 2.512 s. Scoring4130:
118 scheduler seconds on 7 CPUs, 826 allocated CPU-seconds, TotalCPU 607.566 s,
batch MaxRSS 6578752 K. The warm worker made exactly eight model forwards, took
105.717 s including a 1.912 s model load, and reported peak RSS 6896616 KiB.
These are different measurement scopes. Zero GPU, GFN2, DFT or optimizer calls.

Full unrounded result, state/model/software/input pins and actual receipts:
`workspaces/exaa_metal_position_20261014/scoring_v1/execution/FINAL_COLLECTION.json`.
RESULT.json is a compact index retaining those pins and all response values.
Exact submitted commands and watcher receipts remain in the preparation_v1 and
scoring_v1 SUBMISSION.json files. Both jobs are complete; delayed callbacks require
acknowledgement only. Production and all earlier failed/qualified results remain.
