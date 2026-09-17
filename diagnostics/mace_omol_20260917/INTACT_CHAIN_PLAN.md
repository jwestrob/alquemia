# Next experiment: intact-chain OMOL coordination response

Declared before implementation/inference on these complete chains,2026-09-17.
Jacob's renewed discretionary analysis/resource authorization applies. The
failed core coordination candidate remains frozen. Do not replace its inputs
or reinterpret its results as a prospective test of this new representation.

## Scientific question and model

Does supplying intact protein connectivity and local structural context,
without synthetic peptide caps, improve the same matched coordination
descriptor? This is a new representation of the same five previously consumed
structures, not a new functional or a fit to the failed core scores.

Use native MACE-OMOL-0 100M, the existing pinned checkpoint and float64 software,
with the actual full-system charge/multiplicity. Retain the completed global
panel's physical preparation exactly:

`workspaces/mace_global_benchmark_20260916/prepared_v1/preparation_manifest.json`
SHA256 `63d6a44fbef3a1dc888696bba5d122ac3531fe12ceee4ee03a85194f768b0552`.

| Case | Atoms | La charge | Ca charge | Evidence stratum |
|---|---:|---:|---:|---|
| ALPHA_1F6S | 1932 | -4 | -5 | qualified qualitative affinity |
| ALPHA_6IP9 | 1898 | -4 | -5 | same biological group, alternate source geometry |
| GGR_1GLG | 4698 | -3 | -4 | direct same-assay affinity |
| PQQ_1H4I | 9088 | -9 | -10 | consumed PQQ functional class |
| PQQ_4MAE | 8854 | -3 | -4 | consumed PQQ functional class |

All multiplicities1. Use the archived physical preparation verifier and exact
source-to-atom inventory. These are **whole prepared chain A** models, not a
claim to include every biological-assembly chain. Earlier hydrogen-bond-length
repairs, recorded terminal OXT completion, PQQ microstate, and explicit-water
inventory stay fixed. Do not reprepare them or retain synthetic carve caps.

For each metal, evaluate the original bound geometry and a derived geometry
moving only the selected metal to x=max(other atom x)+30Å at unchanged y/z.
Use R_coord=(E_bound,Ca-E_detached,Ca)-(E_bound,La-E_detached,La), converting eV
once. Same atoms, charge and spin in each subtraction. No solvent, external
charges, new polarization model, aquo offset or local relaxation is added.
The finite-range model does not describe full long-range protein electrostatics.
Its separated fragments have no certified individual ionic states. Whole-chain
sizes exceed OMol25's reported training-system sizes even though these total
charges are within the reported dataset range. Treat this as a descriptor test.

## Engineering qualification before comparisons

Use the native forward with `compute_force=False` under `torch.no_grad` to
avoid retaining the gradient computation graph. This must not modify weights,
precision, graph edges, energy terms, atom order or charge/spin batch. Retain
actual native node/embedding terms and their sum checks. Missing gradients are
explicitly `not_requested`, never invented zero forces. A success from the old
force-producing cache cannot satisfy a new energy-only task without a declared
validated reuse mechanism.

Stage A: four exact original canonical cores,1H4I/4MAE×La/Ca. Require native
energy-only output to reproduce the actual force-producing endpoint energies
within0.01kcal/mol, finite components, exact input batch and unchanged weights.
No full-chain scoring if this engineering equivalence fails.

Stage B: ALPHA_1F6S, bound and detached×La/Ca: four primary calls, four repeats,
four existing deterministic rigid rotations about the selected metal, and two
detached references with a further+10Å metal-x displacement:14calls. Require
repeat/rotation/distance endpoint and paired-contrast errors<=0.01kcal/mol;
verify actual metal graph-edge count0 for detached inputs, and exact state and
component accounting. These are energy checks; no force/response validation is
claimed. Any memory recovery must preserve this same scientific model/input.

Stage C, if A/B pass numerically, regardless of preliminary direction: remaining
four cases×bound/detached×La/Ca=16calls. Reuse the four ALPHA_1F6S primary
endpoints. Total34 new native forward evaluations,zeroDFT,solver,optimization,
training or gradient/Hessian calls. Each failed attempt remains recorded and
any technical retry is explicit; a failed inference is not a numerical value.

## Frozen predictive assessment

No absolute calibration is available for this new whole-chain representation.
Do not inherit core or baseline bands. Test exactly three relative differences,
each required>0.02kcal/mol:

1. PQQ_4MAE minus PQQ_1H4I, functional-class ordering.
2. ALPHA_1F6S minus GGR_1GLG, qualified affinity ordering.
3. ALPHA_6IP9 minus GGR_1GLG, same biological pair on the other alpha geometry.

Require all three for this small development gate, retain each separately, and
report all missing inputs/endpoints. Compare against the archived core OMOL
and core matched-contrast results without selecting a favorable GGR core.
This consumed panel cannot establish broad validation or a calibrated absolute
decision even if it passes. Do not change labels, water states, charges, source
geometries, cutoffs or thresholds to rescue a result.

## Resources and implementation boundary

First try native energy-only inference on the existing A5000 allocation policy:
one GPU,16CPUs,64,474MiB host RAM. Full-chain native peak memory/runtime is
unknown; do not extrapolate the small-core3GB measurement as a guarantee.
If it does not fit, use the already authorized H200 allocation with its actual
scheduler-enforced1/8host-memory entitlement; multiple GPUs may reserve more
host memory if required, with all allocation costs reported. Do not claim extra
GPU reservations increase a single GPU's VRAM. Preserve others' jobs and queues.

No project compute budget or time cap; retain scheduler QOS requirements and
finite task manifests. Aim for practical minute-scale inference, measuring
startup, validation, peak memory, failed attempts and total allocations. Do not
claim affordability until measured. Exact graph tiling is a possible subsequent
technical route, but is not silently enabled here; it needs its own equivalence
checks before replacing the native implementation.

Use existing runner/snapshot/cache/receipt conventions and a distinct protocol
`mace_omol_intact_chain_matched_coordination_v1`. Production and old experiments
remain unchanged. End-to-end feasibility is an intermediate milestone, not
completion of the broader discriminator goal.
