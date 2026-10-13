# Two-Ca sources: unchanged local scoring

Research-only source-state extension. Implementation may proceed against the
agreed preparation schema; executable source manifests require actual admitted
preparations. No calculations or new calibration are implied by this plan.

The implementation is `NATIVE.py` (metadata preparation, seven-CPU lane execution,
and collection) plus `SOLVENT.py` (metadata configuration and unchanged corrected
solvent stage). The three Slurm wrappers take explicit manifest/config paths.
Root submits each native lane, the one-CPU collector with `afterany` on every lane,
then solvent with `afterok` on that collector. Collection preserves terminal
source records and explicitly lists missing source records after failed lanes;
its terminal flag denotes a completed collection snapshot, never scientific
success. No native worker or optimizer is automatically restarted.

## Inputs and identity

Curie's complete `FINAL_PREPARATION.json` retains every declared sample, including
failures. A prepared row supplies `case_id` (new `conditioned_2Ca_v1` source),
`original_source_id`, `selection_id`, source class/group, context preparation and
paired endpoints, maps, state signature/key, and spectator-admission receipt.
Require the original canonical fragment union and frozen state anchor, admitted
physical maps, and the explicit spectator boundary check. The full source keeps
both actual Ca ions; the local model represents only its admitted atoms. It is
not an electronic or mechanical model of the omitted spectator.

Q9Z4J7/ExaA and Q88JH5/PedE retain the original Ca reference label; C5AXV8/ExaF
retains La. Recover the exact group and selection from the original125-source
scope through `original_source_id`; do not infer them from a filename. A new fold
is an alternative source hypothesis, not a repaired copy of the old coordinates
or an independent protein. Retain every declared sample; do not select a winner
or introduce a new triple/ensemble policy.

## Native interface and execution

Reuse the qualified `CPU_BATCH.WarmCPU`, `run_source`, `proposal`, and `native_cell`
and pinned archived `motion_envelope_run.build_tasks` / SLSQP implementation.
Generalize only the disulfide metadata and collection adapter. Two or three
independent seven-CPU Slurm jobs each hold one calculator warm across their
assigned sources. Assignment follows manifest order; each lane owns its own
worker manifest, source directories, locks, attempts, and receipts. Preserve the
two real force/energy qualification fixtures per lane and their existing gates.

Each source obtains two paired origins, the same common four-mode force-selected
subspace, one bounded search for each metal, and both cross-scores. Retain the
origin/Ca-proposal/La-proposal pool, exact coordinate aliasing, paired charge and
physical-map checks, optimizer/head/checkpoint/float64 settings, and missing-cell
behavior. Existing worker authorization covers root/proposal/cross manifests;
the historical GPU worker does not accept that interface unchanged. No new GPU
transport or backend qualification campaign is part of this recovery.

## Corrected solvent and frozen decisions

Reuse the pinned12-worker `GFN2_STAGE.run`, `GFN2.py` component verification and
collection, and corrected native xTB / independent tblite component-v2 settings.
Up to12 distinct inputs per backend per source are needed (two metals, three
geometries, vacuum/water); deduplicate only exact compatible input identities.
Do not claim molecular reuse from another fold or from matching protein IDs.
Preserve any interrupted/failed point directory, and never restart it silently.

Use the exact canonical record with SHA256
`929dedab7ba46c7088de7a4a695a8f19fc5af60066ecd3599645f687bd4e6896`.
Operational/mathematical bands are Ca<=-405457.564455343 and
La>=-405454.092284296 kcal/mol; origin bands are Ca<=-405465.7535534819 and
La>=-405455.45587457. Values between a policy's bands are inconclusive.
The original absolute-v1 backend parity failure remains separate and failed.
No threshold fitting, biological relabeling, or unavailable-to-success fallback.

## Output contract

Keep original125 outcomes unchanged. Link each new result by `original_source_id`
and retain its new source ID, preparation/state pins, selection, all six cells,
individual accommodation works, origin/operational/mathematical decisions,
failures, and actual compute receipts. Preserve the existing three-source MxaF /
MMOL recovery overlay separately. Root integrates actual results into a clearly
named alternative-source overlay; no new fold/preparation status implies a score.
Report runtime from actual lane/stage receipts, not GPU extrapolation or a single
source timing multiplied into a benchmark. Production remains unchanged.
