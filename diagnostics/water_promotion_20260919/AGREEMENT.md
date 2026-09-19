# Promoted contextual water preparation

Jacob explicitly authorized promotion on2026-09-19: “go ahead and promote it.
i authorize you to go ahead and undergo the next round of experiments at your
discretion ... make sure you have an eye towards improving the classifier.”

Promote the demonstrated contextual water-H preparation into the versioned
baseline workflow for supported explicit prepared sites. Its default is
`contextual_if_supported`; `original` explicitly retains the historical path.
Both original and prepared quantum scores remain separately recorded. The
electronic method stays native r2SCAN-3c/CPCM(Water)/DefGrid3 single points.

Eligibility requires source-backed peptide-amide-v3 cores, an explicit validated
protein assembly and complete existing site waters. The frozen preparation uses
source-defined3.5Åpolar neighbors and the demonstrated native MACE optimizer.
Only inner-water H coordinates transfer to the original DFT core. No new waters,
oxygen motion, occupancy, proton transfer, entropy or changed cofactor chemistry.
Dry inputs are exact identities. Unsupported wet cofactor/legacy representations
are explicit unsupported states, not successful preparation or silent fallback.

New workflow ID: `baseline_contextual_water_v1`. Prepared DFT protocol remains
`amide_v3_native_r2scan3c_context_prepared_water_H_v1`; original protocol IDs,
inputs, outputs and published references stay accessible. Changed wet geometries
do not inherit old absolute decision bands. Unchanged canonical dry PQQ retains
its own released decision policy.

The supported prepared-site workflow is exposed through the existing
`scripts/affordable_workflow.py baseline` entry point and existing GPU/ORCA
manifest executors. It does not restart watchers, mutate running jobs or rewrite
historical inbox workspaces. Raw legacy inbox cores need compatible source-backed
preparation before this component can operate; promotion does not imply that
unsupported inputs acquired a new score.

Verify real archived alpha inputs/proposals/DFT outputs end to end, retain
original/prepared fields, check dry PQQ identities and actual unsupported/failure
behavior. Reuse successful chemistry; runnable new-case execution must not require
editing source code. Document practical scope and scientifically untested transfer.
