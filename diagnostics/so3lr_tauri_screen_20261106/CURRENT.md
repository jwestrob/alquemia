# Bounded Tauri screen — closed 6 November 2026 UTC

Decision: unresolved reference; no MD pilot. Read REPORT.md, frozen PLAN.md and
SCIENCE_STOP.json. Preserve M1 and end automatic expansion. No numerical job,
reference retry, alternate solver, peptide or whole-protein continuation is
released by this checkpoint.

All24 exact archived configurations were evaluated with L/M-2606:48 primary
predictions,12 identical-coordinate La/Dy pairs per model. Both models pass
native/ASE, physical-state propagation, fragment metadata and energy/force checks.
Raw float64 arrays, control/sensitivity data, inputs and hashes are retained.
L CW/CAc−1.113/+3.963; M+2.560/+10.769 kcal/mol, electronic contrasts only.
Distant spectator and capacity sensitivity remain material concerns.

Exactly six reference invocations49309–49314. La49309/49310/49313 finished with
finite gradients, correct ECP/counts/S2 and positive stability eigenvalues, but
failed frozen printed Tight density/DIIS criteria. Root declared SCIENCE_STOP
before cancelling only unfinished Dy49311/49312/49314 after2686s each. Continuing
Dy could not qualify either contrast because both require the unqualified LaW8.
Their preserved partial outputs contain no admitted final energy/gradient. This
is an unresolved reference screen, not evidence of Tauri accuracy or failure.
No restart or changed criteria. All23 owned Slurm jobs are terminal.

Native audit verifies599/600 at1M using trans's alternate row and fallback576/597
at.1M, but neither common experimental medium is admitted: missing exact anion/
source provenance, with fallback adjustment mismatch. No third pair search or
reserved outcomes. Full native rows stay in the ignored workspace.

Cost150814 allocatedCPU-s (41.89h),146636.412 measuredCPU-s (40.73h), zeroGPU,
including two failed writers, stopped references and the successful saved replay.
Warm29-atom CPU evaluationL.1262s/M.0644s; no H200 or bulk-solvent throughput claim.
GPU49306 was cancelled pending QOSMaxGRESPerUser. Other-session jobs are untouched.

The saved-data reproduction49465/49466/49467 reproduces all five numerical tables,
reference rows and decisions, preserving original artifacts. See
REPRODUCTION_CHECKS.json and ARTIFACTS.json. Delayed completion events for GPU
scout49306 or six-reference watcher are ACK only: WAKE_ACK.json already records
root consumption and closure. Do not start another allocation from those events.
