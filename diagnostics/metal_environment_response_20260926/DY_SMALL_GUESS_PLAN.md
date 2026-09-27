# Two-start isolated native Dy capability diagnostic

Declared27September before new energies. Reuse exact consumed50atomHansEF3Dy
workspaces/lanm_series_followup_20260923/prepared_v1/endpoints/Hans_EF3__Dy/core.xyz.
No geometry/proton/water/charge changes; this older carve has its existing H
preparation and is a local engine diagnostic, not the new195atom model.
Native r2SCAN-3c/DefGrid3/TightSCF/EnGrad, no solvation or external pointcharges,
charge−1, physicalsextet, nativeDyECP28. Two explicit initialguesses PModel/HCore.
No basis/ECP/spin workaround, PAtom (known unavailable), iterations increase,
wholeprotein vacuum GFN2 retry, relaxation, calibration or biologicalscore.

Question: can this native electronic state converge on an actual much smaller
local core, and does initialization change the outcome? Both195atomDy attempts
had catastrophic density excursions beforeTRAH. AtomicDy PModel construction
succeeded; PAtom is unavailable in actuallegacyORCA. HCore is a supported but
crude alternative, not a promised remedy. Success here would still confound
region-size and embedding changes; it cannot uniquely identify their effects.
Failedattempts and both declaredstarts remain visible; no selecting by desiredsign.

Two endpointcalls only, existingnativeexecutor and strictstate/gradient parser.
Record actualRHF/UHF, ECP, electrons, S²/localspin evidence, perstartenergies/forces,
convergence/cost. Agreement does not prove physicalstate correctness; disagreement
requires diagnosis, not choosing the lower favorableclassification.
No dependent expansion until inspection. No arbitrary walltime/projectbudget;
live numerical-health monitor must wake root during pathology, in addition to
terminalcollector/wake. Use available shared24CPU slots,2x12MPI, --mem=0,normal
priority; all currently free slots on the mixed48CPU node were24 at preparation.
Never claim fullnode utilization. No touchingothers' allocations orPQQ priority.
Stop starts after Monday10America/Los_Angeles cutoff. Preserve allartifacts.
