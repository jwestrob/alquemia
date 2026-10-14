# Completed ExaA source/packing investigation

**Both tests are complete; no owned jobs remain.** Read REPORT.md. Crystal1FLG A
classifies Ca under the unchanged accommodated bands; B remains unavailable.
The crystal's accommodation contribution is+1.929997kcal/mol versus+16.691741 in
the ambiguous two-Ca fold. Expanded packing leaves the latter at+16.316953.
This establishes model source sensitivity and weakens the finite missing-packing
hypothesis; it does not establish that AF3 is wrong or that MACE response is right.
An independent same-coordinate electronic-response check is recommended next,
but no further experiment/DFT campaign is submitted by this completed round.

Authorized14October2026 after the inconclusive-case discussion. Read PLAN.md.
No owned Slurm jobs were active at startup. Previous cancellation notifications
do not request reruns.

- Crystal agent: audit both 1FLG chains, compatible state and runnable preparation.
- Packing agent: prepare exact saved two-Ca ExaA0/3 pools with complete added
  all-heavy contact fragments; no optimization.
- Review agent: audit reuse of current MACE/frozen-charge solvent machinery.
- Root: finite shared scoring, submissions, collection and interpretation.

Packing preparation4089 passed both sources and all12 endpoints:192→249atoms,
same Ca−3/La−2 singlets, exact source/motion preservation, common expanded state.
See packing/REPORT.md. New fragments include two complete peptide groups and
Phe442/Leu485/Trp591 sidechains. The experiment tests total context response;
it cannot uniquely attribute a change to Leu or steric resistance.

Saved-artifact regression4090 passed all nine real target matrices and54
explicit missing-cell fixture views. The optional near-threshold origin-retention
example was unrun because no such case exists in the saved panel. No synthetic
energies or new scientific evaluations were used in these implementation tests.

**Molecular job4091 completed**,681s/7CPUs on the standard partition. Its finite
12-cell common pool uses two real warm-MACE qualification controls, then all
12 new MACE endpoints,24 native xTB inputs and24 independent tblite inputs,
followed by the unchanged native frozen-charge solvent replay. No optimization.
The existing component/parity policy remains exact; absolute xTB/tblite offsets
are diagnostic, not silently adopted as an acceptance criterion. Both pools and
all12cells are available. WakePID235236; exact command/config/test pins
are in `workspaces/exaa_source_packing_20261014/packing_score_v1/SUBMISSION.json`.
Comparison4092 completed on oneCPU, retaining all candidate works,
component differences and unavailable results. Both terminal events were inspected;
see `packing_comparison_v1/SUBMISSION.json`. It cannot classify the enlarged
context using the current reference bands. In control0 the selected Ca-minus-La
accommodation work changes5.310100→5.063378kcal/mol; target3 changes
16.691741→16.316953. Added context restrains both metals similarly and preserves
both selected geometries. This finite test does not support omitted packing as
the main cause of the target's large differential response. It does not establish
an expanded-context classification or validate the MACE response independently.

Crystal preparation4088 completed but both chains remained unsupported by its
strict no-missing-atom gate: each crystal lacks one terminal heavy atom and has
zero missing internal heavy atoms. Crystal-only recovery4093 reuses those exact
outputs and admits chainA: the sole addition is terminalTRP616/OXT,41.494Å from
the context, with every deposited heavy coordinate unchanged. No protonation was
repeated. ChainB remains unavailable because its mapped Asp-partner contact is
3.511Å versus the existing3.5Å preparation cutoff; that is a workflow gate, not
proof of a defective crystal. Both rows and the original4088 failures remain.
The actual v3 preparation is in crystal_v3/execution/FINAL_PREPARATION.json.
**Crystal native4095 and solvent4096 completed on7CPUs.**
All6 A matrix cells are complete; both unchanged searches converge and their
proposals pass physical geometry checks without active-bound flags. Runtime,
source mapping and both real model qualification gates pass. The native job
made33 actual forwards in427.602s (scheduler428s), including the two controls.
Solvent4096 takes128 scheduler seconds, reuses all MACE receipts and completes
all6 solvent cells with the same published-parameter GFN2/native frozen ALPB
method. Both metals select the La-generated geometry; their works are−5.169804
(Ca) and−7.099801(La)kcal/mol. The raw change is La-directed, not Ca-directed.
Static inconclusive→accommodated Ca uses distinct previously frozen canonical
bands and is not evidence of a Ca-directed physical improvement. B remains in
the denominator without molecular calls. All component gates pass; see
crystal_scoring/RESULT_REVIEW.md. No calibrations change.
Oxidized-PQQ redox remains an explicit conditional model assumption, not an
experimentally established matched state. Source disagreement cannot uniquely
diagnose a bad fold or correct model.

Exact native config/submission/wakePID254783 are under `crystal_native_v1/`;
the solvent submission/wakePID254785 is `CRYSTAL_SOLVENT_SUBMISSION_v1.json`.
Both use the actual current root thread and separate durable completion events.
All terminal receipts/results acknowledged; do not rerun from delayed wakes.
Read-only collector4106 completed on oneCPU, with wakePID267070;
exact command is `crystal_saved_comparison_v1/SUBMISSION.json`. It retains all10
current modeled ExaA sources and2crystal chains;5historical single-Ca exclusions
remain separate recovery provenance, not additional current observations.
Authoritative source comparison:
`crystal_saved_comparison_v1/execution/REPORT.md` and RESULT.json. Ten current
modeled sources give6Ca/4inconclusive; crystalA givesCa/Bunavailable. Five historical
single-Ca failures are separate, not additional current observations. The
source/packing hypotheses are now evaluated; no automatic further expansion,
refolding, gate relaxation, redox change, refit or production promotion follows.

Preparation/test jobs and commands are in
`workspaces/exaa_source_packing_20261014/{crystal_v2,packing_v1,matrix_review_v1}/`.
Root inspected their actual results; delayed wakes do not request repeats.
Production and all historical results remain unchanged.
