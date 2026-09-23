# Reusable three-source experimental command

**The same experimental scorer now accepts explicit complete triples without a
two-protein whitelist.** New plans use execution-interface v2; all scientific
settings, the actual envelope reference and existing warm/scalar executors remain
unchanged. The completed v1 integration remains historical and replayable. No new
molecular execution, source reconstruction, protonation or score was performed.

| Actual preflight | Sources/groups | Origin MACE/scalar inputs | Maximum later searches/cross-MACE/scalar total |
|---|---:|---:|---:|
|07ab triple|3/1|6/12|6/6/36|
|8344 triple|3/1|6/12|6/6/36|
|Both triples|6/2|12/24|12/12/72|

These are alternative unexecuted plans, not three submitted experiments. All
counts derive from the explicit source list. Any incomplete group is rejected;
missing result members cannot yield a reduced-member median or baseline fallback.
Arbitrary safe protein IDs are metadata, not biological labels. Actual source
sequence/site/assembly and complete chemical-state validation remain mandatory.

The code change is restricted to `pqq_three_source_envelope_execution.py`: select
v1/v2 interface behavior from the plan, derive counts, and keep complete-group
guards. No shared workflow or optimizer was replaced. The operative optimizer
tolerance is1e-8 on the Hartree-equivalent scaled objective; checkpoint, source
coordinates, charge/spin, maps, scalar recipes and all settings replay exactly
against the successful v1 integration. There is no altered Hamiltonian or new
calibration. Molecular cache keys remain scoped to the actual pinned request,
plan, method, state and coordinates through the existing executors.

Six focused API tests pass in28.726s and seven completed-integration regressions
pass in13.300s, zero skips. Tests cover separate/joint plans, identical actual
coordinates/input text, an explicitly named metadata alias of a real prepared
protein, rejection of two genuinely different source sequences, corrupt-state and
duplicate-group rejection, strict missing-member handling, and the old v1 matrix
and decision replay. No successful scientific output was fabricated. The actual
CLI also regenerated a report from the unchanged v1 result without recalculation.

Reference and benchmark qualifications are unchanged: canonical25 and3crystal
calls are retained, but A0A3Ca3 is an inconclusive probe; full100-triple transfer
is a separate qualification task. This interface adds usability, not accuracy
evidence or a production recommendation. Historical DFT, expression/provenance
joins, released defaults and the earlier v1 outputs remain untouched.

[Commands](COMMANDS.md) provide actual request/preparation, preflight, explicit
finite execution and collection/report paths. [Artifact pins](ARTIFACTS.json)
identify all three unexecuted plans and the unchanged historical result. No new
job was submitted, and no GPU/solver allocation was consumed in this task.
