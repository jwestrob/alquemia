# Exact canonical source audit

All28 sources pass:25 designated calibration proteins,1H4I/4MAE structural
transfer and1KB0 external class transfer. No new scientific output was created
by the audit. It reuses the existing fixed-core score_task verifier for54DFT
endpoints and the existing manifested execution verifier for the two1KB0 legs.
It checks normal/converged receipts, exact inputs/outputs, released baseline
contrasts, paired coordinates, fragment charge closure, singlet states, PQQ(3−)
and the frozen dry policy. PQQ cores are unchanged; no peptide repair is applied.

Authoritative inventory:
workspaces/mace_canonical_20260916/audit_v2/inventory.json
SHA256 f51a8641fa50a34087fca9ba76d0e5ca960ac85c1ea5da65a3fa19bed9d5cb21.
Its implementation copy is immutable. audit_v1 contains the same scientific
inventory with its original implementation preserved separately; v2 fixes source
pinning so later edits cannot invalidate that preparation record.

Two initial audit reads exposed legacy serialization details (tasks are a list,
and some input paths are relative to their source manifest). The adapter now
resolves those exact saved paths; no scientific input or output was regenerated.
The initial MACE runner package omitted the existing curvature/gradient parser
helpers needed by the receipt reader. An isolated-environment dry-run caught it
before any Slurm submission. mace_v2 adds those helper files; all task inputs
and model settings are unchanged from the preserved mace_v1 preparations.
Failed preparation/dry-run CPU costs were not independently profiled; no GPU or
DFT calculation ran in those failures.

## Evidence grouping

PDB DBREF records in the pinned original files give chainA1H4I→P16027 and
chainA4MAE→I0JWN7. Both accessions occur in the25 calibration panel. Group both
as structural replicates of calibration proteins, without claiming exact
sequence identity from DBREF alone.1KB0 maps to the existing Q46444 control.
All three transfer cases are retrospective; no new blind case is consumed.
Homologous calibration entries are not assumed independent biological samples.

The1KB0 construction manifest says orca_executed:false because it predates the
completed baseline run. Its actual26-endpoint campaign collection1199299 and
verified per-leg receipts supply the completed1KB0 result. Do not interpret the
old construction-time flag as the current job state. Original bands and energies
remain intact in their released records.

Medium reuses only the four exact archived-core1H4I/4MAE endpoints from
MACE1200717 and GB1200719. Coordinates, charge, spin, checkpoint, software and
all effective model/solvent settings match; the preparation-policy metadata name
is the sole model-record difference. Large has no matching local cache and
runs all56logical endpoints. No DFT cache is accepted as MACE or GB evidence.

## Execution and tests at this checkpoint

MACE medium1200776 and large1200779 use oneA5000,16CPU and64474MiB each, existing
finite runner and per-task receipts. Source tasks52/56 plus4medium reused rows;
zero newDFT. Solver preparation follows only complete actual MACE collections.

Five real source/runner/calibration tests passed8.933s. One additional actual
MACE+GB archive test passed1.348s, reproducing1H4I/4MAE component contrasts with
correct Ca−La sign and single eV conversion. These parser/algebra checks are
separate from the actual scientific model calculations. No synthetic successful
scientific output is substituted for an executable.

## Completed checkpoint

All four jobs completed; solvent IDs 1200781/1200782. Both checkpoints failed
calibration separation; see REPORT.md and result.json. Seven canonical tests
now pass, including the actual failed-calibration regression. Baseline runner
regression: 10 passed, two explicit environment-dependent skips. No scientific
failure or new DFT. Email delivery of the original test was confirmed directly
by Jacob; a requested progress email was submitted successfully. Private
notification receipts remain under workspaces, not in this diagnostic.
