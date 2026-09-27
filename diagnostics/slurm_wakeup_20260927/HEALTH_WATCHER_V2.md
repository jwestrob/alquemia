# Live SCF watcher v2 — early DIIS pathology

Version `native_scf_health_v2` adds an alert before a problematic calculation reaches TRAH. v1policy/evidence files remain unchanged. The updated script retains all v1notification behavior; newv2receipts are required so existingv1identity is not silently migrated. Running oldwatchers remain as launched.

## Frozen rule

For the latest12DIIS iterations, require every DIISErr>0.01, RMS density change>0.1 and maximum density change>1; additionally require energy span>1Eh and latestDIISerror not less than half the first. This describes gross self-consistency failure/nonprogress, not a biological classification or an automatic endpoint rejection. It is alert-only, with no walltime cap or automatic cancellation. Values are ORCA's printed diagnostics.

A DIIS alert is permitted only while that table is the current solver phase. Historical DIIS rows cannot trigger it after enteringTRAH orSOSCF. Actual latest rows are retained in the receipt for root's decision.

## Actual validation

Six real-log tests pass. Prefixes throughDIISiteration16 of both oldDy attempts trigger the new alert; prefixes of both oldLa successes at that sameiteration do not. CompleteDy logs transition toTRAH and therefore no longer label historical DIIS behavior as the current phase. All prior actualfailure/success/progress tests pass. These consumed logs are development fixtures, not independent validation of the thresholds.

Exact policy/timestamp/log hashes/parsedprefixes: HEALTH_POLICY_V2_AND_TESTS.json. No new molecular calculations or Slurmjobs. Root must arm the updated script on explicitlyowned recovery jobs as documented inHEALTH_WATCHER.md, using a freshv2receipt.
