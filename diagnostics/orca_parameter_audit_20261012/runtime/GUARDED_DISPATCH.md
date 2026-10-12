# Guarded native-GFN2 dispatch

`scripts/guarded_native_gfn2_dispatch.py` is a new callable execution path. It uses
the shared `native_gfn2_parameters` guard for staging, preflight and postflight.
It has no fallback to native defaults. Existing scientific state, convergence
and population checks run separately from parameter admission.

Prepare a **fresh** task directory with `endpoint.inp`, its XYZ file and an
`expected.json` record containing the independently derived charge,
multiplicity, valence-electron count, atom count and basis-function count. The
input must request serial `nprocs 1`, `NoAutostart`, and the exact local import:

```text
%method
 XTBParamFile "control.xtb.json"
 ReadXTBParam true
 WriteXTBParam false
end
```

For example, after preparing those actual inputs, use the following from the
project root. Replace `NEW_TASK` with the fresh directory and supply the real
state/medium/source; the values below describe the existing calcium vacuum
control, not a generic protein state:

```bash
python scripts/guarded_native_gfn2_dispatch.py \
  --input NEW_TASK/endpoint.inp --xyz NEW_TASK/core.xyz \
  --expected NEW_TASK/expected.json \
  --parameters params/native_gfn2_parameters_v1/corrected_la.xtb.json \
  --policy params/native_gfn2_parameters_v1/la_shell_three_arrays_v1.json \
  --approval diagnostics/orca_parameter_audit_20261012/parameter_guard/runtime_registry_ca_la_corrected_20261012.json \
  --runtime-binding workspaces/orca_parameter_audit_20261012/runtime_qualification_v1/RUNTIME_BINDING.json \
  --observed-runtime workspaces/orca_parameter_audit_20261012/runtime_qualification_v1/OBSERVED_RUNTIME_FINAL.json \
  --protocol-id YOUR_DECLARED_PROTOCOL \
  --source YOUR_REAL_SOURCE_MANIFEST \
  --charge -6 --multiplicity 1 --medium vacuum
```

The default performs staging and preflight only. An explicitly requested
`--execute` in the initial invocation runs one bounded serial endpoint and
postflight, then records independent scientific admission. A staged dry-run
directory is deliberately not reused for execution; prepare a fresh directory.
The current research registry supports the tested corrected Ca/La profile in
vacuum and `ALPB(Water)`, with multiplicity-one bookkeeping. Unsupported profile,
element, medium, runtime, import bytes or state combinations fail closed.

The dispatcher includes the parameter manifest/cache identity, full source and
implementation pins, preflight and postflight in its receipt. A parameter check
alone cannot accept scientific output. The receipt records `null` scientific
acceptance for a dry run. The finite original 60-control runner remains
unchanged; this new path does not retroactively claim those controls were
guard-dispatched.

## Actual validation

`GUARDED_DISPATCH_VALIDATION_V2.json` pins four successful stage/preflights using
the real corrected Ca/La inputs in both media, rejection of an original erroneous
full file and a mismatched expected state, and an injected parser failure. The
failure check verifies that the owned executor stop flag is set, original signal
handlers are restored, and a failed receipt survives. Its executor/postflight
are explicitly mocked; no simulated scientific value is admitted. The earlier
v1 dry-run records are retained as superseded evidence from before the reviewed
cleanup repair. **Zero additional molecular calls** were made. The new dispatch
execution branch is implemented but has not itself run a molecular endpoint;
the original 60 controls and the shared guard's independent tests retain their
separate evidence scopes.

This opt-in path does not change historical or concurrent production callers.
The campaign's shared integration owner must route any new native work through
the guard rather than claiming every pre-existing entry point is already gated.
