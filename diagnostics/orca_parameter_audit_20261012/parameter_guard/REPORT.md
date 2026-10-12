# Native GFN2 parameter guard — 12 October 2026

## Result

`scripts/native_gfn2_parameters.py` provides an opt-in, fail-closed parameter
contract for **new research tasks**. It does not change a production runner,
historical result, calibration or default. The full pinned La correction changes
exactly nine scalar values in `La.kcn`, `La.shpoly` and `La.lgam`; all other JSON
values remain exactly equal to the original export.

The helper audits `kcn`, `shpoly`, `lgam` and `refocc` by each declared full shell
identity, such as `5d`, against the author TOML. It never applies a global array
permutation. The predeclared comparison tolerance is `atol=rtol=1e-12`, solely
for author/export decimal roundoff. Unrelated values have **no tolerance**.

After the La-only correction, 46 of 86 exported elements match these four author
fields and 40 retain a mismatch. This is a parameter-schema observation, not
runtime qualification of 46 elements or causal proof for all 40 others. In
particular, Dy and Cu are **not qualified** by the corrected-La profile.

## Runtime admission

The runtime registry references actual completed October 12 controls from the
independent runtime task, not fabricated test outputs. It binds the complete
parameter-file digest, ORCA driver, native helper and shared-library hashes,
recorded environment and original-GFN2 state convention. Each control must retain
successful execution, normal termination, accepted convergence/state/population
checks, exact input/output hashes and real parameter-import stdout evidence.

`runtime_registry_ca_la_corrected_20261012.json` covers CHNOS+Ca+La,
vacuum and ALPB(Water), and multiplicity-1 bookkeeping. It includes the actually
observed `orca_guess` helper and `dash` process in addition to the pre-run library
and helper inventory. `OBSERVED_RUNTIME_SNAPSHOT.json` is an immutable byte-exact
copy of the actual observation inventory; it does not depend on an agent's
still-updating collector output. The earlier La-only registry is explicitly superseded and
cannot admit new research tasks. The registry does not establish
physical spin energetics, arbitrary oxidation states, affinity, source-structure
transfer or biological prediction. New geometries reuse the same runtime registry;
they do **not** need a previously successful calculation of themselves. Their
actual input, XYZ, state, source, method and parameter pins enter the cache key.

Other elements/media/multiplicity conventions require actual controls. Default
Ca results cannot silently qualify a different imported full-file profile.
The two separately authorized corrected-profile Ca controls provide that direct
import evidence. Their optional historical `job.parameter_source` field is null
and `allowed_input_changes` remains an empty list; those metadata carryovers are
preserved. Actual full input-file pins, unchanged post-run parameter hashes and
the exact runtime import lines establish which file was used. A non-null source
descriptor, when supplied, must also match; no contradictory metadata is ignored.

## Caller contract

1. Preserve preparation/state and render the ORCA input with one explicit import
   block: `ReadXTBParam true`, `WriteXTBParam false`, and an exact task-local JSON
   filename. No native-default fallback is permitted in research mode.
2. Build the task binding with `build_binding(...)` or CLI `bind`. Include actual
   source-input hashes and the runtime binding from qualification. Supply the
   intended total charge, multiplicity and solvent medium explicitly.
3. Call `stage(...)`. It copies the **entire file byte-for-byte**, refuses existing
   parameter/manifest outputs, and emits the manifest pin and cache identity.
4. Set the recorded runtime environment and call `preflight(...)` immediately
   before execution. Retain its receipt and require success. Do not reuse a cache
   produced under a different guard identity or parameter profile.
5. Execute through the existing authorized runner, then call `postflight(...)`.
   It rehashes parameters, executable/runtime files, geometry and sources, and
   requires actual exact-file import evidence outside the echoed input.
6. **Still run the scientific collector.** Parameter postflight deliberately
   returns `scientific_result_accepted: null`; it cannot replace SCF, gradient,
   state, charge or scoring acceptance. Failure at either layer stays unavailable.

The registry is a reviewed qualification artifact. Its runtime file list must
come from the actual runtime; this helper does not discover dynamically loaded
dependencies or independently prove the list exhaustive. SHA-256 pins provide
identity/change detection, not a defense against a malicious operator rewriting
every approval and digest.

## Historical controls

`mode="diagnostic"` is mandatory for `variant="original"` and
`variant="native-default"`. Both remain explicitly unqualified. A native-default
diagnostic uses `ReadXTBParam false`, `WriteXTBParam true`, no import filename, and
must produce the pinned original full export at `<input-stem>.xtb.json`.
Corrected diagnostic tasks are also supported so a new runtime can be qualified
without circularly requiring its own completed qualification first.

## Files and verification

Pinned full data and correction policies live in
`params/native_gfn2_parameters_v1/`; `SOURCES.json` records exact snapshot origins.
The helper is standard-library Python 3.11+ and starts no subprocess or molecular
calculation. The tests use those real parameter files, actual archived/current
control outputs, and explicitly corrupted copies for rejection tests. Missing
scientific fixtures are reported as skips, never replaced by invented output.

The recorded run in `TESTS.txt` passed 38 tests in 21.469 seconds, zero skips.
That is implementation verification, not 38 new molecular evaluations. This
agent launched **zero ORCA/GFN2/MACE calculations** and used one Python thread.
Fresh scientific runtime controls belong to the separate runtime report.

## Runnable checks

```bash
python scripts/native_gfn2_parameters.py audit \
  --parameters params/native_gfn2_parameters_v1/corrected_la.xtb.json \
  --policy params/native_gfn2_parameters_v1/la_shell_three_arrays_v1.json

python scripts/native_gfn2_parameters.py directives

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -m unittest discover -s tests -p test_native_gfn2_parameters.py -v

python scripts/native_gfn2_parameters.py bind --help
python scripts/native_gfn2_parameters.py stage --help
python scripts/native_gfn2_parameters.py preflight --help
python scripts/native_gfn2_parameters.py postflight --help
```

### Python integration

```python
from native_gfn2_parameters import build_binding, pin, stage, preflight, postflight

# Arguments come from the caller's existing explicit task manifest.
binding = build_binding(
    input_path=task_input, xyz=task_xyz, parameters=full_parameters,
    runtime_binding=qualified_runtime_binding, source_inputs=source_files,
    observed_runtime=qualified_observed_runtime,
    protocol_id=protocol_id, charge=charge, multiplicity=multiplicity,
    medium=medium, parameter_filename="native_gfn2_parameters.json",
)
prepared = stage(binding, full_parameters, correction_policy,
                 approval=pin(runtime_registry), mode="research")
before = preflight(prepared["manifest"]["path"])
# Existing runner executes here; preserve its complete execution receipt.
after = postflight(prepared["manifest"]["path"], actual_stdout, actual_exit_code)
# Independent scientific collection/admission follows.
```

No new workflow system, scheduler policy or approval per geometry is introduced.
