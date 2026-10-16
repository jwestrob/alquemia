# Dy GFN2 parameter and runtime readiness

16 October 2026. Bounded read-only audit; no parameter edits, molecular calls,
package imports, or jobs. The eight-cell LBT3 MACE-vacuum scout remains separate.

## Finding

**The current corrected-La native profile does not qualify Dy.** Its Dy shell
labels are correct, but three coefficient arrays disagree with the author data
when indexed by those exact labels. Both files declare `[5d, 6s, 6p]`:

| Dy field | Current native full profile | Author coefficients in declared shell order |
|---|---|---|
| `kcn` | `[0.1449494, -0.0197233, 0.012143]` | `[0.012143, 0.1449494, -0.0197233]` |
| `shpoly` | `[-0.05567124, 0.41914074, -0.36367686]` | `[-0.36367686, -0.05567124, 0.41914074]` |
| `lgam` | `[1, 0.8846158, 0.5075342]` | `[0.5075342, 1, 0.8846158]` |
| `refocc` | `[1, 1, 1]` | `[1, 1, 1]` |

Exact sources: `params/native_gfn2_parameters_v1/corrected_la.xtb.json`
(Dy begins near line 3317), `author_gfn2.toml` (Dy near line 1282), and
`diagnostics/orca_parameter_audit_20261012/parameter_guard/PARAMETER_AUDIT.json`.
`scripts/native_gfn2_parameters.py::audit_parameters` audits these four fields
by full shell identity with the existing `atol=rtol=1e-12` parameter-roundoff
tolerance. This tolerance is not an energy or molecular-acceptance tolerance.

The current full-profile digest is
`a4f83ca72d58d6e80f634e61bb8e99b7af0f26c38cc4a9295e85215ecf72d05c`.
Original export and author TOML digests are respectively
`f3dac846a59238947d35c23ba5bdb40c19a34482fe9860e8e11bf9c00d719107` and
`76239337c5f7f3285b0b4db06d607490d14cb38a9db90334f212e370f328082a`.
A proposed La+Dy profile must preserve the existing La correction and change
only the three explicitly named Dy arrays; declare both elements in its policy.
No such new profile was created by this review.

## Actual current admission boundary

The machine-readable registry
`diagnostics/orca_parameter_audit_20261012/parameter_guard/runtime_registry_ca_la_corrected_20261012.json`
admits **C, Ca, H, La, N, O**, vacuum and ALPB(Water), multiplicity-1 bookkeeping,
under one exact parameter-file/runtime identity. Its prose scope and REPORT say
CHNOS, but the actual `qualified_elements` list omits S; follow the list.
Neither Dy nor a different full-file digest is admitted.

Native GFN2's original spin-independent, f-in-core bookkeeping is distinct from
the MACE Dy sextet input. Preserve the separate physical Dy(III) hypothesis;
do not transfer multiplicity 6 into this native GFN2 convention or infer
physical spin energetics from successful closed-shell bookkeeping.

## Reusable implementation and actual controls

| Purpose | Exact existing path / entry point | Limit |
|---|---|---|
| Native full-file import guard | `scripts/native_gfn2_parameters.py`: `build_binding`, `stage`, `preflight`, `postflight` | New profile uses diagnostic qualification before any research registry admission |
| Native executor pattern | `workspaces/orca_parameter_audit_20261012/runtime_qualification_v1/RUNTIME_BINDING.json`; guard API above | The existing `guarded_native_gfn2_dispatch.py` hardcodes research approval and cannot run fresh Dy qualification unchanged; use direct diagnostic staging |
| Native observed runtime | `diagnostics/orca_parameter_audit_20261012/parameter_guard/OBSERVED_RUNTIME_SNAPSHOT.json` | Includes actual helper processes; no new runtime inferred from executable name |
| Native completed controls | `diagnostics/orca_parameter_audit_20261012/runtime/RESULTS.json`; `runtime_review/RESULT.v2.json` | 60 completed / 52 admitted; retained default/roundtrip density failures and Zn-water charge-parity failure |
| Standalone executable | `/home/jwestrob/miniconda3/envs/lanm_bench/bin/xtb` | Installed xTB 6.7.1; pin actual bytes and runtime |
| Standalone parameters | `/home/jwestrob/miniconda3/envs/lanm_bench/share/xtb/param_gfn2-xtb.txt` | Use original standalone author parameters; do not import ORCA JSON here |
| Standalone task pattern | `scripts/standalone_xtb.py::command`, `execute_cell`, `parse_task`, `collect` | Existing finite accuracy/state/solvent checks; original script is PQQ-scoped |
| Shared backend runtime pins | `workspaces/orca_parameter_audit_20261012/context_gfn2_v1/RUNTIME.json` | Includes standalone and independent backend identities |
| Independent tblite | `workspaces/orca_parameter_audit_20261012/tblite_0_7_0/venv/bin/python`; `scripts/context_factorial_gfn2.py` | Its present element/valence dictionaries omit Dy; not an already-ready Dy adapter |

Actual standalone-versus-tblite component precedent is
`workspaces/orca_parameter_audit_20261012/context_gfn2_component_v2/FINAL_COLLECTION.json`:
40/40 Ca/La solvent components and atomic populations pass. Absolute total
energies failed the original parity comparison, traced to dispersion differences.
Do not silently equate absolute totals or transplant a total-energy threshold
without the declared decomposition. These receipts supply comparison patterns;
none establishes Dy parity.

### The documented dispersion mismatch is a different backend pair

The context-factorial discrepancy is specifically **standalone xTB 6.7.1 versus
tblite 0.7.0**. Its four-call causal diagnostic changed only three-body D4 `s9`
and located the geometry-dependent residual in that non-self-consistent term;
it did not test ORCA native versus standalone or establish a universal
element-reference offset. See `context_GFN2/SCOUT_REPORT.md` and
`DISPERSION_DIAGNOSTIC_PLAN.md` under the October 12 diagnostic.

There is direct corrected-native/standalone precedent on a larger 222-atom
fixed state. `runtime/RESULTS.json` → `La_Ca_upstream_comparisons` records
corrected La vacuum/water residuals of −2.1278771100696758e−7 and
−4.5639524159923894e−8 Eh, with maximum atomic-charge residuals
3.672771312412948e−6 and 4.37719767243383e−6 e. Both pass the unchanged
1e−6 Eh / 1e−5 e gates; Ca controls pass as well. These are recorded values,
not newly calculated here. Exact original interpretation and upstream identity:
`/groups/banfield/users/jwestrob/obsidian-vault/projects/lanthanide-binding/orca_parameter_audit_2026-10-12/la_followup/native_reciprocal_controls/REPORT.md`
and its `endpoint_comparison.tsv`; the enclosing `la_followup/REPORT.md`
explicitly identifies the independent author-parameter standalone xTB route.

Therefore retain the declared absolute gates for the new 48-atom native/xTB
pair. No known mismatch justifies pre-excusing a failure as dispersion. If an
absolute gate fails, preserve that failure and report atomic populations,
same-geometry ALPB-minus-vacuum agreement, and available component diagnostics
separately. Do not fit an offset, silently change the qualification target,
or promote component agreement to absolute parity.

Reusable parsing pieces:

- `diagnostics/orca_parameter_audit_20261012/runtime_review/REVIEW.py::parse_raw`
  reads native property `Single_Point_Data.FinalEnergy`, full-precision atomic
  charges, stdout/state and all three printed SCF criteria. Use the property
  populations rather than four-decimal printed charges for a 1e−5-e gate.
- `workspaces/orca_parameter_audit_20261012/runtime_qualification_v1/source/run_native_controls.py::parse_output`
  supplies native state, termination and density checks; its four-decimal Mayer
  charges alone are insufficient for that population comparison.
- `scripts/standalone_xtb.py::parse_task` checks actual GFN2, charge, unpaired
  count, effective electrons, restart/settings, JSON energy and ALPB reference.
- `diagnostics/pede_solvent_response_20261013/solvent_audit/DECOMPOSE.py::scalar`
  and `stdout_record` parse standalone SCC/repulsion/restraint, dispersion,
  `Gsolv/Gelec/Gsasa/Ghb/Gshift` and closure. Its inventory wrapper is
  PedE-specific. Printed total dispersion is not an isolated three-body term.
- The inspected native ORCA outputs/property file do not expose an equivalent
  component-resolved three-body D4 scalar through these parsers. Thus a future
  discrepancy cannot be causally assigned to that term merely by subtracting
  unmatched printed decompositions. Any additional causal comparison requires
  a separately declared, bounded experiment.

The existing standalone command requests GFN2, exact charge, `--uhf 0`,
electronic temperature 300 K, maximum 500 SCC iterations, explicit accuracy,
`--norestart`, isolated parameter path and JSON populations; ALPB adds water
with the declared solution reference convention. Its solvent parameter
temperature is 298.15 K, separate from electronic smearing temperature.
Use the actual declared accuracy and acceptance from the next root-owned plan,
not whichever result later appears most favorable.

The reusable serial native command pattern is
`[ORCA, 'endpoint.inp', '--bind-to none']`; exact runtime/environment comes from
the pinned binding. New qualification stages the imported profile through
`stage(..., mode='diagnostic', variant='corrected')`, with actual preflight and
postflight receipts. `scripts/context_factorial_gfn2.py::clean_env` and
`execute_point` provide the one-thread standalone command/environment pattern;
its broader preparation function is not Dy-ready.

## Native fixed-charge ALPB source API

Reusable source root:
`workspaces/pede_solvent_response_20261013/fixed_charge_replay_v3/execution_v1/sources/xtb/`.
The pinned source is xTB 6.7.1, commit
`26b28010e805f7d1aeeef39813feb473e69cc4be`.

- `src/api/calculator.f90`: GFN2 loading and solvent API.
- `src/xtb/gfn2.f90`, `src/solv/model.f90`, `src/solv/gbsa.f90`.
- Static library: same execution root, `build/xtb/libxtb.a`.
- Existing wrapper: `workspaces/pede_solvent_response_20261013/fixed_charge_replay_v3/implementation/REPLAY.f90`;
  correct `mctc_io_convert` input conversion, `newBornModel`,
  `born%getEnergyParts`. This wrapper hardcodes its 192-atom scope.
- Actual component controls:
  `diagnostics/pede_solvent_response_20261013/fixed_charge_replay/EVIDENCE.json`
  and linked `FINAL_COLLECTION.json`: eight evaluations / sixteen component
  checks on Ca/La. They do not supply qualified Dy vacuum charges or Dy ALPB parity.

## Smallest next qualification, separate from the MACE scout

Root's intended finite matrix is **La/Dy × vacuum/ALPB × corrected-native/standalone**:
eight fresh fixed evaluations using the exact 48-atom 7CCO matched-core pair at
`workspaces/lanm_coupled_response_20261016/lbt3_preparation_v1/execution/7CCO/{matched_core_La,matched_core_Dy}.xyz`.
Both representations here have charge −1; physical Dy multiplicity and GFN2's
effective bookkeeping must remain separate records. No geometry change, spin
substitution, failed-SCF retry, old whole-protein GFN2 rerun, or new solvent score
for the LBT3 matrix follows from this readiness note.

Use one new pinned La+Dy full native profile; the La calls directly establish
that importing that full file preserves its unaffected La behavior. The Dy
calls test actual parameter import, solver/state/population acceptance and
independent-backend behavior in both media. Root will declare exact numerical
tolerances from the existing receipts before execution. Only accepted results
can support a new runtime registry or subsequent solvent use. Do not reuse the
old uncorrected Dy native solvent components or describe an unrun correction as
a completed composite score.

## Source-only staging adapter

`diagnostics/lanm_coupled_response_20261016/DY_QUALIFY.py` is ready for the
root-owned plan and executor; SHA256
`71c16ea9d4f7978c2c4e3f6cc855d8522b7ebc0daecbf2a5fb98eb2d9fa5c8a8`.
It stages the exact eight cells, task-local full La+Dy profile, diagnostic
parameter policies, commands, state records and source/runtime pins. It takes
`--config ROOT_CONFIG.json --output FRESH_DIRECTORY`; the required config keys
are documented in its module header. It does not execute engines or qualify a
runtime. The only static review blocker was a file-location-dependent repository
root; the final source now pins the absolute repository path for snapshot use.
Neither adapter execution nor numerical work was performed in this review.
Staging is sequential; concurrent native guard checks require a constant pinned
process environment, rather than changing it independently in worker threads.
