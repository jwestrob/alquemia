# Root handoff: GPU parity adapter

The adapter and saved-only tests are prepared. **No native GPU result exists.**
Root owns completion of the isolated build, executable binding and submission.

## Files

- `SPEC.md`: unchanged numerical gates and finite run inventory.
- `RUN.py`: runtime admission, native output parsing, parity and failure records.
- `CONFIG_TEMPLATE.json`: pinned scientific inputs; final build/executable fields
  deliberately unbound. Preserve this file and write `CONFIG_BOUND.json` separately.
- `RUN.sbatch`: one H200 / four CPUs / `standard_expand`, no memory/time override.
- `TEST_RUN.py`: actual saved CPU-output replay and refusal-path tests; zero native calls.

The matching workspace preserves `draft_before_resume_v1/`,
`saved_preflight_v1/`, and `validation_v1/`.

## Binding contract

After a successful reviewed build, root copies the template to a new config and
sets `binding_status` to `root_bound`. Each pin is exactly
`{"path": "absolute existing path", "sha256": "full file SHA256"}`.

Fill `executables.all_tests` and `executables.tinker9` with the actual binaries.
Fill `completed_build.cmake_cache`, `completed_build.compile_commands`, and
`completed_build.receipt`. Put actual logs, build decisions and compatibility
patch records in `completed_build.additional_provenance` as pins.

For each source repository (`.`, `tinker`, `ext/ext/catch2v3`), fill
`source.working_tree_diff_sha256` with SHA256 of the exact bytes returned by:

```text
git -C REPOSITORY diff --binary --no-ext-diff --ignore-submodules=all HEAD --
```

Submodule HEADs are checked separately. The empty diff has the usual SHA256 of
empty bytes. Known Catch2 header patches are preserved and pinned, not ignored.
If any scientific source artifact named in the template changes, stop for
explicit review instead of silently regenerating its pin.

The passed-build receipt can be a separate root binding record pointing to the
actual build evidence. The adapter requires this minimal structure:

```json
{
  "status": "build_passed",
  "source_commit": "44bcd7c898e68f6827c112df79b8361e443730d2",
  "source_working_tree_diff_sha256": {
    ".": "ROOT_DIFF_SHA256",
    "tinker": "TINKER_DIFF_SHA256",
    "ext/ext/catch2v3": "CATCH2_DIFF_SHA256"
  },
  "executables": {
    "all_tests": {"path": "/actual/all.tests", "sha256": "ACTUAL_SHA256"},
    "tinker9": {"path": "/actual/tinker9", "sha256": "ACTUAL_SHA256"}
  }
}
```

The build record's binary/diff pins must equal the bound configuration. The
adapter checks the cache settings CUDA, GPU host mode off, double precision,
and custom `COMPUTE_CAPABILITY=90`. Root also reviews actual compiler commands
and loader requirements; this preparation does not certify build success.

Only saved inspection is permitted on login. `RUN.py --saved-preflight-only`
verifies science inputs and source artifacts without loading either binary.
It is not completed-build or GPU qualification. Any bound runtime execution
requires the actual H200 Slurm allocation. Native errors, suite failures,
missing arrays, failed numeric gates, or changed pins stop the attempt.

## Successful execution scope

A passing result has all four native POLPAIR cases and all four control rows,
with no missing values. It qualifies this engine/build on these tests only.
It does not authorize parameter promotion, protein dynamics, a new damping
model, or a claim of physical accuracy. Root retains those scientific decisions.
