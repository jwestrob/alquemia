# Allocation-aware runtime wrapper

`scripts/metal_environment_runtime.py` wraps the unchanged renderer and adds
allocation-derived `%maxcore` only to runtime inputs. It preserves the original
schema, template path/hash, parser method identity and MPI provenance. The final
runtime hash and sidecar include the memory policy and both implementation hashes.

Pinned executor layout:

- Copy wrapper as `implementation/render_orca_runtime_input.py`.
- Copy unchanged renderer as `implementation/_base_render_orca_runtime_input.py`.
- Pin both files in the implementation inventory; the existing runner's renderer
  hash then identifies the actual wrapper it imports.
- Executor worker count must match `METAL_ENV_WORKERS` (default 4). Ranks are
  supplied through the unchanged `nprocs` argument (planned default 16).
- Call `allocation_policy()` first if adapting concurrency to the allocation.
  Rendering rejects an overcommitted executor layout rather than silently
  changing the executor's worker count.

Memory comes only from positive `SLURM_MEM_PER_NODE`, or from
`SLURM_MEM_PER_CPU * SLURM_CPUS_ON_NODE`. No `/proc`/whole-node guesses. Explicit
local Slurm task slots and allocated CPUs constrain concurrency. Reserve 25%;
convert MiB to decimal MB; divide by all concurrently active ranks; round down
to a multiple of 100 MB. `%maxcore` is ORCA's working-memory setting, not a claim
that actual total process memory is strictly bounded by it.

At 64 ranks, `--mem=256000` MiB yields 3145.728 decimal MB/rank before rounding,
so `%maxcore 3100`. **256 GiB instead means 262144 MiB**, yielding 3221.225472
decimal MB/rank and `%maxcore 3200`. These are distinct allocations.

Scientific templates must omit `%maxcore`; existing declarations reject rather
than silently overriding a potentially intentional policy. Tests adapt a copy
of an archived real 1H4I input by removing only its resource line; method,
geometry, charge, multiplicity and embedding references remain untouched.

```bash
python -m unittest discover -s tests -p test_metal_environment_runtime.py -v
```

Four tests passed, zero skips; no molecular run or new allocation.
