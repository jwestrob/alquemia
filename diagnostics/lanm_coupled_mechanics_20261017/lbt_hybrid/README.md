# Local LBT hybrid adapter: root execution handoff

Ready for Slurm preflight; no new molecular results are implied. Source received
AST validation and independent static review only. `CONFIG.json` was produced by
the metadata-only `bind` stage. All numerical stages remain root-owned.

## Commands inside the allocation

From the repository root, with a fresh output directory:

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python \
  diagnostics/lanm_coupled_mechanics_20261017/lbt_hybrid/EVALUATOR.py \
  --stage preflight \
  --config diagnostics/lanm_coupled_mechanics_20261017/lbt_hybrid/CONFIG.json \
  --output /absolute/new/workspace/preflight
```

The preflight uses no model forward, OpenMM Context or new solvent/GFN call.
It checks existing real artifacts and writes the exact actual direction IDs and
vectors to `evaluator/FD_DIRECTIONS.json`. Review its `RESULT.json` before the
separate numerical scout. Preflight needs only one Slurm CPU; the interpreter
remains the pinned MACE venv so software joins match the numerical stage.

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python \
  diagnostics/lanm_coupled_mechanics_20261017/lbt_hybrid/EVALUATOR.py \
  --stage scout \
  --config diagnostics/lanm_coupled_mechanics_20261017/lbt_hybrid/CONFIG.json \
  --output /absolute/new/workspace/scout
```

Use one GPU and12 Slurm CPUs for `scout`. No memory or time override. Nine
composites are declared, with the central ALPB exactly reused. Two GFN calls
per composite run concurrently on one thread each after the MACE call. Native
MM runs in its separately pinned Python subprocess within the same allocation.
Every component uses current physical coordinates. Root attaches its existing
completion watcher; this adapter does not create background jobs.

For `origins`, copy the same config into a new workspace and add only
`scout_result: {path, sha256}` pointing to the actual admitted scout RESULT.
Run `--stage origins --config ... --output ...` with one GPU/12 CPUs. It reuses
the admitted central local point and all four compatible original global-ALPB
cells, evaluating at most seven new composite origins. Full8-cell availability
and within-region source responses are recorded. Root reviews measured scout
costs before submitting this stage.

The same evaluated model is callable through
`HybridEvaluator(cfg, fresh_output).evaluate(source, region, metal, positions)`.
Positions are all256 real atoms, with physical target index255. Returned local
MACE/GFN core index0 is mapped back to255; caps have no independent coordinates.
`difference(new, old)` returns stable native-component work in kcal/mol.
Returned provenance pins the complete cell and every component.

## Scope and remaining gates

Baseline and expanded partitions have172/210 electronic atoms and84/46 external
charges respectively. Neither has the old48-atom boundary. The four same-source
origin pairs use one consistent ledger and the qualified full corrected La/Dy
lineage. The La baseline scout qualifies only its declared directions; Dy and
expanded-region derivatives still need separate finite checks before optimization.

Independent review caught and fixed a provenance-only issue before submission:
displaced components now cite current PHYSICAL.json/core.xyz, with original
preparation retained explicitly. Original DFT ECP metadata is kept as unexecuted
preparation metadata, separate from the actual spin-independent GFN2 state.

No optimization, new DFT, source mixing, changed parameter, production change or
experimental preference claim is implemented. See `SPEC.md` for the frozen
interaction accounting, dimensions, tolerances and complete finite denominator.
