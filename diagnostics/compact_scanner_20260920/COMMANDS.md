# Prepared-input composite scanner operations

The successful four-site pilot is complete: job1203319, `pilot_v2`. Read and
validate it without rerunning science. Production defaults remain unchanged.

```bash
SCANNER_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
SCANNER_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
SCANNER_GPU_PY="$SCANNER_ROOT/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$SCANNER_ROOT"
"$SCANNER_PY" scripts/compact_solvation_scanner.py validate \
 --manifest workspaces/compact_scanner_20260920/pilot_v2/manifest.json
cat diagnostics/compact_scanner_20260920/REPORT.md
"$SCANNER_PY" -m unittest discover -s tests -p 'test_compact_solvation_scanner.py' -v
sacct -n -P -j 1203299,1203319 \
 --format=JobID,State,ElapsedRaw,AllocCPUS,CPUTimeRAW,TotalCPU,MaxRSS,ReqMem,AllocTRES,NodeList
```

Actual result: `workspaces/compact_scanner_20260920/pilot_v2/result_1203319.json`.
Reviewed collector replay: same directory, `reviewed_collection_v2.json`.
Both exist and are immutable. Collection needs a new explicit output filename;
it does not calculate missing components or substitute baseline scores.

## Explicit prepared pair: preparation only

This real1H4I source example contains **no old MACE energies**. The following
command builds a new finite manifest but runs no scientific executable. Its
output directory must not already exist. To use another site, supply a different
pair JSON path through `--pairs`; no source-code editing is needed.

```bash
"$SCANNER_PY" scripts/compact_solvation_scanner.py prepare-pairs \
 --pairs diagnostics/compact_scanner_20260920/PREPARED_PAIR_1H4I.json \
 --comparison workspaces/compact_solvation_20260920/full_v1/comparison_v1.json \
 --agreement diagnostics/compact_scanner_20260920/PLAN.md \
 --output workspaces/compact_scanner_20260920/example_1H4I_explicit_pairs_v1 \
 --cpu-python "$SCANNER_PY" --gpu-python "$SCANNER_GPU_PY"
"$SCANNER_PY" scripts/compact_solvation_scanner.py validate \
 --manifest workspaces/compact_scanner_20260920/example_1H4I_explicit_pairs_v1/manifest.json
```

Pair JSON records `model`, `software`, and `cases`. Each case has identity,
biological group, role, label scope, expected class (null for predictions), and
`representations.context.{protocol_id,preparation,endpoints}`. Each Ca/La endpoint
has a hash-pinned `xyz`, `charge`, `multiplicity`, `metal_index`. Context provenance
must come from a supported complete-context preparation and preserve source
geometry/water inventory. Archived energies/receipts are optional. Supported
PQQ bands require the canonical PQQ complete-context protocol **and** PQQ scope;
generic preparations receive raw contrasts, not PQQ calls. Use
`PQQ_prediction` for unlabeled compatible PQQ. The operation never fits bands to
the incoming cases. Arbitrary raw XYZ/folds are not automatically compatible.

## Fresh execution, when intentionally requested

This example would repeat the real1H4I source pair; it was **not** executed as an
additional scientific test. The four-site measurement above is already complete.
The corrected batch reserves1H200/32CPUs/200000MiB with32×1CPU Slurm tasks and
an explicit1×32CPU MACE step. Each site uses2MACE+4nativeGFN2 calls.

```bash
sbatch --parsable \
 --output="$SCANNER_ROOT/workspaces/compact_scanner_20260920/example_1H4I_explicit_pairs_v1/slurm_%j.out" \
 --error="$SCANNER_ROOT/workspaces/compact_scanner_20260920/example_1H4I_explicit_pairs_v1/slurm_%j.err" \
 diagnostics/compact_scanner_20260920/run_recovery.sbatch \
 "$SCANNER_ROOT/workspaces/compact_scanner_20260920/example_1H4I_explicit_pairs_v1/manifest.json"
```

The batch runs its immutable implementation snapshot, executes fresh endpoints,
and writes the score/partial result automatically. MACE uses the pinned venv
invocation separately from its symlink target; ORCA uses the established MPI
runner. A partial failure stays null; there is no fallback and no automatic
scientific retry. Re-collect a later partial result with `collect --manifest`
and a new explicit `--output` filename.

## Historical failed layout

`SUBMISSION.json` and `run.sbatch` document1203299/`pilot_v1`. They are preserved
for provenance and are **not current execution instructions**. The1×32CPU layout
exposed only one MPI slot; all16GFN2 attempts stopped beforeSCF. Recovery and its
complete actual cost are recorded in `RECOVERY.md`, `RECOVERY_SUBMISSION.json`
and `SACCT.txt`. The current executable batch is `run_recovery.sbatch`.
