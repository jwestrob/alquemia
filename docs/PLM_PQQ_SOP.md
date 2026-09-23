# Nikasha SOP: score PLM PQQ proteins with local accommodation

**Recommended manuscript candidate:** the explicit three-source, 4.3 Å
motion-envelope protocol. This is an opt-in command; `nikasha standard` still
uses the released static scorer. Original DFT results remain available.

This SOP starts with existing holo folds and verified PQQ-site selectors. It
does not launch folds or score the cohort automatically. The complete
[packaging record](../diagnostics/plm_pqq_delivery_20260923/REPORT.md) distinguishes
checks run here from the earlier molecular validation.

## 1. Inputs and fixed selection rule

For each protein, supply **three distinct La-conditioned structures**, the
same full sequence, residue numbering, assembly and PQQ chemical state. Use a
declared seed/sample rule before inspecting scores. For the existing PLM AF3
batch, this is all saved **seed 101, samples 0/1/2**, without filtering by
classification, motif expectation or which fold gives a preferred result.

Each source request contains the actual CIF, metal/PQQ selectors and these five
homologous roles: anchor glutamate, anchor asparagine, catalytic aspartate,
catalytic-Asp cationic partner, and the D+2 homolog. The last role need not be
acidic; its identity must be supplied correctly, not changed to match an expected
metal. Selectors must come from the existing sequence/structural mapping.
Do not copy residue numbers from an unrelated protein.

The supported chemistry is the existing **dry canonical PQQ core** and complete
polar side-chain context. PQQ is the fixed `C14H3N2O8`, charge −3 state; standard
protein protonation uses the frozen seeded pH 7 preparation. Water occupancy,
redox, alternate proton inventories, missing heavy atoms and unsupported donors
are not repaired by this scorer. Unsupported inputs remain explicit exclusions.

### Existing XoxF cohort

The original inventory contains 176 proteins; its compatible candidate manifest
contains 137. Those numbers are historical XoxF scope, **not all PLM PQQ proteins**.
New Exa/Ped/ADH-like proteins must enter through their own source/role manifests;
family names alone are neither eligibility nor biochemical labels.

The existing source builder checks the exact full sequence, actual chain/ligand
identities and all three fold files. Example using two already consumed proteins:

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
PLM_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

"$PLM_PY" diagnostics/accommodation_goal_20260920/prepare_plm_sources.py \
  --candidate-manifest workspaces/plm_xoxf_all_fixed_core_20260916/retry_1200785/prepared_batch/candidate_manifest.json \
  --release-sources diagnostics/pqq_fast_release_20260920/SOURCES.json \
  --targets PQQSEQ_07ab500e3df76b30d71c PQQSEQ_83440678cbbd658047c9 \
  --output workspaces/plm_pqq_scan_example_v1/sources
```

For another batch, change `--targets` to its explicit target list and use a new
output directory. This builder accepts the historical AF3 input contract above;
it will reject another seed, assembly or sample count. For other supported fold
sources, supply the explicit three-source JSON directly to the next step, using
the emitted files as the schema example. Do not modify code or forge historical
sample IDs to bypass a mismatch. The generic preparer checks actual sequence,
numbering, site roles and atom membership independently.

Retain each source's `source_provenance.native_AF3_input` pin when using the PLM
biological-table join below; that join verifies the full input sequence. A
different folding-input format needs an explicit verified sequence crosswalk,
not identifier-only matching.

## 2. Validate and declare the batch (no protonation or scoring)

```bash
"$PLM_PY" scripts/pqq_plm_prepare.py request \
  --sources \
    workspaces/plm_pqq_scan_example_v1/sources/PQQSEQ_07ab500e3df76b30d71c_sources.json \
    workspaces/plm_pqq_scan_example_v1/sources/PQQSEQ_83440678cbbd658047c9_sources.json \
  --protein-ids PQQSEQ_07ab500e3df76b30d71c PQQSEQ_83440678cbbd658047c9 \
  --package-path params/pqq_plm_envelope_v1.json \
  --agreement diagnostics/plm_pqq_delivery_20260923/PLAN.md \
  --output workspaces/plm_pqq_scan_example_v1/run

"$PLM_PY" scripts/pqq_plm_prepare.py dry-run \
  --request workspaces/plm_pqq_scan_example_v1/run/REQUEST.json
```

`--sources` and `--protein-ids` are paired in order. Any number of complete
triples can be declared; model execution stays warm within a batch. Preparation
is fresh by default. To reuse exact existing protonation, explicitly supply
one `--source-preparations` path per protein. Sources/configuration must match;
protein ID alone never establishes a cache hit.

Per complete triple, the planned molecular work is six origin force evaluations,
six bounded MACE searches, up to six cross-MACE evaluations, and 36 native GFN2
single points. Optimizer-internal MACE calls vary. These are finite task counts,
not a total project compute cap. There are no new DFT calls or folds.

## 3. Submit the declared batch

The wrapper performs fresh preparation, builds the complete groups' existing
execution plan, then invokes the unchanged scoring/collection wrapper. Incomplete
groups are listed in `prepared/READY.json`; other complete groups can proceed.
No two-source median or baseline substitute is created.

```bash
PLM_ROOT=$PWD
PLM_RUN="$PLM_ROOT/workspaces/plm_pqq_scan_example_v1/run"
sbatch --parsable --partition=gpu --nodelist=node-224-2t-8gpu-1 \
  --nodes=1 --ntasks=32 --cpus-per-task=1 --gres=gpu:1 --mem=200000M \
  --job-name=nikasha-plm-pqq \
  --output="$PLM_RUN/slurm_%j.out" --error="$PLM_RUN/slurm_%j.err" \
  "$PLM_ROOT/diagnostics/plm_pqq_delivery_20260923/run.sbatch" \
  "$PLM_RUN/REQUEST.json"
```

This is the measured resource profile: one H200, 32 allocated CPUs and
200000 MiB (about 195 GiB) host RAM. Use one scalar rank per task and the
wrapper's thread settings. The prior six-source scoring batch took 251 seconds
on this allocation, **excluding fresh protonation**. It is a measured batch
example, not a guarantee for a different pocket size or isolated-site latency.

For preparation alone, run the following instead of submitting. This uses the
existing one-thread CPU protonator, creates no endpoint energies, and leaves an
execution plan ready for the existing envelope wrapper:

```bash
"$PLM_PY" scripts/pqq_plm_prepare.py prepare \
  --request "$PLM_RUN/REQUEST.json" --output "$PLM_RUN/prepared"
"$PLM_PY" scripts/pqq_three_source_envelope_execution.py dry-run \
  --plan "$PLM_RUN/prepared/scoring/plan.json"
```

If preparation was already run, submit
`diagnostics/pqq_three_source_envelope_execution_20260923/run.sbatch` with
`$PLM_RUN/prepared/scoring/plan.json` under the same resources. Do not use the
combined wrapper again: output directories are immutable and repeat execution
is deliberately rejected.

## 4. Read results and handle failures

Each submitted batch writes:

- `prepared/READY.json`: every requested protein, including preparation failures.
- `prepared/scoring/RESULT_JOBID.json` and `.md`: per-source original/accommodated
  contrasts, decisions, common-pool choices and accommodation works; strict
  three-source medians and ranges; unavailable components and reasons.
- `prepared/scoring/execution_*.json`, endpoint receipts and logs: actual work,
  timing, source/implementation hashes and component results.

The actual filename includes the returned Slurm job ID. To re-collect a completed
or failed run, call the **pinned** script in `prepared/scoring/implementation/`
with `collect --plan ... --output ...`, using a new output filename, followed by
`report --result ... --output ...`. Collection does not rerun chemistry. Example
against the completed historical integration, producing a new report only:

```bash
"$PLM_PY" scripts/pqq_three_source_envelope_execution.py report \
  --result workspaces/pqq_three_source_envelope_execution_20260923/run_v1/RESULT_1211626.json \
  --output workspaces/plm_pqq_scan_example_v1/historical_report.md
```

**Use the operational median decision as the primary protein result.** Retain
all individual-source scores, their range, and the mathematical-minimum variant.
If one required member or matrix cell is unavailable, the protein result is
unavailable. An inconclusive member can coexist with a supported complete median;
report that disagreement. A boundary-reaching search is a bounded accommodation
proposal, not a claim of complete relaxation or an equilibrium state.

An SCF failure is not an inconclusive classification. Preserve the failed cell
and report unavailable. There is a separately validated two-start recovery for
one reference case, but it is **not an automatic general retry policy** in this
package. Any recovery result stays separately named, with its actual numerical
checks; never loosen stopping, swap a fold or select an electronic start by class.

## 5. Join batches to the biological evidence

Use the original protein table and a list of disjoint result files. The exporter
rejects duplicate protein results, mixed candidate references and sequence
mismatches. It preserves every original field, including DFT, motif, phylogeny,
genome provenance and missing transcript values. Preparation failures can be
supplied alongside the molecular results. Example on the real completed batch:

```bash
"$PLM_PY" scripts/pqq_plm_export.py \
  --proteins workspaces/nikasha_plm_export_20260922/export_v1/proteins.tsv \
  --results workspaces/pqq_three_source_envelope_execution_20260923/run_v1/RESULT_1211626.json \
  --output workspaces/plm_pqq_scan_example_v1/export
```

For new runs add their explicit JSON paths after `--results` and their
`prepared/READY.json` paths after `--preparations`. Avoid a recursive glob that
could include superseded retries or duplicate proteins. The output is
`proteins.tsv` plus `RECEIPT.json`. Unscored proteins remain `candidate_not_scored`;
preparation failures are `candidate_unavailable_preparation`; missing scores are
empty, never zero. Actual source-level details are in
`candidate_source_details_json` and the original result JSON.

The existing 176-row table is the old XoxF cohort. Before joining additional
PQQ families, extend the existing biological export using their actual full
sequence/gene/domain crosswalks. The exporter intentionally rejects absent IDs;
do not fabricate gene, genome or transcript associations. See
[the established joins and normalization definitions](../diagnostics/nikasha_plm_sources_20260922/SOURCES.md).

## 6. Interpretation and manuscript reporting

`R = E_Ca − E_La`; larger values are more La-like **on this protocol's scale**.
The operational frozen bands are:

| Raw composite R (model kcal/mol) | Interpretation |
|---|---|
| ≤ −405463.5325440529 | Ca-supported |
| ≥ −405456.4611478038 | La-supported |
| Between those values | Inconclusive |
| Missing required calculation | Unavailable |

The package reads these from the pinned canonical-only reference; the table is
documentation, not a second calibration. The large raw offset includes element
reference energies and is not a binding free energy. Do not compare it numerically
with old aquo-referenced DFT S, infer Kd/probability/physiological occupancy, or
assign catalytic substrate from the metal call.

Report the declared cohort denominator, preparation and numerical exclusions,
complete/inconclusive protein decisions, source spread and discordance. Preserve
the original DFT results as a separate column. The strongest tenfold benchmark
and practical threefold test are different context policies; use the correct
benchmark claim for this command. Ready-to-edit paper material and exact settings:
[manuscript draft](../diagnostics/plm_pqq_delivery_20260923/MANUSCRIPT.md) and
[technical supplement](../diagnostics/plm_pqq_delivery_20260923/SUPPLEMENT.md).
