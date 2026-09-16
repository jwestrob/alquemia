# PLM ADH9 fixed-core preparation — 16 September 2026


## Completed result — 16 September 2026

Slurm1199976 completed successfully in10m55s: all four endpoints validated.
Both compatible sites are **indeterminate** under the frozen fixed-core bands:

- PLM2_30 `32301_3` / PQQSEQ_48f861015fad150af40a: S=15.503189972320033 kcal/mol.
- PLM0_60 `4380_6` / PQQSEQ_13d74836d4b7a3e02140: S=21.45008326919374 kcal/mol.

Ca-supported≤14.857129202922806; La-supported≥23.460061205609236.
The other four selected models remain unscored. No inference of metal dependence
or substrate specificity follows from these two indeterminate scores.

Detailed result/validation records:
`/groups/banfield/users/jwestrob/EastRiver/EastRiver_PLM/revision_analysis/2026-09-11_PQQ_ADH/energetics_queue/adh9/execution/`
(`results.tsv`, `results.json`, `completion.json`, `terminal_review.json`).
The watcher delivered completion and stopped. The four-protein AF3 proposal
remains pending; no new calculations were launched by this completion review.


Preparation completed: **two La/Ca pairs ready; four selected models unsupported**.
No ORCA calculation, new fold, control rerun, alternate-model selection, or job
submission occurred in this preparation step.

The user supplied `adh9_itol_ids.txt` and approved resuming ADH9 energetics.
The parent session recorded that scope in PLM
`revision_analysis/2026-09-11_PQQ_ADH/energetics_queue/adh9/authorization.json`.
The six existing selected models and sequence-supported residue mappings were
reviewed in `existing_site_review.json` in that same directory. Their source
hashes, reviewed roles, original selection/qualification records and summary
confidence files are frozen in `candidate_manifest.json`. The manifest does
not substitute another model when the selected model fails.

## Current outputs

All candidate products are under:

`/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/plm_adh9_fixed_core_20260916/`

`prepared_pairs.json` is the execution handoff (`plm.adh9.prepared_pairs.v1`).
It retains all six cases, their source hashes, status, failure reason or the
hash-linked carve manifest. Each ready case directory contains its normalized
and protonated coordinates, provenance, complete atom ledger, and two native
ORCA input/XYZ pairs. No electronic energies are present.

| Target / selected rank | Outcome | Fixed protein core | Detail |
|---|---|---|---|
| PQQSEQ_48f861015fad150af40a / 0 | Ready | E205 N270 D321 D323 K345 | CN8; protein–La confidence 0.990677; raw Asp–Lys 3.487501 Å |
| PQQSEQ_13d74836d4b7a3e02140 / 0 | Ready | E201 N269 D320 D322 K344 | CN8; protein–La confidence 0.991584; raw Asp–Lys 3.351242 Å |
| PQQSEQ_2de92e200673b3ce1bb7 / 1 | Unsupported | Retained in manifest | Asp–Lys 3.602184 Å exceeds 3.5 Å |
| PQQSEQ_1a5e88d41a1e9303c550 / 1 | Unsupported | Retained in manifest | Asp–Lys 3.536213 Å exceeds 3.5 Å |
| PQQSEQ_cd0303ddb2d784459166 / 2 | Unsupported | Retained in manifest | Asp–Lys 3.541260 Å exceeds 3.5 Å |
| PQQSEQ_91654f14790b0a6e438d / 1 | Unsupported | Retained in manifest | CN6 fails the production CN7 gate |

Both prepared clusters contain **78 atoms**, including all **27 PQQ atoms**.
Their derived La/Ca total charges are **−2/−3**, with singlet multiplicity.
Zero heavy atoms were repaired, removed, or added. Heavy-coordinate changes
are limited to the frozen PDB serialization tolerance of 0.001 Å. La/Ca
nonmetal coordinates are byte-identical.

## Protocol and implementation

Scientific protocol:
`pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3`.

`prepare_candidates.py` supplies candidate-specific source and role authority
while directly importing the **unchanged, hash-verified** calibration helpers:
`fixed_core_carver.py`, `protonate_standard_only.py`, and their pinned canonical
dependencies. It does not modify the 25-protein calibration map or the inbox.
The original `implementation_pins.json` and `result.json` under
`diagnostics/pqq_pmdh_fixed_core_calibration_20260914/` remain calibration
authority. All 13 pinned helper hashes matched at preparation and validation.

Preparation is pH 7 standard-residue-only protonation, with the original
deterministic seed, package versions and one-thread OpenMM CPU configuration.
The complete PQQ(3−) microstate and side-chain/cap atom construction come from
the original helpers. Raw source models already use CCD PQQ and LA identifiers;
normalization here only serializes those coordinates to PDB. Unsupported
source identities, alternate conformers, extra chains, existing hydrogens,
heavy-atom repairs, or source/hash discrepancies fail explicitly.

All five core roles must be supplied explicitly. The anchor Glu/Asn,
catalytic Asp and its Arg/Lys partner are retained; the D+2 homolog is included
when acidic, exactly as in fixed-core v3. Core membership is independent of a
metal-distance radius. The catalytic-Asp/cation O···N gate remains 3.5 Å.
Source and prepared donor ledgers must agree, and no qualifying donor may
come from outside the fixed core. Production admission remains typed CN≥7
at 3.1 Å, ≤2 direct N donors, and protein–La confidence≥0.9; the calibration-only
CN6 exception is not used.

Electronic inputs use the unchanged ORCA 6.1.1 native
`r2SCAN-3c NoAutostart CPCM(Water) DefGrid3` recipe. No water, point charges,
geometry relaxation, custom La basis, apo arm or aquo rerun is added.

The released supported bands are read from the original calibration result:
Ca-supported `S ≤ 14.857129202922806`, La-supported
`S ≥ 23.460061205609236` kcal/mol, otherwise indeterminate. The aquo term is
a common reporting gauge; positive S alone is not a biological preference
classification. Motif/charge confounding and model-conditioned geometry still
limit interpretation. No score has been assigned during preparation.

## Interface and checks

The input is an explicit `plm.adh9.candidate_manifest.v1` JSON with `targets`.
Each target specifies `case_id`, `target_id`, `rank`, hash-linked `source_cif`
and `summary_confidence`, the protein chain, exact metal/PQQ selectors, all five
role selectors, and original source provenance. Output is a new workspace;
the wrapper refuses to overwrite an existing run.

Historical command for this completed preparation (do not repeat into the
existing workspace):

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
 /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
 diagnostics/plm_adh9_fixed_core_20260916/prepare_candidates.py \
 --manifest diagnostics/plm_adh9_fixed_core_20260916/candidate_manifest.json \
 --implementation-pins diagnostics/plm_adh9_fixed_core_20260916/implementation_pins.json \
 --output workspaces/plm_adh9_fixed_core_20260916
```

`preparation_validation.json` records PASS from `validate_preparation.py`:
all six cases retained; source/helper/input/coordinate hashes checked;
zero-repair and coordinate invariants checked; all 77 nonmetal atoms matched
to their fragment provenance in each pair; complete PQQ and charge ledgers
confirmed; bad source hashes rejected. The separate PLM execution agent owns
job packaging, runtime resource normalization, execution and score collection.
