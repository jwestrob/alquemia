# Experimental accommodation results joined to the existing PLM evidence

The existing176-protein export now has a separate candidate overlay for the two
actually computed PLM triples. All34 original columns and176 rows retain their
original values and order. Nineteen appended fields retain the candidate's own
protocol/reference, three source scores, median, spread, decisions, selected
Ca/La geometries and accommodation works. The174 other proteins explicitly say
`candidate_not_scored`; missing energies are empty, never zero.

| Protein | Same-context origin range | Accommodated range | Three-source prediction |
|---|---:|---:|---|
|PQQSEQ_07ab500e3df76b30d71c|40.138064|2.307280|Ca-supported|
|PQQSEQ_83440678cbbd658047c9|72.416250|17.230067|Ca-supported|

Ranges are modelkcal/mol. These are unknown biological preferences, not
classification tests. The8344 sample1 remains inconclusive individually. The
envelope candidate remains experimental: its25canonical/3crystal calls pass,
but the A0A3Ca3 benchmark probe is inconclusive and full transfer is separate.
The exact [molecular report](../pqq_three_source_envelope_execution_20260923/REPORT.md)
retains all sources, components, search limits and costs.

## Existing biological evidence is retained

Both proteins are unbinned in the current export and have the existing ENDD motif
assignment. The07ab RNA profile is missing;8344 has a profile. Those facts and
the original source-gene/domain/tree identities are copied unchanged. Existing
sample/genome normalization remains in the original relational export; this
join creates no new transcript values, genome associations or biological labels.

Every candidate source's pinned AF3 input supplies the full protein sequence.
Its exact SHA256 and length agree with the existing protein export. The join
also checks the actual plan/reference/request pins, all declared source rows and
three-source group membership. It does not match only abbreviated identifiers.

OriginalDFT `S_kcal_mol` and candidate raw `R_model_kcal_mol` use different gauges;
do not subtract them. The original DFT predictions remain separately accessible.

## Artifacts and checks

Output: `workspaces/plm_candidate_overlay_20260923/export_v1/proteins.tsv`,
SHA256`3a6eadbd22485c97e6f5c37bc9ff4bd2db9b7558f8422dc82fc21206a28a5cd1`.
Adjacent`RECEIPT.json` records input/result/reference/code hashes and identities.
The script refuses to overwrite an existing output directory.

Four real-artifact tests passed in0.466s, zero skips: all original values and
missingness, actual source/score copies, rejection of a deliberately corrupted
copy of the full-sequence hash and a real result with one deleted source row.
No fabricated scientific result, new molecular calculation or production change.

Authorization: Jacob's September23 discretionary overnight work, plus the
existing PLM delivery scope; this is an export of completed authorized work.

## Runnable reproduction

From the repository root, the command below creates a new export directory;
it launches no molecular calculations. The executed directory was`export_v1`.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/plm_candidate_overlay.py \
  --proteins workspaces/nikasha_plm_export_20260922/export_v1/proteins.tsv \
  --result workspaces/pqq_three_source_envelope_execution_20260923/run_v1/RESULT_1211626.json \
  --output workspaces/plm_candidate_overlay_20260923/export_replay_v1
```
