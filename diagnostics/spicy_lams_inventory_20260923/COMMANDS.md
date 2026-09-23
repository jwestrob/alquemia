# Metadata-only commands

Run from the repository root with the already-installed system Python
(`/home/jwestrob/.pyenv/versions/3.12.0/bin/python`,openpyxl+gemmi). No installation
or Slurm allocation is required. Existing outputs are write-once; use a new output
path only for an explicitly needed repeat inventory.

The completed commands were:

```bash
python diagnostics/spicy_lams_inventory_20260923/build_inventory.py \
  --output workspaces/spicy_lams_inventory_20260923/inventory_v2
python diagnostics/spicy_lams_inventory_20260923/assemble_delivery.py \
  --inventory workspaces/spicy_lams_inventory_20260923/inventory_v2 \
  --output workspaces/spicy_lams_inventory_20260923/delivery_v1
```

Existing results: read`inventory_v2/ORTHOLOGS.tsv`,`MATURE_SEQUENCES.faa` and
`delivery_v1/STRUCTURE_COVERAGE.tsv`. Join by accession+sequence_sha256 or exact
ortholog_id; preserve reserved_panel and construct scope. This supplies no
molecular execution command, new prediction, per-ion outcome, or inferred label.

The methods-only source capture is
`workspaces/spicy_lams_inventory_20260923/PRIMARY_METHODS_ONLY.json`; public
assay-method text was extracted separately from outcome sections. Existing
reserved seals remain unchanged and verified in`SEAL_CHECK.txt`.
