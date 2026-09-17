# Saved-component replay

Run from the repository root; this reads actual completed outputs and submits
no calculation. Use a fresh output directory (existing results are immutable).

```bash
python scripts/mace_short_range.py --medium workspaces/mace_global_benchmark_20260916/mace_v1/medium/collection_job_1200701.json --large workspaces/mace_global_benchmark_20260916/mace_v1/large/collection_job_1200702.json --local workspaces/mace_local_correction_20260916/mace_v2/collection_job_1200717.json --plan diagnostics/mace_short_range_20260916/PLAN.md --output workspaces/mace_short_range_20260916/replay_v1
python -m unittest discover -s tests -p test_mace_short_range.py -v
```

One test passes on all actual saved component pairs, including preservation of
both failed alpha/GGR comparisons and unavailable classifications. This is a
parser/accounting regression; no new scientific executable call occurred.
