# Numerical parameter path: prepared, not submitted

The source-only preparation passed and rendered 18 exact windows for two fixed fixtures. Fourteen saved controls pass. Root owns the 20-call submission; no native program was executed here.

`PATH.py` resolves native endpoint pair parameters before interpolation. Water atom types 36/37 use VDW classes 34/35. Ordinary VDW diameters are halved before CUBIC-MEAN; explicit VDWPR is already the full pair Rmin and is written directly. HHG well depths, minimum-nonzero Thole and explicit pair overrides match the pinned native source. Polarizability changes continuously; metal charge stays +3 and mass stays 138.91. Every nonmetal input line and parameter remains unchanged, including hydrogen reduction factors.

`RUN.py` requires exact source/binary/runtime/CPU-control/GPU-parity pins and re-renders the frozen manifest before execution. It evaluates nine GPU windows per fixture, with spatial finite differences at the two midpoints, then the same two midpoint calculations on CPU. Actual native output must contain complete finite energy/gradient arrays and terminal norm summaries; zero exit alone is insufficient. Endpoint equivalence, CPU/GPU parity, spatial FD and two-step lambda FD gates remain those declared in `SPEC.md`. A failure stops the finite run and is retained.

## Root-ready artifacts

- `CONFIG_RUN.json`: frozen, runnable binding; `RUN.sbatch`: one H200, four allocated CPUs, one numerical thread, standard_expand, default memory/time, shared qualified CUDA12 runtime.
- Matching workspace `prepared_v1/MANIFEST.json`: `0c9587ee9907b63a2f334cd93f5f3d60ba39aa6f30eb35ce888aca139e2a9700`.
- `CONFIG_RUN.json`: `d724fec00eed7bb159b35064ffd46d44bc004d95c5793646a1a19636de5972f1`.
- `READY.json` binds preparation, source, tests and launch files. `TESTS_v1.json` records the actual 14-test saved run.

The initial restricted-parser draft and its free-text-header parsing limitation are preserved in matching workspace `draft_parser_v1/`; the released parser filters native keywords before parsing quoted fields. No scientific parameter was changed in that repair.

The mass cancellation statement applies to the same classical unconstrained degrees of freedom and conventions in both legs. Nuclear quantum/isotope effects and mismatched mass-dependent constraints are outside this test. The candidate parameters remain physically unqualified. This path test does not release dynamics, sampling, parameter fitting or a selectivity claim.
