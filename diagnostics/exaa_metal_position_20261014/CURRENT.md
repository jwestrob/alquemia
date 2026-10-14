# ExaA metal placement: implementation in progress

The complete eight-cell DFT comparison supports the large target and smaller
control response. Root selected the finite reciprocal metal-position test in
PLAN.md to isolate one constraint in the donor-only search. No new preparation
or scoring has executed yet. No new production, calibration or physical state.

Curie owns mapped preparation; Heisenberg owns reuse of the pinned seven-thread
warm MACE executor; Gauss owns independent review. Root alone submits jobs and
integrates results. First execute preparation/physical checks on Slurm, then
evaluate the eight declared vacuum cells only if their unchanged gates admit
them. Existing original energies are reused by exact scientific identity.

Do not interpret implementation readiness as molecular execution. All outputs
belong in the matching workspace. See the preceding completed reference at
`../exaa_donor_response_20261014/FINAL_REPORT.md` and its actual result pins.
