# Saved post-relaxation contacts

The runner is ready and actual real-artifact scout **4859 passed** all four checks
(`0:0`, one second/one CPU, zero energy calls); do not submit it again. Its receipt is
`workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_contact_scout_v1/SUBMISSION.json`.

Scout CONFIG SHA256:
`f7b1ddd29a75ae3ae2d278cf7b29309d1bf34805c22529627060c5e762ecd8b5`.
The same analyzer checks actual4831's completed Ca0 four geometries and nine
real unexecuted source rows. It requires paired coordinates, exact named atom
identities, unchanged4626 baseline donor distances and explicit missing proposals.
This tests data handling; it adds no molecular evaluation.

The source passed AST and shell syntax checks. Independent static review found
no remaining schema or missing-candidate blocker. Actual Slurm verification is
recorded in `SCOUT_REPORT.md`; no numerical analysis was run on the login host.
The nine real unexecuted rows remain explicitly missing their relaxed proposals.
The completed Ca0 control does not resolve sample4's named H–H contact.

## Final ten-source analysis

After the actual remaining-stage collection and terminal EXECUTION receipt exist,
freeze them with this metadata-only command from the repository root:

```bash
/usr/bin/python3 diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/contact_review/ANALYZE.py freeze \
  --collection workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_remaining_v1/execution/FINAL_COLLECTION.json \
  --output workspaces/pqq_failure_causal_synthesis_20261015/asn_relaxation_contact_review_v1
```

The command prints the exact new CONFIG path and SHA256. Root then submits the
copied `RUN.sbatch` with those two arguments, on one standard CPU with no memory
or time override, and arms the normal completion wake. The script refuses to
overwrite an existing output directory. Do not reuse the scout CONFIG for this
final collection or bind an unfinished collection.

Outputs are `execution/RESULT.json`, `TABLE.csv` (all ten sources × four
geometries) and `EXECUTION.json`. Track the **HD22–Gly292 backbone-H pair itself**,
with HD21 and the Asn backbone-H pairs retained separately. The output also keeps
OD1/ND2–metal distances, direct-bond-excluded physical H-contact minima, original
optimizer/final-source gates and boundary flags. Missing candidates remain
unavailable. Contact lengthening is descriptive; it supplies no new energy,
biological label, or unconstrained-minimum claim. See `METHOD.md` for the frozen
definitions.
