# Mex saved-gradient replay diagnosis

Job12673 failed before any MACE or ALPB calls at exact dictionary equality of
`mapped_physical_gradient`. Preserve that job and all four accepted12570 GFN2
outputs. No gate or scientific model changes are authorized by this diagnostic.

Run one CPU Slurm job using the actual Python environment used by12673. Verify
the failed result, unchanged static configuration and all original local-field
source/runtime/raw-output pins. Invoke the unchanged saved-output parser and
force scatter for all four12570 rows. Compare energies, populations, raw QM and
point-charge gradients exactly. For each mapped component, total, force and
translation diagnostic, record differing value counts, maximum absolute error,
and the largest differing indexed pairs. Independently scatter the saved raw
gradients and require equality with the parser's fresh scatter.

Record actual interpreter path/prefix, package versions/locations and NumPy build
configuration. A shared Python binary hash does not establish identical installed
packages across environment paths. Numerical differences are findings to inspect,
not evidence to assume before the job. No molecular calls, new energies, optimizer
steps, source changes, relaxed tolerance or automatic static retry occur.
