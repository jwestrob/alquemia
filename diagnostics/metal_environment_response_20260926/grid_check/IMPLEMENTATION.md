# Eight-cell numerical follow-up — 27 September 2026

The original native derivative experiment passed selected finite differences,
step refinement and repeat checks, but failed its frozen rigid-transform gates.
This remains the original result. The new implementation prepares exactly eight
endpoints to distinguish translation sensitivity and numerical quadrature
sensitivity. No new endpoints ran during this implementation subtask.

Two original-grid translation-only A endpoints, plus refined A/B/combined-rigid-A
for each Ca and La. All elements, states, charges, source inventories and actual
embedding charges remain unchanged. Translation is (0.173,0.117,0.231) Å;
combined transformation additionally rotates +0.37 rad about z. A/B use exact
archived bytes. No force differences are evaluated in this new refined profile.

The only supported additional numerical profile is
`angular7_intacc7_unpruned_no_h_reduction_v1`:

```
%method
 AngularGrid 7
 IntAcc 7.0
 GridPruning Unpruned
 HGridReduced false
 DoEQ false
end
```

The existing native r2SCAN-3c/DefGrid3/TightSCF/EnGrad recipe otherwise stays
unchanged. `mace_omol_vacuum.parse_endpoint` accepts this explicitly named
optional grid profile, defaults unchanged, and requires exact scientific input
plus actual grid-generation header evidence. It rejects old-grid output under
the refined profile. Successful loading or echoed input alone is insufficient.
Actual refined-header parsing remains unexecuted until outputs exist; an
unexpected header is an explicit unavailable status, not fabricated success.

Collector retains original failed gate records, separate original-translation
checks, refined rigid residuals for energy and both core/MM gradients, each
metal's refined A/B response and its change from original, and the change in
metal-dependent double difference. Frozen rigid tolerances are 1e-5 Eh and
1e-4 Eh/bohr. Response changes use the original 0.05 kcal/mol uncertainty target.
Even passing results do not establish refined finite-difference qualification,
a complete hybrid Hamiltonian, or biological accuracy.

## Operations

`python scripts/metal_environment_grid_check.py prepare` requires explicit
`--source-manifest`, `--source-collection`, `--force-collection`, `--agreement`,
`--output`, `--workers` and `--mpi-ranks`. Root owns the new agreement and
submission. Use the actual completed scout and force collections, workers 8,
ranks 43 on the specified 344-CPU allocation; memory derives from actual Slurm
allocation through the existing runtime wrapper. Prepared implementation offers
`dry-run`, `execute`, and `collect`, each with `--manifest`; `collect --output`
writes an immutable result. No submission occurs in this script.

```bash
python -m unittest discover -s tests -p test_metal_environment_grid_check.py -v
```

Four real-artifact tests passed, zero skips, 10.33 seconds. These cover unchanged
old-grid parsing, rejection of an old real output as refined, actual geometry
invariance, temporary eight-cell preparation/dry-run, and missing-cell rejection.
No fabricated molecular outputs or chemical labels are used.
