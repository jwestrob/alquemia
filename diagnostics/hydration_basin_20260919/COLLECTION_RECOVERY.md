# Collection recovery; no scientific rerun

DFT job1202052 completed all20 actual native EnGrad calculations normally.
Its original collection retained16 records and marked four unavailable with
`changed validation geometry`, from exact Python tuple equality when regenerating
rigid rotations on the compute node. The original collection remains untouched.

The collector now uses the same atom-order check and 1e-12 Å regeneration tolerance
already used by the validation manifest checker. Exact source XYZ hashes remain
mandatory. All20 source XYZ files also regenerate exactly on the login host;
no compute-node numerical difference magnitude was retained, so none is claimed.
Native output geometries, analytic gradients, input hashes and receipts are
independently parsed and verified. No coordinate, scientific threshold, Hamiltonian,
endpoint calculation or biological label changed. No additional DFT run occurred.

`validation_v1/dft/collection_recovered_v1.json` contains all20 actual results.
Frozen `analysis_implementation_v2` performs recovery/comparison; the later v3
snapshot adds the report's conditional thermal-extent diagnostic. Earlier snapshots
and the incomplete first collection remain available. This is a collection defect,
not a scientific nonconvergence or successful substituted calculation.
