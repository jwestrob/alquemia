# Execution adjustments

The initial one-node128CPU request was rejected by Slurm before allocation;
compatible nonexcluded nodes provide64CPUs. The same eight scientific tasks were
split by structure into two64CPU manifests:1201824/1201825, four16-rank workers
each. No scientific input or run count changed. The original orientation_v1
manifest was never executed and is superseded by the two case manifests.

Initial reference job1201830 failed during MPI startup because the runner used
allocated-node CPU count rather than the requested single-rank policy. Explicit
execution_resources (mpi_ranks1/concurrent_tasks1) fixed the launcher. The default
runner policy and active protein jobs are unchanged. Fresh reference_v2 completed
as1201831, with old failed attempts retained. See WATER_REFERENCE.md for actual
cost and reference results. The cluster allocated a whole64CPU node even though
one rank was requested/used; the cost record includes all allocated CPUs.

Source preflight initially compared every chain/terminus in the two deposits.
1F6S contains additional chains and two extra modeled terminal residues. The
corrected integrity check uses their common site-chain residue mapping:120/120
shared residues identical. Every selected union fragment exists in both
structures. No source coordinates, sequence or scientific sample was changed.

Native Opt selects TightSCF (energy tolerance1e-8), confirmed in real output;
old square SPs used1e-6. State SPs with no mobile waters explicitly use TightSCF
to match optimizer electronic settings. Old/new contrast changes include the
representation and SCF-policy difference, which are reported separately from
within-network water reorientation.
