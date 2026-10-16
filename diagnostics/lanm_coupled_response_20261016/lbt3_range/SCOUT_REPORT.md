# Range-only scout is admissible; torsion now limits the motion

7CCO/La continuation12544 completed in33s on12CPUs/oneH200. It reused the old
terminal point and original-source coordinates/metric; allfour actual saved
start/gradient/guard preflights passed. Two new composite evaluations supplied
two accepted steps:2MACE/4GFN calls, no failures or repeated origin calls.

Increasing only source displacement .20→.40Å gives another−0.266482kcal/mol
of own composite work (total−3.146811 from original). Final source extent
.225220Å, target .088740Å and torsion .198922rad pass all physical guards.
The four subsequent trials violate the retained .20rad torsion bound; no
stationary minimum is claimed. Active projected maximum force remains4.963665
kcal/mol/Å. No geometry/score threshold was changed after seeing this result.

Scheduler allocation:396CPU-s/33GPU-s; measured TotalCPU40.233s, MaxRSS1448008K.
Application32.295742s, warm-model loading6.562158s, MACE1.108977s and summed
concurrent GFN work31.580613s are distinct measures. Actual receipts and sources
are in SCOUT_RESULT.json/SCOUT_ACCOUNTING.tsv.

The execution/physical gate passes irrespective of small energy gain. Root
submitted the unchanged otherthree as12545,12CPUs/oneGPU, completionwatcher1947155;
configf3d4c206 under lbt3_range_remaining_v1. No complete paired range result yet.
