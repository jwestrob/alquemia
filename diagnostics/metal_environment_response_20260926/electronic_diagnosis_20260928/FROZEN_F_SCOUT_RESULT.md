# First completed Dy frozen-f analytic endpoint

The real Hans50atom isolated origin converged in20SCF cycles with a matched native
analytic SCF/ECP/D4 gradient. Energy−1407.723125602833Eh; ECP55,208explicit electrons,
restricted valence multiplicity1, physicalDyIII4f9 sextet metadata retained.
This proves executable energy/gradient output, not force consistency or affinity.

Worker1220294 completed549s×40CPUs; collector1220295 took1CPU-second:
21961allocatedCPU-seconds, zeroGPU. ORCA molecular wall544.032s. This is a
reference-scout cost, not an inference throughput claim.

Initial frozen collector expected a named ECP; actual custom basis prints an
unnamed ECP55 header. The old collection error is retained. Parser correction
b6b3e64 accepts the actual header while still enforcing element/core count;
QUALIFIED_COLLECTION_v2.json records collector hash and successfully parses the
same immutable molecular output. No new chemistry or relaxed convergence.

Declared four physical metal displacements now run in worker1220300 with
collector1220301:10MPI×4tasks,40sharedCPUs,mem0. Both original and half-step central
differences must pass frozen0.1kcal/mol/A derivative error and0.05refinement
criteria. No claim of general force qualification before those results.
