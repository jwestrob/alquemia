# Continue the declared matched Ca/La reference

The two target-La scouts reproduce the large saved MACE response closely:
native vacuum DFT work is −26.70294112086986kcal/mol, versus MACE vacuum
−26.782803599660028kcal/mol. The difference is +0.07986247879016872kcal/mol.
These are the admitted completed endpoints joined by Slurm4111, not in-flight
SCF energies or cross-Hamiltonian absolute totals.

Both endpoints pass the17 declared native checks, including final TRAH residuals
8.804145e−7 and5.466705e−7 below their printed1e−5 threshold. Charges, electron
count, basis/ECP, coordinates, executable identity and physical state are fixed.
Near-final small/negative-gap warnings and level shifting remain explicit:
neither solution is independently stability-qualified. This agreement supports
the particular conditional local La response; it does not prove a ground-state
surface, real protein conformational population or metal preference.

The origin took6026s and the proposal6494s, each with12 allocated CPUs:
150240 allocatedCPU-s combined, zeroGPU. The scalar join took less than one
reported wall second on oneCPU (scheduler TotalCPU0.077s). Prior preflight adds
72allocatedCPU-s. These are one-time reference costs, not scanner throughput.

## Next execution, within the original eight-cell plan

The reference is executable and informative enough to finish the declared
conditional comparison. Submit the remaining six original cells without changing
method, state, coordinates, tolerances or calibration:

1. sample3 Ca origin and adaptive_La: establish the metal-dependent work.
2. sample0 La origin and adaptive_La: matched control's La response.
3. sample0 Ca origin and adaptive_La: control differential and source contrast.

Each is a separate12-MPI shared standard job with allocation-derived memory,
no priority/memory/time override. The scheduler enforces the current user limit;
no whole-node or bypass request. No scouts repeat. Use task-local locks and the
same immutable manifest d37733923c3203f41b73d5a3c1789d396e22f8da69e624df2c2d724d0141bcf1.
Arm individual completion/health watchers and one afterany matrix collector.
Missing or failed cells remain visible. No automatic restart, new optimization,
CPCM calculation or third geometry is added.

This weakens a MACE-only explanation of the large La work. Calcium is still
needed: a correct La response does not establish a correct differential response.
The source/control comparison can then distinguish a generally large response
from the extra response of the ambiguous fold. Only the complete matched result
should guide the next method-improvement decision. No classifier improvement or
production promotion is claimed by this scout.
