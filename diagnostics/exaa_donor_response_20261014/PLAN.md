# ExaA: test the saved donor response against independent vacuum DFT

Authorized14October2026: Jacob requested autonomous continuation and subsequent
experiments. Root owns execution and integration. Production, reference bands,
old outcomes and other sessions remain unchanged.

## Question and declared test

Does r2SCAN-3c support the large MACE-vacuum stabilization along the saved ExaA
lanthanum proposal? Added packing changed the ambiguous fold's differential work
only16.692→16.317kcal/mol. The conditionally prepared crystal has a much smaller
response, but that alone cannot distinguish a fold-specific response from a learned
surface error. This test isolates the electronic response at identical coordinates.

Use the already consumed two-Ca ExaA/Q9 sample3 target and sample0 control from the
previous packing experiment. No newly selected favorable structure. For each,
evaluate the exact saved192-atom origin and adaptive_La geometry under both Ca
and La: eight finite singlepoints. Both metals traverse the same physical path.
All atoms, protonation, dry inventory, disulfide, full donors and14caps remain
unchanged. The remote spectator remains outside the local region, exactly as in
the existing comparison. This does not test global assembly or redox assignment.

For each metal report W_M = E_M(q_La) − E_M(q0), then D = W_Ca − W_La.
Positive D is a La-directed change in the paired contrast. Compare native DFT
with the saved MACE vacuum works and report target-minus-control differences.
The operational scorer selects adaptive_Ca for calcium and adaptive_La for
lanthanum; this common-path diagnostic is deliberately distinct from its full
three-candidate row minima. The adaptive_Ca cells remain inventoried but unrun.
No DFT absolute score is assigned a MACE calibration band.

## Reference and admission

Freeze native ORCA6.1.1 `r2SCAN-3c NoAutostart DefGrid3 TightSCF`, vacuum,
Ca charge−3 / La−2, physical singlets. Native basis/auxiliary basis and La46-core
ECP must be verified in actual outputs; no added dispersion/composite correction.
Derive explicit electron counts from these coordinates. This is an independent
electronic treatment, not the MACE training Hamiltonian or an exact oracle.
No native-GFN2 parameter file belongs in these DFT inputs. Any later GFN2 work
must retain the corrected parameter guard and its qualified solver route.

Reuse the native task runner, allocation-aware renderer and corrected PedE TRAH
collector, all copied and pinned. Require native normal termination and convergence,
correct executable/input/output identity, basis/ECP/electron inventory, full
populations and retained GBW. For TRAH, require the actual final macro gradient
residual below the actual printed TolG/TRAH threshold. Keep other density residuals,
frontier orbital warnings and state limitations visible: convergence is not
electronic-state uniqueness. Preserve failures; no automatic seed or solver sweep.

No numerical agreement threshold is being fit to these outcomes. Report signed
individual and differential discrepancies continuously. Use1kcal/mol as a declared
resolution floor for qualitative interpretation, not a demonstrated DFT error bar
or scientific acceptance bound. A reversed or strongly reduced response warrants
checking reference state before changing the learned scorer. Agreement supports
only this conditional local finite motion. Further correction/promotion requires
the existing reference population, not a repair tuned to these two structures.

## Stages and cost

1. Run Slurm-only real-input/parser preflight on all eight planned cells.
2. Launch two target-La scout endpoints, origin and adaptive_La, as independent
   12-MPI jobs. Their pair supplies a meaningful individual response and actual
   convergence/runtime evidence. Arm completion and health wakes.
3. After inspecting both scouts, execute the remaining six cells if this reference
   is executable and interpretable. A failed endpoint stays missing; method changes
   require a named diagnostic and do not silently replace it. No new permission is
   needed within this authorized scope.
4. Join admitted results to saved MACE works on Slurm and interpret the matched
   response. CPCM, analytic-force evaluation, other geometries and optimization
   are separate possible follow-ups, not hidden additions to this matrix.

The prior PedE192-atom reference took about1hour for easier endpoints and over
12hours for its slow Ca-vacuum endpoint. This is targeted development validation,
not routine scoring cost. Shared standard queue, current12CPU/job60CPU/user limits,
no priority/memory/time override. Memory follows actual allocation with existing
operational headroom. No arbitrary project cap; record real costs and inspect
specific nonprogress rather than repeatedly spending on an unchanged failure.

## Reuse and delivery

`inventory/REPORT.md` and `method_review/REPORT.md` pin the actual source inventory
and bounded cache search. No exact compatible completed ExaA DFT cells or GBWs
were found. Older73-atom origin states and PedE GBWs are incompatible. Do not
reuse them just because a protein name or atom count resembles this input.

Retain each endpoint's unrounded energy, state/populations, warnings, execution
receipt, runtime and GBW. Deliver available/missing cells, common-path works and
method discrepancies, with no classifier or affinity claim from two structures.
