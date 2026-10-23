# Anchored-strain refinement: source-only preparation

Root declaration54936c17 fixes three original diagonal strains and two half-steps,1.25e−5 and6.25e−6 Å. Exactly12 new inactive baseline gradient endpoints, no central rerun, build, initializer or factorization. Reuse the admitted cutoff12 baseline binary (not the term-observer binary), exact saved central gradient and original water/ion model, parameters, key, PME64³/order5/alpha0.4, cutoff12, VDW9/taper8.1, PCG and boundary settings.

The raw virial-only projection and the corrected derivative are both reported. For each original strain S about fixed metal position a, q has units Å and L0=1Å:

    dE/dq = [V:S − (S a)·sum_i grad_i E]/L0.

NUMERICS uses the original raw projection and literal-print arithmetic. Corrected analytic printing allowance adds sum_j |(S a)_j/L0| sum_i roundoff(print(G_ij)). Central differences and two-step tests use the unchanged original slope_checks code object and0.001kcal/mol/Å gate, with only the two prospectively declared steps supplied in private globals. Two-step allowance retains both complete endpoint+analytic print allowances. No imported helper/source/global is changed.

PREPARE reuses the frozen original parse_xyz/write_xyz methods, .17g serialization and exact stored tangent arrays. Each state follows the same arithmetic r′=r+qv, cell′_axis=cell_axis*(1+q), retaining the metal anchor. Actual saved/readback geometry roundtrip and all topology tokens are checked. Runtime validates pins and literal key/parameter identity; RUN.py is byte-identical to the admitted cutoff12 native runner. Each immutable state output is permitted once; no automatic retry or fallback. Build/root-attribution authority and independent source review plus per-state root release are required before native execution.

The collector requires all12 terminal state receipts and original center, reports9 raw and9 corrected checks, and retains incomplete/failing results. Scheduler cost sums only12 unique completed parent allocations; build and central are inherited, with no new cost. Historical failed9Å/12Å/term comparisons remain immutable. Passing this finite panel would establish only sampled baseline numerical strain consistency under the corrected chain rule; it would not qualify the active correction, physical accuracy, full periodic dynamics or a production model.

Resource envelope: original standard/standard_expand/banfield1CPU wrapper, default memory/time, no requeue. Every native child remains capped at3600s and3GiB with all numerical thread variables1. Root owns submission and subsequent interpretation. No source preparation action invokes the model.
