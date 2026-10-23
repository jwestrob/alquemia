# Saved electronic redistribution diagnostic

Question: does the electronic redistribution produced by a complete water shell
differ from the sum of isolated metal–water interactions, and does that difference
change between La/Dy and eight/nine waters?

Use only the admitted common-parent-basis 68 new plus 12 reused raw ORCA outputs.
Extract the printed Mulliken and Loewdin atomic charges and electronic, nuclear,
and total dipole vectors. Retain every ghost-center population. Use the same
linear combinations as the energy decomposition, separately within each parent:

- C = complex − ion − shell.
- P = sum(metal–water) − N ion − sum(water).
- B = C − P.

Report complete per-center vectors and the metal-center term under both population
schemes. Report dipole C/P/B in atomic units. Check the expanded B expression,
zero total charge of each difference and nuclear-vector cancellation. This avoids
assuming that ghost populations vanish or cancel between unlike parent geometries.

Before execution, require one complete indexed population table per scheme and
one dipole block per output, matching center labels and the full parent size.
Printed atomic charges must sum to the known molecular charge within 5e-5 e;
derived difference sums within 5e-4 e. Require dipole total = electronic + nuclear
within 5e-9 atomic units, nuclear cancellation within 1e-6 atomic units and
algebraic closure within 1e-10. These are printed-output reconstruction checks,
not scientific accuracy tolerances. Preserve failures instead of substituting data.

Population partitions are basis-dependent diagnostics, not observable atomic
charges or a unique charge-transfer energy. Compare both schemes and retain
disagreement. Dipole differences are charge-neutral and origin-independent when
the coordinate/nuclear cancellation checks pass. Their energy decomposition does
not uniquely separate polarization, Pauli effects, or charge transfer. No new
model parameters or amplitudes will be inferred from these values alone.

One local thread, zero molecular calls, no geometry changes or new outcome reads.
Keep the original outputs and common-basis comparison immutable. Record source,
raw-output and collector pins and actual elapsed/CPU cost. The goal authorizes
this contained saved-output diagnostic; it does not release a production model.
