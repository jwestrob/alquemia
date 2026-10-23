# Anchored12 refinement source review

Passed; no correction requested. Reviewed PREPARE/BINDING/RUN/NUMERICS/COLLECT,
the frozen declaration/templates,33 saved controls and saved geometry preparation.
No controls, endpoint arithmetic, molecular programs, builds or factors were rerun.

The analytic projection correctly subtracts `(S a)·Σg/L0` from `V:S/L0` for the
fixed-metal deformation. This is applied to all three original diagonal strains.
Positive dE/dr tokens and the exact binary64 anchor/tensor constants are converted
to70-digit Decimal arithmetic. The corrected printing allowance is the original
virial allowance plus `Σ_j |(S a)_j/L0| Σ_i rounding(g_ij)`, a conservative triangle
bound. Raw and corrected projections/checks remain separate. No new arithmetic
or physical acceptance margin is introduced.

The inherited slope-check code is reused with only private STEPS changed to the
declared1.25e−5 and6.25e−6Å. The original0.001kcal/mol/Å gate, endpoint printing
allowance and two-step rule remain unchanged, including both analytic print
allowances in the inherited two-step bound. All9 corrected and9 raw checks are
reported; complete collection and numerical success are distinct. A complete
numerical failure returns a completed diagnostic with failed numerical status,
not an accuracy admission. Historical results are retained.

Exactly12 native slots cover three axes×two steps×two signs. The original stored
tangents, metal anchor, .17g serializer, full topology and orthorhombic/cutoff
checks are retained. Saved preparation records exact readback and all checks
passed. The central gradient is reused from the admitted cutoff12 baseline.
The runner is byte-identical to that baseline and binds its existing gradient_old
executable; keys/parameters remain literal center pins. No new central program,
build, initializer, active response certificate or factorization is permitted.

Each endpoint has one immutable output, one-program release and exact source,
build, input and receipt joins. Original1CPU/default-allocation scheduling and
3600s/3GiB child limits remain. Collector accounting includes exactly12 parent
allocations once, with central and historical build costs separate. Prospective
success would qualify only this sampled baseline derivative/convergence test.

Reviewer contribution disclosure: I authored the earlier source-only mesh/anchor
hypothesis note. I did not author these refinement implementation files or rerun
their preparation, controls or actual numerical test.
