# Terminal-candidate La BIO real operator preparation

Prepare one immutable assembly input for the selected component-fit vector
e41fcd3d, the composed La BIO box, and its admitted initialized PME snapshot.
No assembly, symbolic analysis, factorization, native call, or molecular energy
is executed during preparation. A separate root review and release must bind the
exact source, configuration, output directory and resource limits before assembly.

The future single-thread worker calls the unchanged qualified ASSEMBLE.py once
for normalized, unshifted real-space B0+H, with intrinsic positive alpha, literal
Thole parameters and candidate POLPAIR overrides. All four MUTUAL scales are
explicitly one in the BIO parameter file; direct and permanent-field scales are
different quantities. The frozen native mutual-field source supplies that mapping.
Site order, coordinates, type/alpha map, orthorhombic cell, 0.4 inverse-A split,
9 A cutoff and PME snapshot are joined literally. No effective-alpha substitution.

H follows periodicresponsekernel/periodicresponsegeometry in the frozen periodic
prototype: nearest image, original normalized water weight up to 20 A, positive
C2 switch to zero at 21 A. In water atom order, q accumulates 2*f_i*prefix_i;
s=sum(f_i), eta=q/(1+s*s), g=2*s/(1+s), sigma=g*(c+d*eta).
The normalized metal diagonal is exp(sigma)-1, including the collective d term.
The native diagonal before normalization is this value divided by intrinsic
metal alpha. The drive coefficients a,b are retained as provenance; they enter
the RHS, not this Hessian. No stiffness value is evaluated during preparation.

Future output is lower CSC float64/int32 .npy data/indices/indptr, unshifted site
upper bound U_B, detailed H inputs/readback and exact existing graph hashes.
Expected dimension 35874, 1815015 pairs, 16371009 lower entries. Changed damping
and H change values, not geometry graph. Original symbolic data values are not
reused. No reciprocal matrix or threshold shift is added here.

Future worker limits: one numerical thread, 8 GiB address space, 4 GiB assembler
estimate ceiling, 180 s CPU and supervisor wall limit; no retries. Pair-query
allocation precedes the assembler's count-based estimate and is covered by the
process cap. Root chooses launch timing. The later 12-thread factor is separate.

The manifest binds nominal tau contract 2aa3ac43 and the La initialized U_R and
delta_odd exactly. e_app remains null, with no established application-error bound.
This is ordinary-double operator preparation: neither full native operator parity,
periodic stability, physical accuracy, force correctness nor MD readiness follows.
The original composed model remains prepared_only_not_executable.

Preparation performs only pin/input/topology/site-map checks. Six finite synthetic
checks cover weight normalization, compact cutoff/half-switch, collective cross
sum, zero-stiffness bypass, and rejection of non-unity mutual scaling. They do not
call ASSEMBLE or any backend, repeat prior qualification, or evaluate a real-box H.
