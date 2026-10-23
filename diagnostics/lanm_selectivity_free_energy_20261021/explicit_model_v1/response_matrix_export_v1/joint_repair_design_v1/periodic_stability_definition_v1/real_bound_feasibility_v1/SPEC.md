# Fixed real-space stability-bound feasibility

Question: does the proposed Ewald split permit a useful inexpensive lower bound
on the uncorrected BIO aqueous polarization operator at the two already saved
geometries? This is an engineering diagnostic, not a periodic energy calculation,
model-accuracy test, or admission of an unequilibrated box.

Use exactly the La and Dy geometry/topology pins in periodic_fixture_handoff_v1.
Use BIO water/salt alpha and Thole values, intrinsic metal alpha/Thole, and the
17676 metal–water POLPAIR values recorded in that handoff. Set the custom h/H
correction to zero for this baseline-only test. The old17605 keys remain intact
and are not executed. No candidate from the live fit is read.

Fix the engineering Ewald split at0.4 inverseAngstrom and real cutoff9Angstrom,
as in the earlier unexecuted periodic proposal. Do not scan either parameter.
These choices do not bind later simulation settings. Cell must be the49.324A
cube; real cutoff is below half its length. All mutual group scales are1 in the
pinned BIO file. Charge-penetration, extra images and vacuum surface terms are
outside this diagnostic. All atom types and bonds must match the declared roles.

Implement the exact real-pair algebra in pinned induce.f/umutual2a and
damping.f/dampewald+dampthole, including the native damping>=50 branch.
Use the default Thole minimum/zero fallback and the explicit metal–O overrides.
For the normalized matrix B=I-sqrt(alpha)(T_real+T_self)sqrt(alpha):

- Atom blocks: diagonal(1-alpha*tself)I, off-diagonal minus normalized pairT.
- Fixed molecular blocks: each intact water's9coordinates; singleton ion blocks.
  Include all within-block pairs in the diagonal block. No repartitioning/search.
- Compute Gershgorin block norm L/U for both partitions. Cross-block spectral
  norms preserve all Cartesian signs inside each block before taking the norm.
- Report the sites/blocks attaining the bounds, pair counts and timings.

If L<=0, this bound fails to certify even the real part. A negative bound does
not show that B or the complete PME operator is indefinite. If L>0, it is only
promising: reciprocal upper norm, application/rounding errors, nonzero H, native
parity and the unchanged relative spectral criterion still require qualification.
Ordinary double arithmetic is explicitly exploratory, not a rigorous certificate.

Implementation controls: Ewald0 vacuum algebra, normalized pair rotational
invariance, and a tiny direct full-matrix eigenvalue enclosure for both fixed
partitions. Controls use artificial matrices/geometries, no real box or native
calls. Freeze source/input pins and obtain an independent source review before
the one two-box evaluation. At most1 numerical thread; require a free slot
within the aggregate16thread login allowance. No native, fit, new reference,
coordinate change, engine change, held-out outcome or scheduler allocation.
