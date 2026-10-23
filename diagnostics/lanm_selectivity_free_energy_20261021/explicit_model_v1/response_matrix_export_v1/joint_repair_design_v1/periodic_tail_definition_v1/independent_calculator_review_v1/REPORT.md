# Independent periodic-tail calculator source review

**Passed for the declared source-only scope; no blocker found.** This review binds frozen `TAIL.py` (`87560816…`) and READY (`10e7b4e7…`). It does not evaluate either actual box, rerun the author's controls, or qualify a periodic engine. The reviewer did not author the calculator and previously reviewed its symbolic definition.

## Arithmetic and infinite sum

The basic interval operations use explicit floor/ceiling contexts and all signed product/division endpoints. Integer powers retain outward multiplication. Monotone exp, ln and sqrt use the nearest context and one representable neighbor on each side; the exact zero/unit exceptions are valid. This relies on the documented Decimal/libmpdec rounding behavior, rather than proving that library. The installed Python fallback sqrt source also explicitly selects HALF_EVEN. [Python Decimal documentation](https://docs.python.org/3.12/library/decimal.html)

Overflow, underflow, subnormal results, clamping, invalid operations and zero-containing divisors fail closed. The post-padding exponent check is necessary because neighbor operations need not raise a subnormal exception. The saved four edge controls exercise this exact path; the broad control with a similar name only rejects an invalid input earlier. Correctly rounded transcendental values compared at higher precision use the same library, so those controls test the enclosure implementation, not an independent transcendental oracle.

The shell count is exactly `(2k+1)^3−(2k−1)^3 = 24k²+2`. Literal representatives give distances between `max(0,ak−D)` and `sqrt(3)bk+D`. The chosen K is strictly above D/a; earlier shells use monotonic F, and the remaining positive polynomial uses the reviewed Stirling/geometric identity for degrees 0–7. All coefficients and summed terms are nonnegative. Unsupported work limits raise an error without truncating the infinite sum. The physical-tail interval correctly has lower endpoint zero; the lower endpoint of the computed upper-bound expression is not advertised as a physical lower bound.

## Metadata and conditional propagation

Both prepared inputs exactly copy the pinned inventory's orthogonal sides, 3,980-water count and literal metal–oxygen representative certificate. Exact rational arithmetic verifies the supplied D² inequality. This review verifies those metadata joins and authority hashes; it does **not** independently recompute the maximum distance from coordinates. The CLI also does not dereference the authority hashes. Root's release wrapper must verify the frozen source/READY/configuration/geometry/inventory and representative certificate before evaluation. No minimum-image or future-configuration bound is substituted.

The h/H and directional-derivative formulas match the independently reviewed definition, including the global coefficient bounds, contact-safe V bounds, displacement-dependent derivative mass, and fixed-dual explicit-gradient contraction. All retained sums, coefficient/dipole bounds and their units remain caller-certified inputs. A conditional call does not establish that these inputs came from a particular geometry or that its T0/T1 were the geometric outputs. The two prepared inputs request geometry only.

Full self-consistent force error remains uncertified: changing the image sum changes both induced channels and the native contractions. This source review supplies no periodic stability, total-force, physical-accuracy or sampling admission.

## Evidence and cost

The 58 saved mathematical controls and four saved exponent-edge controls all pass and bind the exact implementation. They were inspected, not rerun. The accompanying `SOURCE_CHECKS.json` records fresh pin/metadata consistency checks and their measured cost. All work used at most one local thread, with zero calculator evaluations, solver calls, native calls, fit-outcome reads or changes to frozen source. `REVIEW.json` contains the complete evidence pins and scope limitations.
