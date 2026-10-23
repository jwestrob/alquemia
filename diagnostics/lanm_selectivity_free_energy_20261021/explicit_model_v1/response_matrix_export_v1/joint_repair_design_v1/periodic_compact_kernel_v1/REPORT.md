# Compact response weight: prototype and uniform bound

**The declared compact kernel passes its finite mathematical implementation checks.** A single outward scalar calculation bounds its difference from the original infinite environment-weight lattice sum throughout the declared domain. This is not a native-engine or total-force qualification.

## Definition and uniform domain

The frozen [declaration](DECLARATION.md) precedes calculation. [COMPACT.py](COMPACT.py) implements Ftilde=F S with the declared quintic switch from 20 to 21 Å, including S′F and all radial product derivatives. Its pure environment prototype uses Ftilde consistently in s, Q, image-excluding q, h, metal beta and their displacement Jacobians. It accepts artificial parameters; no fitted or selected coefficients were read. The same h is to enter both mutual right-hand sides, and beta enters H once per primary metal. This mathematical prototype does not replace the qualified native O(Q) implementation or perform polarization solves.

For orthorhombic cells with shortest side≥44 Å and ≤3980 primary waters, the support lies strictly below half the shortest side. Each primary water has at most one contributing image. Nearest-image evaluation is therefore exactly the complete lattice sum of the compact Hamiltonian. Nearest-image branch changes occur where the value and first two derivatives are zero. This is a different, explicitly approximate kernel from the original noncompact F; minimum image has not been substituted for its infinite or central+26 sum.

## Bound proof and result

With D(r)=3+6r+4r², F′/F=−8r(1+2r)/D. Hence F decreases and |F′|≤4F because 4D−8r(1+2r)=12+16r≥0. Also rF decreases for r≥20: the sign is set by −[16r³+4r²−6r−3], negative already for r≥1. For the switch, S′=−30t²(1−t)², so |S′|≤15/8. Its positive factorization S=(1−t)³(1+3t+6t²), with matched first/second endpoint derivatives, establishes nonnegativity, monotonicity and C² joins.

Nearest representatives give the componentwise nonnearest shell bound r≥44(k−1/2), with exactly 24k²+2 indices at maxnorm k. This does not use an invalid vector-radius D=22 assumption. Monotonic rF permits the sharper weighted shell bound l_k F(l_k), with l_k=44(k−1/2). Setting j=k−1 yields positive polynomials of degree6/7 times exp(−88) exp(−176j). The pinned, previously reviewed Decimal/Stirling arithmetic evaluates the infinite expressions with outward rounding at precision80; no lattice truncation or coordinate sampling is involved.

| Uniform scalar error | Upper bound, rounded upward for display | Units |
| --- | ---: | --- |
| Sum of absolute weight changes T0 | 2.572417e−24 | dimensionless |
| Sum of absolute radial-derivative changes T1 | 1.505293e−23 | Å⁻¹ |
| Sum of radius-weighted weight changes V0 | 5.151237e−23 | Å |
| Sum of radius-weighted derivative changes V1 | 3.013147e−22 | dimensionless |

V1 applies to a homogeneous strain of operator norm at most one. Exact outward endpoints, separate nearest and nonnearest terms, and units are in the workspace `actual_v1/RESULT.json` pinned by [READY.json](READY.json). Physical scalar errors lie between zero and those upper endpoints; the narrow positive-expression intervals are not lower bounds on the physical errors. No upper cell-size bound or minimum contact distance is needed for these scalar bounds.

## Finite controls and costs

[TEST_COMPACT.py](TEST_COMPACT.py) passed all141 scalar checks from a small fixed set of artificial inputs. They cover exact rational switch joins and slope identities, radial monotonic inequalities, first/second product derivatives, Cartesian/strain derivatives within the switch, artificial periodic/reimage identity and shell counts, invalid cell rejection, and all Cartesian h/beta chain derivatives for three artificial waters. The maximum artificial h-derivative difference was 2.174e−11; beta was 3.218e−12. The predeclared relative allowance was 2e−7, with an absolute floor1e−12 for ordinary artificial environment quantities and no absolute floor for the tiny transition-kernel checks. The control values and finite differences are not physical error criteria or a formal floating-program proof.

Controls: 0.01529045 s wall / 0.01450332 s process CPU. One scalar bound run: 0.10292194 s wall / 0.09522015 s process CPU. Both used one local thread and zero native calls, physical-box evaluations, response solves, optimizer calls or fit-outcome reads. The recorded costs cover the control/calculation bodies; source inspection and launcher overhead are excluded.

## Integration boundary

These results qualify the declared weight kernel and its uniform scalar approximation bound within this prototype. Retained-water collective factors change along with omitted weights; h/H error propagation must include those changes. Coefficient-specific h/H and explicit-gradient bounds, response stability, induced-dipole/native-derivative bounds, solver/PME arithmetic and total self-consistent energy/force/stress admission are still separate. The existing switch contact regularization, reaction forces, two mutual channels and strain convention must be retained when porting the formulas. No periodic native flag, sampling, fit change or physical-accuracy claim follows from this result. The historical 27-image plan remains unchanged.
