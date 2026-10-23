# Orthorhombic periodic-tail calculator

This source implements the approved infinite-image definition's geometric upper bounds and conditional propagation formulas. It changes no geometry, coefficient, native kernel or physical cutoff. The two actual fixture inputs are prepared from frozen metadata but **have not been evaluated**. Root must release that saved mathematical evaluation after source review. No fit output is read.

## Input and image convention

`TAIL.py calculate(config)` accepts exact decimal strings, three positive orthorhombic side lengths, a positive integer primary-water count and a certified upper bound D on the literal primary metal–oxygen distances. The representative convention must be `literal_primary_coordinates`. It checks D² against the supplied exact upper squared distance using rational arithmetic. This check establishes internal certificate consistency; it does not independently reread coordinates.

`PREPARE_INPUTS.py` binds the two rows of fixture inventory `ae407a66…`. It copies their exact side strings, 3980-water counts and literal D certificate, preserving geometry and inventory pins. The oxygen representatives slightly straddle a conventional primary cell. A common fundamental-cube assumption and a minimum-image half-diagonal are therefore not used. The future-configuration domain is unbound.

The geometry-only inputs do not admit their older water parameter keys or a periodic Hamiltonian. In particular, no `EWALD-BOUNDARY` keyword is rendered: the existing native parser selects VACUUM when that keyword is present regardless of its text. The tail function does not choose an electrostatic boundary condition.

## Mathematical output

For N=1,2,3, the exact definition retained is

```
T_p = Nw sum_(k>N) (24 k²+2) (sqrt(3) b k+D)^p
                      F(max(0,a k-D)), p=0,1;
F(r) = [9 exp(10)/1849] (1+2r+(2r)²/3)² exp(-4r).
```

For orthorhombic cells a=min(side), b=max(side). T0 bounds the omitted sum of F; T1 bounds the omitted sum of F times distance. The finite algorithm evaluates early shells directly and upper-bounds the remaining positive shell polynomial with the approved finite Stirling/geometric-series expression. K is chosen strictly above D/a and at least N+1. More than 10,000 early shells or K above 100,000 is an unsupported domain error, never a substituted truncation.

Each `positive_bound_expression_enclosure` encloses the explicitly computed **upper-bound expression**, which may itself exceed the exact T_p expression. Its lower endpoint is not a lower bound on the unknown physical tail. `physical_tail_enclosure` is always [0, upper]. The scalar bounds may be conservative; no physical or numerical acceptance tolerance is chosen.

## Arithmetic guarantee and failure policy

Inputs are exact finite decimal strings; binary floats and nonfinite values are rejected. Decimal arithmetic has 80 digits by default (allowed 40–250), with normal exponent domain −999999 to +999999. Addition, subtraction, multiplication and division use ROUND_FLOOR for lower endpoints and ROUND_CEILING for upper endpoints. Products/divisions examine all signed endpoint combinations. Integer powers use repeated interval multiplication. All shell-series terms are nonnegative.

CPython Decimal/libmpdec documents exp, ln and sqrt as correctly rounded to ROUND_HALF_EVEN, independent of the requested rounding mode. The implementation therefore evaluates each monotone endpoint and expands it to the preceding/following representable value. Exact sqrt(0), ln(1) and exp(0) are handled directly. Positive rational power kappa^(2/3) uses these outward ln/exp operations and outward multiplication/division. This guarantee relies on the documented correctness of the installed libmpdec implementation; the code is not a formal proof of that library.

Overflow, underflow, subnormal values, clamping, invalid operations and zero-interval division fail. Padded endpoints are also checked against the normal exponent range because next-minus/next-plus need not signal subnormal results. An exponential cannot return a zero tail. Unsupported/unbounded arithmetic raises an error; callers must preserve stderr and must not reinterpret failure as a successful zero or infinite bound. No cancellation-prone alternating series or approximate floating-point power is used.

## Conditional field, H and derivative API

`bound_functions(ar, inputs)` implements all four propagation expressions in DEFINITION.md section 4 and the explicit fixed-dipole gradient expression in section 5. Required exact-string **upper bounds** are:

- `T0`: omitted sum f.
- `T_displacement`: omitted sum f |dr| for the declared directional derivative.
- `retained_abs_df`: retained sum |df| across the complete metal environment.
- `retained_B0`, `retained_B1`: retained sums from section 4, aggregated across primary waters for a global field-gradient bound.
- `abs_a`, `abs_b`, `abs_c`, `abs_d`: coefficient magnitude bounds; `kappa`, `alpha`: positive fixed values.
- `exp_sigma_upper`: a certified exponential bound, checked against exp(2(|c|+|d|)). This prevents silently assuming the frozen coefficient box for arbitrary supplied amplitudes.

The implementation derives V0=kappa^(2/3), V1=5 kappa and omitted sum |df| ≤4 T_displacement. For a unit Cartesian derivative, T_displacement≤T0 is sufficient. For a unit-operator-norm homogeneous strain use the geometric T1 instead. A collective coordinate needs its declared displacement bound. Field/H units follow the native e/Å² and Å⁻³ convention; directional derivative units depend on the supplied coordinate/strain definition.

Optional `abs_K`, `metal_u_norm`, `metal_v_norm`, and `max_water_u_plus_v_norm` must all be supplied together to report an explicit fixed-dipole correction-gradient bound. They are conditional norm inputs, not computed or admitted here. The caller must certify these retained sums and bounds; the calculator only propagates them. For multiple metals, sum the per-metal triangle bounds under a separately admitted topology.

The two real input configurations currently request **geometry only**. No retained sums, coefficients, inverse-response bound or periodic dipoles are bound. Full self-consistent force error remains uncertified: changing images changes u/v and native frame/mutual contractions as well as the explicit h/H terms. No periodic stability, physical accuracy, sampling or production qualification follows.

## Controls and execution

The source controls use synthetic mathematical inputs only. They check degree-0…7 geometric moments against an independent exact rational generating-function recurrence; directed basic and fractional-power arithmetic; radial normalization; finite positive shell comparisons; early-shell handling; exact-zero response; independent rational propagation; and invalid/unbounded inputs. A focused exponent test separately verifies rejection when a valid normal exp result becomes subnormal only after outward padding. The original broad control's similarly named case rejects an out-of-domain input before sqrt; it is not counted as that focused post-exp proof.

The 200-digit transcendental comparisons use the same library and test enclosure plumbing, not an independent transcendental implementation. No actual box, native executable, geometry coordinate file or fit output is evaluated by the controls.

READY.json binds exact sources, controls, input configurations and future argv. The CLI refuses to overwrite an existing output. Each result records its exact input/source pins and arithmetic library version. Run serially with the declared one-thread environment only after root release. There is no automatic continuation to image experiments or a native implementation.
