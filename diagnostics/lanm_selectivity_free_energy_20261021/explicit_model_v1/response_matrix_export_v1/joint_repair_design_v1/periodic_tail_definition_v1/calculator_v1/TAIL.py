"""Certified scalar tail bounds for the frozen water-response lattice definition.

Pure Decimal arithmetic; no geometry/native/fit access. See SPEC.md for scope.
All input decimal quantities are exact strings. Output intervals enclose the
positive upper-bound expression, not the unknown physical tail itself.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import decimal as dc
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path


class BoundError(ValueError):
    pass


@dataclass(frozen=True)
class Interval:
    lo: dc.Decimal
    hi: dc.Decimal

    def json(self):
        return {"lower": str(self.lo), "upper": str(self.hi)}


def exact_decimal(value):
    if not isinstance(value, str) or len(value) > 500:
        raise BoundError("Decimal input must be an exact string of at most 500 characters")
    try:
        x = dc.Decimal(value)
    except dc.DecimalException as e:
        raise BoundError("Invalid decimal string") from e
    if not x.is_finite():
        raise BoundError("Nonfinite input")
    return x


class Arithmetic:
    """Outward enclosure using libmpdec's correctly rounded operations.

    exp/ln/sqrt always round HALF_EVEN: pad by one representable neighbor in
    each direction. All basic operations use FLOOR/CEILING directly. Underflow,
    subnormal results, exponent clamping and overflow are errors, never zero.
    """
    def __init__(self, precision=80, emin=-999999, emax=999999):
        if type(precision) is not int or not 40 <= precision <= 250:
            raise BoundError("Precision must be an integer in [40,250]")
        if type(emin) is not int or type(emax) is not int or not emin < 0 < emax:
            raise BoundError("Finite decimal exponent domain required")
        self.precision = precision
        self.down = dc.Context(prec=precision, rounding=dc.ROUND_FLOOR, Emin=emin, Emax=emax)
        self.up = dc.Context(prec=precision, rounding=dc.ROUND_CEILING, Emin=emin, Emax=emax)
        self.near = dc.Context(prec=precision, rounding=dc.ROUND_HALF_EVEN, Emin=emin, Emax=emax)
        for ctx in (self.down, self.up, self.near):
            for signal in (dc.InvalidOperation, dc.DivisionByZero, dc.Overflow,
                           dc.Underflow, dc.Subnormal, dc.Clamped):
                ctx.traps[signal] = True
        self.zero = self.point("0")
        self.one = self.point("1")

    def checked(self, lo, hi):
        if not lo.is_finite() or not hi.is_finite() or lo > hi:
            raise BoundError("Nonfinite or reversed enclosure")
        for x in (lo, hi):
            if x and not self.near.Emin <= x.adjusted() <= self.near.Emax:
                raise BoundError("Result outside normal decimal exponent domain")
        return Interval(lo, hi)

    def point(self, s):
        x = exact_decimal(s)
        # Exact input need not fit precision; directed operations round it later.
        return self.checked(x, x)

    def integer(self, n):
        if type(n) is not int:
            raise BoundError("Integer required")
        return self.point(str(n))

    def add(self, a, b):
        return self.checked(self.down.add(a.lo, b.lo), self.up.add(a.hi, b.hi))

    def sub(self, a, b):
        return self.checked(self.down.subtract(a.lo, b.hi), self.up.subtract(a.hi, b.lo))

    def mul(self, a, b):
        pairs = [(x, y) for x in (a.lo, a.hi) for y in (b.lo, b.hi)]
        return self.checked(min(self.down.multiply(x, y) for x, y in pairs),
                            max(self.up.multiply(x, y) for x, y in pairs))

    def div(self, a, b):
        if b.lo <= 0 <= b.hi:
            raise BoundError("Division interval includes zero")
        pairs = [(x, y) for x in (a.lo, a.hi) for y in (b.lo, b.hi)]
        return self.checked(min(self.down.divide(x, y) for x, y in pairs),
                            max(self.up.divide(x, y) for x, y in pairs))

    def power(self, a, n):
        if type(n) is not int or not 0 <= n <= 100000:
            raise BoundError("Only bounded nonnegative integer powers supported")
        out = self.one
        while n:
            if n & 1:
                out = self.mul(out, a)
            n //= 2
            if n:
                a = self.mul(a, a)
        return out

    def _transcendental(self, a, name):
        if name == "sqrt" and a.lo < 0:
            raise BoundError("Negative square-root interval")
        if name == "ln" and a.lo <= 0:
            raise BoundError("Nonpositive logarithm interval")
        fn = getattr(self.near, name)
        low, high = fn(a.lo), fn(a.hi)
        if name == "exp" and (low <= 0 or high <= 0):
            raise BoundError("Exponential underflow cannot become zero")
        # sqrt(0), ln(1), exp(0) are exactly known, avoiding subnormal padding.
        def padded(value, argument, direction):
            if (name == "sqrt" and argument == 0) or (name == "ln" and argument == 1):
                return dc.Decimal(0)
            if name == "exp" and argument == 0:
                return dc.Decimal(1)
            return self.near.next_minus(value) if direction < 0 else self.near.next_plus(value)
        return self.checked(padded(low, a.lo, -1), padded(high, a.hi, 1))

    def exp(self, a):
        return self._transcendental(a, "exp")

    def sqrt(self, a):
        return self._transcendental(a, "sqrt")

    def ln(self, a):
        return self._transcendental(a, "ln")

    def rational_power_positive(self, a, numerator, denominator):
        if a.lo <= 0 or type(numerator) is not int or type(denominator) is not int or denominator <= 0:
            raise BoundError("Positive base and integer rational exponent required")
        return self.exp(self.mul(self.div(self.integer(numerator), self.integer(denominator)), self.ln(a)))

    def total(self, xs):
        out = self.zero
        for x in xs:
            out = self.add(out, x)
        return out

    def nonnegative(self, value, name, positive=False):
        x = self.point(value)
        if x.lo < 0 or (positive and x.lo <= 0):
            raise BoundError(name + " must be " + ("positive" if positive else "nonnegative"))
        return x


def poly_mul(ar, a, b):
    out = [ar.zero] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] = ar.add(out[i+j], ar.mul(x, y))
    return out


def stirling2(n, k):
    if n == 0:
        return int(k == 0)
    if k == 0 or k > n:
        return 0
    row = [1]
    for _ in range(n):
        row = [0] + [j * (row[j] if j < len(row) else 0) + row[j-1]
                     for j in range(1, len(row)+1)]
    return row[k]


def geometric_moment(ar, rho, K, degree):
    """Enclose sum(k**degree*rho**k, k=K..infinity), all terms positive."""
    if not 0 < rho.lo <= rho.hi < 1 or type(K) is not int or not 0 <= K <= 100000:
        raise BoundError("Geometric series must have 0<rho<1 and bounded integer K")
    if type(degree) is not int or not 0 <= degree <= 7:
        raise BoundError("Only shell degrees 0..7 supported")
    denominator = ar.sub(ar.one, rho)
    outer = []
    for j in range(degree+1):
        inner = []
        for ell in range(j+1):
            factor = stirling2(j, ell) * math.factorial(ell)
            if factor:
                inner.append(ar.div(ar.mul(ar.integer(factor), ar.power(rho, ell)),
                                    ar.power(denominator, ell+1)))
        outer.append(ar.mul(ar.integer(math.comb(degree, j) * K**(degree-j)), ar.total(inner)))
    return ar.mul(ar.power(rho, K), ar.total(outer))


def normalization(ar):
    # C = [P(5)*exp(-5)]^-2 = 9 exp(10)/1849 exactly.
    return ar.div(ar.mul(ar.integer(9), ar.exp(ar.integer(10))), ar.integer(1849))


def radial_f(ar, r):
    if r.lo < 0:
        raise BoundError("Negative radial argument")
    t = ar.mul(ar.integer(2), r)
    p = ar.add(ar.add(ar.one, t), ar.div(ar.power(t, 2), ar.integer(3)))
    return ar.mul(ar.mul(normalization(ar), ar.power(p, 2)), ar.exp(ar.mul(ar.integer(-2), t)))


def shell_tail(ar, water_count, a, b, D, N, p):
    """Bound omitted weighted f for literal primary representatives and cube N."""
    if type(water_count) is not int or water_count < 1:
        raise BoundError("Positive finite primary water count required")
    if N not in (1, 2, 3) or p not in (0, 1):
        raise BoundError("Only declared cubes 1,2,3 and moments 0,1")
    if not (0 < a.lo <= a.hi <= b.lo and D.lo >= 0):
        raise BoundError("Invalid orthorhombic geometry bounds")
    # K > D/a. Floor of an outward upper quotient plus one is safe.
    quotient = ar.div(D, a)
    K = max(N+1, int(quotient.hi.to_integral_value(rounding=dc.ROUND_FLOOR))+1)
    if K > 100000 or K-(N+1) > 10000:
        raise BoundError("Unbounded/unsupported initial-shell work; no truncation substituted")
    if ar.mul(a, ar.integer(K)).lo <= D.hi:
        raise BoundError("Could not certify positive tail radius")
    root3b = ar.mul(ar.sqrt(ar.integer(3)), b)
    finite = []
    for k in range(N+1, K):
        rlower = max(dc.Decimal(0), ar.sub(ar.mul(a, ar.integer(k)), D).lo)
        rupper = ar.add(ar.mul(root3b, ar.integer(k)), D)
        term = ar.mul(ar.integer(24*k*k+2), ar.power(rupper, p))
        finite.append(ar.mul(term, radial_f(ar, ar.point(str(rlower)))))
    # Positive polynomial (2+24k²)*(sqrt(3)*b*k+D)^p*P(2*a*k)^2.
    la = ar.mul(ar.integer(2), a)
    P = [ar.one, la, ar.div(ar.power(la, 2), ar.integer(3))]
    polynomial = poly_mul(ar, [ar.integer(2), ar.zero, ar.integer(24)], poly_mul(ar, P, P))
    if p:
        polynomial = poly_mul(ar, polynomial, [D, root3b])
    rho = ar.exp(ar.mul(ar.integer(-4), a))
    moments = [geometric_moment(ar, rho, K, d) for d in range(len(polynomial))]
    tail = ar.total(ar.mul(c, moment) for c, moment in zip(polynomial, moments))
    tail = ar.mul(ar.mul(normalization(ar), ar.exp(ar.mul(ar.integer(4), D))), tail)
    expression = ar.mul(ar.integer(water_count), ar.add(ar.total(finite), tail))
    return {"N": N, "p": p, "units": "dimensionless" if p == 0 else "angstrom",
            "physical_tail_enclosure": {"lower": "0", "upper": str(expression.hi)},
            "positive_bound_expression_enclosure": expression.json(),
            "first_infinite_shell_K": K, "explicit_initial_shells": K-(N+1),
            "rho_enclosure": rho.json(), "polynomial_degree": len(polynomial)-1}


def bound_functions(ar, inputs):
    """Pure conditional h/H/explicit-gradient bound propagation, one metal.

    Caller certifies retained sums, omitted f|dr|, coefficient magnitudes and
    dipole norms. These are not obtained or qualified by this scalar function.
    h sums must aggregate primary-water norms for the final gradient formula.
    """
    required = {"T0", "T_displacement", "retained_abs_df", "retained_B0", "retained_B1",
                "abs_a", "abs_b", "abs_c", "abs_d", "kappa", "alpha", "exp_sigma_upper"}
    if not required <= inputs.keys():
        raise BoundError("Missing conditional bound input: " + str(sorted(required-inputs.keys())))
    x = {k: ar.nonnegative(inputs[k], k, positive=k in ("kappa", "alpha", "exp_sigma_upper")) for k in required}
    # exp(sigma)<=exp(2*(|c|+|d|)) for nonnegative environments.
    # Verify supplied constant rather than silently assuming the coefficient box.
    coefficient_bound = ar.exp(ar.mul(ar.integer(2), ar.add(x["abs_c"], x["abs_d"])))
    if x["exp_sigma_upper"].lo < coefficient_bound.hi:
        raise BoundError("Supplied exp_sigma_upper does not enclose the declared coefficient bound")
    V0 = ar.rational_power_positive(x["kappa"], 2, 3)
    V1 = ar.mul(ar.integer(5), x["kappa"])
    T0, displacement, DR = x["T0"], x["T_displacement"], x["retained_abs_df"]
    B0R, B1R, b = x["retained_B0"], x["retained_B1"], x["abs_b"]
    Tdf = ar.mul(ar.integer(4), displacement)
    B0T = ar.mul(V0, T0)
    B1T = ar.mul(ar.add(ar.mul(ar.integer(4), V0), V1), displacement)
    amplitude = ar.add(x["abs_a"], b)
    Eh = ar.mul(ar.integer(3), ar.add(ar.mul(amplitude, B0T), ar.mul(ar.mul(b, T0), B0R)))
    Edh = ar.mul(ar.integer(3), ar.total([
        ar.mul(amplitude, B1T), ar.mul(ar.mul(b, B0T), ar.add(DR, Tdf)),
        ar.mul(ar.mul(b, T0), B1R),
        ar.mul(ar.mul(b, B0R), ar.add(ar.mul(ar.mul(ar.integer(2), T0), DR), Tdf))]))
    L1 = ar.add(ar.mul(ar.integer(2), x["abs_c"]), ar.mul(ar.integer(6), x["abs_d"]))
    L2 = ar.add(ar.mul(ar.integer(4), x["abs_c"]), ar.mul(ar.integer(28), x["abs_d"]))
    Bbeta1 = ar.div(ar.mul(x["exp_sigma_upper"], L1), x["alpha"])
    Bbeta2 = ar.div(ar.mul(x["exp_sigma_upper"], ar.add(L2, ar.power(L1, 2))), x["alpha"])
    Ebeta = ar.mul(Bbeta1, T0)
    Edbeta = ar.add(ar.mul(ar.mul(Bbeta2, T0), DR), ar.mul(Bbeta1, Tdf))
    out = {"status": "conditional_bounds_only", "field_tail": Eh.json(),
           "field_directional_derivative_tail": Edh.json(), "H_operator_tail": Ebeta.json(),
           "H_directional_derivative_tail": Edbeta.json(), "V0_enclosure": V0.json(),
           "omitted_abs_df_upper": str(Tdf.hi),
           "full_self_consistent_force_error_certified": False}
    norms = {"abs_K", "metal_u_norm", "metal_v_norm", "max_water_u_plus_v_norm"}
    if norms & inputs.keys():
        if not norms <= inputs.keys():
            raise BoundError("Incomplete explicit-gradient dipole norm inputs")
        y = {k: ar.nonnegative(inputs[k], k) for k in norms}
        explicit = ar.mul(ar.div(y["abs_K"], ar.integer(2)), ar.add(
            ar.mul(ar.mul(y["metal_u_norm"], y["metal_v_norm"]), Edbeta),
            ar.mul(y["max_water_u_plus_v_norm"], Edh)))
        out["explicit_fixed_dipole_gradient_tail"] = explicit.json()
    return out


def calculate(config):
    if config.get("schema") != "orthorhombic_literal_tail_input_v1":
        raise BoundError("Unsupported input schema")
    ar = Arithmetic(config.get("precision", 80))
    sides = config["cell_sides_A"]
    if not isinstance(sides, list) or len(sides) != 3:
        raise BoundError("Exactly three orthorhombic side lengths required")
    lengths = [ar.nonnegative(s, "cell side", positive=True) for s in sides]
    a, b = ar.point(str(min(x.lo for x in lengths))), ar.point(str(max(x.hi for x in lengths)))
    reps = config["representatives"]
    if reps.get("convention") != "literal_primary_coordinates":
        raise BoundError("Only literal primary representatives are supported")
    D = ar.nonnegative(reps["displacement_upper_A"], "D")
    squared = exact_decimal(reps["squared_displacement_upper_A2"])
    if squared < 0 or Fraction(D.lo)**2 < Fraction(squared):
        raise BoundError("D does not bound the supplied exact squared displacement")
    water_count = config["water_count"]
    result = {"status": "certified_geometric_tail_bounds", "box_id": config.get("box_id"),
              "input": config, "arithmetic": {"precision": ar.precision,
                  "libmpdec_version": dc.__libmpdec_version__, "rounding": "outward Decimal intervals",
                  "underflow_policy": "error; never zero"},
              "geometric_bounds": {"a_A": str(a.lo), "b_A": str(b.hi), "D_A": str(D.hi)},
              "cubes": {str(N): {"T0": shell_tail(ar, water_count, a, b, D, N, 0),
                                  "T1": shell_tail(ar, water_count, a, b, D, N, 1)} for N in (1, 2, 3)},
              "scope": "geometry-only; no native key, periodic Hamiltonian, stability or accuracy admission",
              "representative_certificate_independently_replayed": False,
              "full_self_consistent_force_error_certified": False,
              "new_native_calls": 0, "fit_outputs_read": 0}
    if "conditional_propagation" in config:
        result["conditional_propagation"] = bound_functions(ar, config["conditional_propagation"])
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = Path(args.out)
    if out.exists():
        raise BoundError("Refusing to overwrite a result")
    raw = Path(args.config).read_bytes()
    config = json.loads(raw)
    result = calculate(config)
    result["configuration"] = {"path": str(Path(args.config).resolve()), "sha256": hashlib.sha256(raw).hexdigest()}
    result["implementation"] = {"path": str(Path(__file__).resolve()), "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


if __name__ == "__main__":
    main()
