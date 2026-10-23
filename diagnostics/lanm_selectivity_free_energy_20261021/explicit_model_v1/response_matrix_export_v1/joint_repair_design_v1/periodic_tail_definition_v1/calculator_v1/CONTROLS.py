"""Small mathematical controls only; never reads actual box/fit/native artifacts."""
import decimal as dc
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("tail_calculator", HERE / "TAIL.py")
t = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = t
spec.loader.exec_module(t)


def exact_geometric_moment(rho, K, degree):
    # Independent rational differentiation of r*d/dr, then finite-prefix removal.
    polynomial = [Q(1)]
    for d in range(1, degree+1):
        deriv = [i*polynomial[i] for i in range(1, len(polynomial))]
        inner = [Q(0)] * max(len(deriv)+1, len(polynomial))
        for i, val in enumerate(deriv):
            inner[i] += val
            inner[i+1] -= val
        for i, val in enumerate(polynomial):
            inner[i] += d*val
        polynomial = [Q(0)] + inner
    whole = sum(coef*rho**i for i, coef in enumerate(polynomial)) / (1-rho)**(degree+1)
    return whole - sum(Q(k**degree)*rho**k for k in range(K))


def main():
    wall, cpu = time.monotonic(), time.process_time()
    checks = []
    def check(name, passed, detail=None):
        checks.append({"name": name, "passed": bool(passed), "detail": detail})
        if not passed:
            raise AssertionError(name)
    def contains(iv, exact):
        return Q(iv.lo) <= exact <= Q(iv.hi)
    def reject(name, fn):
        try:
            fn()
        except (t.BoundError, dc.DecimalException):
            check(name, True)
        else:
            check(name, False)
    ar = t.Arithmetic(80)
    check("signed_basic_arithmetic_encloses_exact_rational", contains(ar.div(ar.point("-2"), ar.point("3")), Q(-2, 3)))
    check("cancellation_preserves_exact_zero", ar.add(ar.point("-1"), ar.point("1")).lo == 0)
    check("integer_power_encloses_exact_rational", contains(ar.power(ar.div(ar.one, ar.integer(3)), 7), Q(1, 3**7)))
    # Exact published-looking digits are not used: high-precision libmpdec oracle
    # checks only padding/integration, not independent correctness of libmpdec.
    oracle = dc.Context(prec=200)
    for name, value in (("exp", "-0.123"), ("sqrt", "2"), ("ln", "0.125")):
        iv = getattr(ar, name)(ar.point(value))
        truth = getattr(oracle, name)(dc.Decimal(value))
        check(name+"_encloses_200digit_value", iv.lo < truth < iv.hi)
    check("zero_and_unit_transcendentals_are_exact",
          ar.sqrt(ar.zero).lo == ar.sqrt(ar.zero).hi == 0 and
          ar.ln(ar.one).lo == ar.ln(ar.one).hi == 0 and
          ar.exp(ar.zero).lo == ar.exp(ar.zero).hi == 1)
    two_thirds = ar.rational_power_positive(ar.point("8"), 2, 3)
    check("positive_fractional_power_encloses_exact_four", contains(two_thirds, Q(4)))
    for degree in range(8):
        for K in (0, 1, 3):
            iv = t.geometric_moment(ar, ar.point("0.5"), K, degree)
            exact = exact_geometric_moment(Q(1, 2), K, degree)
            check(f"rational_series_d{degree}_K{K}", contains(iv, exact))
    check("stirling_known_rows", [t.stirling2(5,k) for k in range(6)] == [0,1,15,25,10,1])
    f0 = t.radial_f(ar, ar.zero)
    check("radial_normalization_at_2_5", contains(t.radial_f(ar, ar.point("2.5")), Q(1)))
    mock = {"schema":"orthorhombic_literal_tail_input_v1", "box_id":"synthetic_math_only",
            "cell_sides_A":["3","4","5"], "water_count":2,
            "representatives":{"convention":"literal_primary_coordinates", "displacement_upper_A":"1.25",
                               "squared_displacement_upper_A2":"1.5"}}
    result = t.calculate(mock)
    # For this domain all N use K=N+1 and the same positive majorant.
    for p in (0,1):
        vals = [dc.Decimal(result["cubes"][str(N)]["T"+str(p)]["physical_tail_enclosure"]["upper"]) for N in (1,2,3)]
        check(f"mock_positive_monotone_T{p}", vals[0] > vals[1] > vals[2] > 0)
    # Independently summed finite part of the same majorant must fit inside
    # the infinite expression. This is an arithmetic control, not tail truth.
    for p in (0,1):
        a,b,D = ar.point("3"),ar.point("5"),ar.point("1.25")
        partial = []
        for k in range(2,14):
            rlo = ar.point(str(3*k-dc.Decimal("1.25")))
            rhi = ar.add(ar.mul(ar.sqrt(ar.integer(3)),ar.integer(5*k)), D)
            partial.append(ar.mul(ar.mul(ar.integer(2*(24*k*k+2)), ar.power(rhi,p)), t.radial_f(ar,rlo)))
        bound = dc.Decimal(result["cubes"]["1"]["T"+str(p)]["physical_tail_enclosure"]["upper"])
        check(f"direct_positive_shell_partial_below_T{p}", ar.total(partial).hi < bound)
    # D>a*(N+1) exercises the finite early shells with zero lower radius.
    early = t.shell_tail(ar,1,ar.point("1"),ar.point("2"),ar.point("5"),1,0)
    check("early_shell_branch", early["first_infinite_shell_K"] == 6 and early["explicit_initial_shells"] == 4)
    base = {"T0":"0.001","T_displacement":"0.002","retained_abs_df":"2",
            "retained_B0":"1","retained_B1":"3","abs_a":"0","abs_b":"0",
            "abs_c":"0","abs_d":"0","kappa":"0.125","alpha":"1","exp_sigma_upper":"1",
            "abs_K":"2","metal_u_norm":"1","metal_v_norm":"2","max_water_u_plus_v_norm":"3"}
    zero = t.bound_functions(ar,base)
    check("zero_response_exact_bypass", all(zero[k]["lower"] == zero[k]["upper"] and dc.Decimal(zero[k]["upper"]) == 0
          for k in ("field_tail","field_directional_derivative_tail","H_operator_tail",
                    "H_directional_derivative_tail","explicit_fixed_dipole_gradient_tail")))
    active = dict(base, abs_a="0.1",abs_b="0.2",abs_c="0.1",abs_d="0.1",exp_sigma_upper="4")
    response = t.bound_functions(ar,active)
    check("conditional_positive_bounds", all(dc.Decimal(response[k]["upper"]) > 0 for k in
          ("field_tail","field_directional_derivative_tail","H_operator_tail","H_directional_derivative_tail","explicit_fixed_dipole_gradient_tail")))
    # Exact algebra for V0=.25,V1=.625, active magnitudes above.
    T0,Tdisp,DR,B0R,B1R,A,b,V0,V1 = map(Q,(".001",".002","2","1","3",".3",".2",".25",".625"))
    Tdf,B0T,B1T = 4*Tdisp,V0*T0,(4*V0+V1)*Tdisp
    Eh = 3*(A*B0T+b*T0*B0R)
    Edh = 3*(A*B1T+b*B0T*(DR+Tdf)+b*T0*B1R+b*B0R*(2*T0*DR+Tdf))
    L1,L2 = Q(".8"),Q("3.2")
    Eb = 4*L1*T0
    Edb = 4*(L2+L1**2)*T0*DR+4*L1*Tdf
    for name, exact in (("field_tail",Eh),("field_directional_derivative_tail",Edh),("H_operator_tail",Eb),
                        ("H_directional_derivative_tail",Edb),("explicit_fixed_dipole_gradient_tail",2*Edb+3*Edh)):
        iv = response[name]
        check(name+"_rational_propagation", Q(iv["lower"]) <= exact <= Q(iv["upper"]))
    reject("reject_float_input",lambda:ar.point(0.1))
    reject("reject_nonfinite",lambda:ar.point("NaN"))
    reject("reject_zero_denominator",lambda:ar.div(ar.one,ar.zero))
    reject("reject_nonconvergent_series",lambda:t.geometric_moment(ar,ar.one,2,7))
    reject("reject_bad_displacement_certificate",lambda:t.calculate(dict(mock,representatives=dict(mock["representatives"],displacement_upper_A="1"))))
    reject("reject_zero_cell_side",lambda:t.calculate(dict(mock,cell_sides_A=["0","4","5"])))
    reject("reject_minimum_image_convention",lambda:t.calculate(dict(mock,representatives=dict(mock["representatives"],convention="minimum_image"))))
    reject("reject_false_exp_bound",lambda:t.bound_functions(ar,dict(active,exp_sigma_upper="1")))
    reject("reject_incomplete_gradient_norms",lambda:t.bound_functions(ar,{k:v for k,v in base.items() if k!="metal_u_norm"}))
    narrow=t.Arithmetic(40,emin=-10,emax=10)
    reject("reject_underflow_not_zero",lambda:narrow.exp(narrow.point("-100")))
    reject("reject_overflow",lambda:narrow.exp(narrow.point("100")))
    reject("reject_subnormal_padding",lambda:narrow.sqrt(narrow.point("1e-20")))
    source = HERE/"TAIL.py"
    output={"status":"passed_mathematical_mock_controls","implementation":{"path":str(source),"sha256":hashlib.sha256(source.read_bytes()).hexdigest()},
            "controls_implementation":{"path":str(Path(__file__).resolve()),"sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
            "checks":checks,"count":len(checks),"wall_seconds":time.monotonic()-wall,"cpu_seconds":time.process_time()-cpu,
            "thread_environment":{k:os.environ.get(k) for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS")},
            "actual_box_evaluations":0,"native_calls":0,"fit_outputs_read":0,
            "limits":"libmpdec correctness assumed; transcendental oracle is not an independent implementation"}
    dest = Path(sys.argv[1])
    dest.parent.mkdir(parents=True,exist_ok=True)
    with dest.open("x") as handle:json.dump(output,handle,indent=2,sort_keys=True);handle.write("\n")
    print(json.dumps({"status":output["status"],"checks":len(checks),"wall_seconds":output["wall_seconds"],"cpu_seconds":output["cpu_seconds"]}))

if __name__=="__main__":main()
