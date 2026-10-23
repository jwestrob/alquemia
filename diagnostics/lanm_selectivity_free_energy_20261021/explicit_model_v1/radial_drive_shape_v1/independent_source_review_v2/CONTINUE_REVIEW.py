#!/usr/bin/env python3
"""Narrow source-only continuation; no author imports, controls or native work."""
import hashlib
import json
from pathlib import Path
import re
import time


def pin(p):
    p = Path(p).resolve()
    return {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


def read(p):
    return json.loads(Path(p).read_text())


HERE = Path(__file__).resolve().parent
BASE = HERE.parent
CHECKS = []
PINS = {}


def check(name, passed, **evidence):
    CHECKS.append({"name": name, "passed": bool(passed), **evidence})
    if not passed:
        raise ValueError(name)


def checked(p):
    actual = pin(p["path"])
    if actual != p:
        raise ValueError("pin changed: " + p["path"])
    PINS[p["path"]] = actual
    return Path(p["path"])


start, cpu = time.monotonic(), time.process_time()
out = {"status": "failed_radial_drive_shape_source_review", "checks": CHECKS,
       "new_native_calls": 0, "new_builds": 0, "new_solver_calls": 0,
       "author_control_reruns": 0, "reviewer_implementation": pin(__file__)}
try:
    rp = pin(BASE / "READY.json")
    check("unchanged_frozen_READY", rp["sha256"] ==
          "cec443d8a5c00fad44912efadf62e7ad114e4e6bdcad2dd33b4e267607f0803d")
    ready = read(rp["path"])
    previous = pin(BASE / "independent_source_review_v1/REVIEW.json")
    check("failed_v1_preserved", previous["sha256"] ==
          "ac6502f59a5d85f6f87a1a79977748d5a7907d8190bbe56c262202e60f4ea71c")
    old_review = read(previous["path"])
    check("v1_scope_exact", old_review["status"] == "failed_radial_drive_shape_source_review"
          and len(old_review["checks"]) == 42
          and all(x["passed"] for x in old_review["checks"][:41])
          and old_review["checks"][-1] == {"name": "same_complete_source_inventory", "passed": False})
    out.update(ready=rp, preserved_failed_review=previous,
               inherited_completed_checks=41, inherited_verified_pins=39)
    for key in ("implementation", "radial_helper", "source_manifest", "inputs", "domain",
                "specification", "qualification_plan", "local_comparison_plan",
                "controls_implementation", "controls"):
        checked(ready[key])
        out[key] = ready[key]
    inputs = read(ready["inputs"]["path"])["sources"]
    new = read(ready["source_manifest"]["path"])
    parent = read(checked(new["parent_source_manifest"]))
    check("unchanged_parent_authority", new["parent_source_manifest"] == inputs["parent_source_manifest"])
    out["parent_source_manifest"] = new["parent_source_manifest"]
    oldroot, newroot = Path(parent["source_root"]), Path(new["source_root"])
    def inventory(manifest, root):
        return {str(Path(p["path"]).relative_to(root)): p
                for p in manifest["source_outputs"] + [manifest["cmake"]]}
    old_files, new_files = inventory(parent, oldroot), inventory(new, newroot)
    omitted = {"cmake/0README", "cmake/TinkerOpenMM.cmake"}
    check("only_two_ancillary_manifest_omissions", set(old_files) - set(new_files) == omitted
          and not set(new_files) - set(old_files), parent_count=len(old_files), new_count=len(new_files))
    supplemental = []
    for rel in sorted(omitted):
        checked(old_files[rel])
        actual = pin(newroot / rel)
        check("actual_ancillary_copy:" + rel, actual["sha256"] == old_files[rel]["sha256"])
        supplemental.append({"relative_path": rel, "parent": old_files[rel], "actual_copy": actual,
                             "part_of_new_manifest": False, "binding": "supplemental review pin"})
    out["supplemental_ancillary_files"] = supplemental
    for p in old_files.values():
        checked(p)
    for p in new_files.values():
        checked(p)
    check("all_manifest_file_bytes_match_their_pins", True,
          unique_parent_files=len(old_files), unique_new_files=len(new_files))
    changed = sorted(rel for rel in new_files if old_files[rel]["sha256"] != new_files[rel]["sha256"])
    check("only_two_native_source_deltas", changed == ["source/induce.f", "source/polar.f"], changed=changed)
    for rel in changed:
        text = (newroot / rel).read_text()
        restored = re.sub(r"^c     RADIAL_SHAPE_BEGIN_([A-Z]+)\n.*?^c     RADIAL_SHAPE_END_\1\n",
                          "", text, flags=re.M | re.S)
        if rel == "source/induce.f":
            restored = restored.removesuffix("\n")
        check("literal_parent_restoration:" + rel, restored.encode() == (oldroot / rel).read_bytes())
    helper = Path(ready["radial_helper"]["path"]).read_text()
    induce = (newroot / "source/induce.f").read_text()
    check("appended_helper_exact", induce.count(helper) == 1)
    check("outer_and_inner_exact_zero_bypass",
          "if (lrshape(slot).ne.0.0d0) then" in induce and "if (e.eq.0.0d0) return" in helper)
    check("shape_setup_precedes_active_return",
          induce.index("call radialdrivesetup") < induce.index("if (.not.lractive) return"))
    check("finite_parser_domain_and_duplicates", all(s in helper for s in (
        "ieee_is_finite(e)", "abs(e).gt.1.0d0", "seen(slot)",
        "ia.ne.500.and.ia.ne.501", "ib.ne.36", "extra record tokens")))
    check("equal_opposite_coordinate_jacobians", all(s in helper for s in (
        "lrdh(:,:,oi,oi)=lrdh(:,:,oi,oi)+jac", "lrdh(:,:,oi,m)=lrdh(:,:,oi,m)-jac")))
    check("no_H_or_T_mutation_in_new_helper", not any(s in helper.lower() for s in
          ("lrbeta", "lrdb", "lrcoef", "uinp", "uind", "polarity")))
    check("unchanged_two_RHS_gradient_and_symmetric_virial_consumers",
          "(uind(k,i)+uinp(k,i))*lrdh(k,l,i,j)" in induce
          and "fac=0.5d0*electric/dielec" in induce
          and all(new_files["source/" + f]["sha256"] == old_files["source/" + f]["sha256"]
                  for f in ("empole1.f", "epolar1.f", "epolar3.f", "moments.f", "gradient.f")))
    cfg = read(BASE / "build_v1/CONFIG_TEMPLATE.json")
    builder = checked(cfg["parent_builder"])
    check("declared_build_explicitly_disables_OpenMM", '"-DTINKER_OPENMM=OFF"' in builder.read_text())
    cmake = Path(new["cmake"]["path"]).read_text()
    check("unchanged_conditional_OpenMM_include", "if(TINKER_OPENMM)" in cmake
          and "include(${TINKER_DIR}/cmake/TinkerOpenMM.cmake)" in cmake)
    out["build_scope"] = {"parent_builder": cfg["parent_builder"], "TINKER_OPENMM": "OFF",
                          "build_executed_by_review": False,
                          "scheduler_wrapper_review": "Separate root review; no native-source modification."}
    controls = read(ready["controls"]["path"])
    check("saved_36_controls_bound_without_rerun", controls["check_count"] == 36
          and len(controls["checks"]) == 36 and all(x["passed"] for x in controls["checks"])
          and controls["implementation"] == ready["controls_implementation"]
          and controls["source_manifest"] == ready["source_manifest"]
          and all(controls[k] == 0 for k in ("native_calls", "builds", "optimizer_calls", "molecular_or_reference_evaluations")))
    domain = read(ready["domain"]["path"])
    check("declared_domain_and_normalized_scales", domain["bounds"] == [[-1, 1], [-1, 1]]
          and domain["center"] == [0, 0] and domain["normalized_probe_step"] == .005
          and domain["actual_e_probe_values"] == [-.01, .01]
          and domain["normalized_trust_half_width"] == .025
          and domain["physical_e_trust_half_width"] == .05)
    plan = Path(ready["local_comparison_plan"]["path"]).read_text()
    check("common_signed_W9_safeguard_both_cases", "identical in both new cases" in plan
          and "|r_z+J_z delta z|≤|r_z(center)|" in plan
          and "Both minimize the SAME t" in plan and "Retain original25 physical gates" in plan)
    check("prospective_caps_unchanged", ready["caps"] == {
        "qualification_native_programs": 49, "qualification_internal_FD_energies": 132,
        "new_derivative_vectors": 4, "new_derivative_native_programs": 112, "affine_comparison_solves": 2})
    out["math_findings"] = [
        "With F=f phi psi, delta h=3e F r. Its Jacobian is 3e[F I+(F'/r)rr^T], with F'=(f'phi+f phi')psi+2r f phi/r0^2, as implemented.",
        "The shared h enters both direct RHSs. For U=-K/2 d^T A^-1 p, its explicit coordinate gradient is -K/2 (u+v)^T dh; unchanged native gradient consumes exactly that dual sum.",
        "Oxygen and metal coordinate derivatives are equal and opposite; original symmetric displacement-gradient virial is retained. New solved u/v feed unchanged permanent-frame forces.",
        "At contact, delta h=0 and J=-3e f(0)kappa I is finite. At ideal r0, delta h=0 while J=6e phi(r0)nn^T. Literal near-r0 coordinates are not clamped.",
        "The added polynomial-times-exponential drive decays at separation. Pair-additive drive does not make the self-consistent polarization energy pair-additive.",
        "The new helper changes neither H nor T. Fixed-geometry susceptibility is independent of field-independent h, so shape alone cannot repair chi.",
        "An ideal-r0 direct shift cannot repair central P/B or W9 moment. The meaningful local question is whether the two force-slope columns free useful motion of the original twelve columns under common safeguards.",
        "The same original25 physical gates, three electrical norm safeguards and signed W9 absolute-residual safeguard apply to both twelve- and fourteen-column comparisons; no new physical dipole tolerance is claimed.",
        "Source restoration establishes exact removed-delta equivalence. Actual new-binary e=0 admission/cache bridge remains required before reuse of the saved24 derivatives at exact feca."
    ]
    out["limitations"] = [
        "Source/math review only; no compiled/native qualification or new numerical experiment.",
        "The failed v1 metadata assumption remains failed and immutable. This continuation resolves its scope with two supplemental exact copied-file pins; neither is part of the new manifest.",
        "OpenMM ancillary module is outside the explicitly OFF build scope; this receipt does not qualify an OpenMM build.",
        "Finite added drive and declared e box do not establish global potential stability or physical accuracy.",
        "Local affine comparison cannot establish global infeasibility, and no coefficients/candidates are selected by this review.",
        "Native49-program qualification, binary/cache bridge, four changed-metal probes and two affine solves require their separate concrete execution authorities."
    ]
    out["status"] = "passed_radial_drive_shape_source_review"
    check("required_READY_contract_exact", all(out[k] == v for k, v in ready["source_review_contract"].items()))
except Exception as exc:
    out["status"] = "failed_radial_drive_shape_source_review"
    out["error"] = repr(exc)
finally:
    out["check_count"] = len(CHECKS)
    out["verified_pin_count"] = len(PINS)
    out["cost"] = {"metadata_review_wall_seconds": time.monotonic() - start,
                   "metadata_review_CPU_seconds": time.process_time() - cpu, "numerical_threads": 1}
    with (HERE / "REVIEW.json").open("x") as stream:
        stream.write(json.dumps(out, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"review": pin(HERE / "REVIEW.json"), "status": out["status"],
                      "check_count": out["check_count"], "cost": out["cost"], "error": out.get("error")}))
