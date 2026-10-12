#!/usr/bin/env python3
"""Fail-closed parameter contract for new ORCA native-GFN2 research tasks.

This module performs no molecular calculations and changes no production runner.
Parameter consistency, runtime qualification, SCF acceptance and biological
validation are separate claims. Runtime approvals can cover new geometries;
each task's geometry/state/source pins still participate in its cache identity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import tomllib

DATA = Path(__file__).resolve().parents[1] / "params/native_gfn2_parameters_v1"
ORIGINAL = DATA / "original_native.xtb.json"
AUTHOR = DATA / "author_gfn2.toml"
ORIGINAL_SHA = "f3dac846a59238947d35c23ba5bdb40c19a34482fe9860e8e11bf9c00d719107"
AUTHOR_SHA = "76239337c5f7f3285b0b4db06d607490d14cb38a9db90334f212e370f328082a"
FIELDS = ("kcn", "shpoly", "lgam", "refocc")
# Predeclared export/decimal roundoff tolerance, never an energy tolerance.
ABS_TOL = REL_TOL = 1e-12
LOCAL_PARAMETERS = "native_gfn2_parameters.json"
LOCAL_MANIFEST = "native_gfn2_guard.json"
METHOD = "Native-GFN2-xTB"
STATE_CONVENTION = "original-gfn2-spin-independent-charge-and-multiplicity-bookkeeping"


class GuardError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise GuardError(message)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def _pairs(pairs):
    out = {}
    for key, value in pairs:
        require(key not in out, f"duplicate JSON key: {key}")
        out[key] = value
    return out


def load_json(path):
    def invalid(value):
        raise GuardError(f"nonfinite JSON constant: {value}")
    return json.loads(Path(path).read_text(), object_pairs_hook=_pairs,
                      parse_constant=invalid)


def pin(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": sha256(path)}


def check_pin(record):
    require(isinstance(record, dict) and isinstance(record.get("path"), str),
            "file pin needs path and sha256")
    require(re.fullmatch(r"[0-9a-f]{64}", str(record.get("sha256", ""))),
            "file pin has invalid SHA-256")
    path = Path(record["path"])
    require(path.is_absolute() and path.is_file(), f"missing absolute pinned file: {path}")
    require(sha256(path) == record["sha256"], f"file digest changed: {path}")
    return path


def _diff(old, new, path=()):
    """Every semantic leaf difference, including insertions and type changes."""
    if type(old) is not type(new):
        return [{"path": list(path), "before": old, "after": new}]
    if isinstance(old, dict):
        require(old.keys() == new.keys(), f"JSON keys changed at {path}")
        return [d for key in old for d in _diff(old[key], new[key], (*path, key))]
    if isinstance(old, list):
        require(len(old) == len(new), f"JSON array length changed at {path}")
        return [d for i, (a, b) in enumerate(zip(old, new)) for d in _diff(a, b, (*path, i))]
    return [] if old == new else [{"path": list(path), "before": old, "after": new}]


def audit_parameters(parameters, policy):
    """Audit a complete JSON against the pinned export and author shell labels.

    policy names individual element/field corrections and pins candidate bytes.
    Unchanged fields compare exactly, NOT with the author roundoff tolerance.
    No global array permutation or implicit element correction is performed.
    """
    policy = load_json(policy) if not isinstance(policy, dict) else policy
    require(policy.get("schema_version") == 1 and policy.get("correction_id"), "invalid correction policy")
    require(policy.get("original_sha256") == ORIGINAL_SHA == sha256(ORIGINAL), "original export pin mismatch")
    require(policy.get("author_sha256") == AUTHOR_SHA == sha256(AUTHOR), "author parameter pin mismatch")
    candidate_sha = sha256(parameters)
    require(policy.get("candidate_sha256") == candidate_sha, "candidate full-file pin mismatch")
    old, candidate = load_json(ORIGINAL), load_json(parameters)
    author = tomllib.loads(AUTHOR.read_text())["element"]
    corrections = policy.get("corrections")
    require(isinstance(corrections, dict), "explicit correction mapping required")
    for element, fields in corrections.items():
        require(element in old["element"] and element in author, f"unknown correction element: {element}")
        require(isinstance(fields, list) and fields and len(fields) == len(set(fields))
                and set(fields) <= set(FIELDS), f"unsupported correction fields for {element}")
    differences = _diff(old, candidate)
    for change in differences:
        path = change["path"]
        require(len(path) == 4 and path[0] == "element" and
                path[2] in corrections.get(path[1], []),
                f"unapproved parameter change: {path}")
    elements = {}
    for symbol, native in candidate["element"].items():
        require(symbol in author, f"element absent from pinned author data: {symbol}")
        trusted = author[symbol]
        shells, reference_shells = native["shells"], trusted["shells"]
        require(len(shells) == len(set(shells)) and set(shells) == set(reference_shells),
                f"shell identities differ for {symbol}")
        fields = {}
        for field in FIELDS:
            values, reference = native[field], trusted[field]
            require(len(values) == len(shells) and len(reference) == len(reference_shells),
                    f"wrong shell coefficient count: {symbol}.{field}")
            # Full shell identity, e.g. 5d, not an index or a blanket s,p,d rule.
            expected = [reference[reference_shells.index(shell)] for shell in shells]
            require(all(type(v) in (int, float) and math.isfinite(v) for v in values),
                    f"nonfinite/nonnumeric parameter: {symbol}.{field}")
            passed = all(math.isclose(x, y, abs_tol=ABS_TOL, rel_tol=REL_TOL)
                         for x, y in zip(values, expected))
            fields[field] = {"passed": passed, "native": values, "author_by_shell": expected,
                             "max_abs_residual": max(abs(x-y) for x, y in zip(values, expected))}
            if field in corrections.get(symbol, []):
                require(passed, f"declared correction does not match author: {symbol}.{field}")
        elements[symbol] = {"shells": shells, "fields": fields,
                            "parameter_consistent": all(v["passed"] for v in fields.values()),
                            "runtime_qualified": False}
    return {"schema_version": 1, "parameter_sha256": candidate_sha,
            "original_sha256": ORIGINAL_SHA, "author_sha256": AUTHOR_SHA,
            "correction_id": policy["correction_id"], "corrections": corrections,
            "float_tolerance": {"absolute": ABS_TOL, "relative": REL_TOL},
            "differences": differences, "elements": elements,
            "scope": "shell-parameter audit only; no runtime or biological qualification"}


def _parameter_name(binding):
    name = binding.get("parameter_filename", LOCAL_PARAMETERS)
    require(isinstance(name, str) and re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*\.json", name),
            "parameter file must be a safe task-local JSON basename")
    return name


def import_directives(filename=LOCAL_PARAMETERS):
    filename = _parameter_name({"parameter_filename": filename})
    return (f'%method\n XTBParamFile "{filename}"\n'
            ' ReadXTBParam true\n WriteXTBParam false\nend\n')


def xyz_elements(path):
    lines = Path(path).read_text().splitlines()
    require(len(lines) >= 2, "XYZ missing count/comment")
    try:
        count = int(lines[0])
    except ValueError as exc:
        raise GuardError("invalid XYZ atom count") from exc
    require(count > 0 and len(lines) == count + 2, "XYZ count/trailing data mismatch")
    elements = set()
    for line in lines[2:]:
        row = line.split()
        require(len(row) == 4 and re.fullmatch(r"[A-Z][a-z]?", row[0]), "unsupported XYZ atom row")
        require(all(math.isfinite(float(x)) for x in row[1:]), "nonfinite XYZ coordinate")
        elements.add(row[0])
    return sorted(elements)


def _clean_input(path):
    # Fail-closed subset of ORCA input syntax; unsupported includes are not parsed.
    text = Path(path).read_text()
    require("$" not in text, "multi-job/include ORCA input unsupported")
    return "\n".join(line.split("#", 1)[0] for line in text.splitlines())


def validate_input(binding, variant):
    path = check_pin(binding["input"])
    text = _clean_input(path)
    simple = " ".join(line.strip()[1:] for line in text.splitlines() if line.lstrip().startswith("!"))
    methods = re.findall(r"\b(?:Native-)?(?:sp)?GFN[012]-xTB\b|\bXTB[012]?\b", simple, re.I)
    require(len(methods) == 1 and methods[0].lower() == METHOD.lower(), "input must select only Native-GFN2-xTB")
    require(binding["model"].get("method") == METHOD, "model method differs")
    solvent = re.findall(r"ALPB\s*\(\s*([^)]*)\)", simple, re.I)
    medium = "vacuum" if not solvent else f"ALPB({solvent[0]})"
    require(len(solvent) <= 1 and not re.search(r"\b(?:CPCM|SMD|GBSA)\b", text, re.I), "unsupported solvent setup")
    require(binding["model"].get("medium", "").lower() == medium.lower(), "medium binding differs")
    for keyword, wanted in (("ReadXTBParam", "false" if variant == "native-default" else "true"),
                            ("WriteXTBParam", "true" if variant == "native-default" else "false")):
        values = re.findall(rf"^\s*{keyword}\s+(\S+)\s*$", text, re.I | re.M)
        require([x.lower() for x in values] == [wanted], f"require exactly one {keyword} {wanted}")
        require(len(re.findall(rf"\b{keyword}\b", text, re.I)) == 1, f"ambiguous {keyword}")
    files = re.findall(r'^\s*XTBParamFile\s+"([^"\n]+)"\s*$', text, re.I | re.M)
    filename = _parameter_name(binding)
    count = len(re.findall(r"\bXTBParamFile\b", text, re.I))
    if variant == "native-default":
        require(files == [] and count == 0 and filename == path.stem + ".xtb.json",
                "native-default diagnostic must use the explicit input-stem export, without import filename")
    else:
        require(files == [filename] and count == 1, f"require exact task-local XTBParamFile {filename}")
    xyz = re.findall(r'^\s*\*\s+xyzfile\s+(-?\d+)\s+(\d+)\s+([^\n]+?)\s*$', text, re.I | re.M)
    require(len(xyz) == 1 and not re.search(r"^\s*\*\s+xyz\s", text, re.I | re.M), "require one xyzfile geometry")
    charge, multiplicity, name = xyz[0]
    require((path.parent / name.strip('"')).resolve() == check_pin(binding["xyz"]).resolve(), "XYZ input path differs")
    require(int(charge) == binding["state"]["charge"] and int(multiplicity) == binding["state"]["multiplicity"],
            "charge/multiplicity input differs")
    return path


def runtime_identity(binding, parameter_sha):
    return {"parameter_sha256": parameter_sha,
            "executable_sha256": binding["executable"]["sha256"],
            "runtime_files": sorted([{"role": p["role"], "sha256": p["sha256"]}
                                     for p in binding["runtime_files"]], key=lambda p: p["role"]),
            "environment": binding["environment"], "method": METHOD,
            "state_convention": binding["state_convention"]}


def build_binding(*, input_path, xyz, parameters, runtime_binding, source_inputs,
                  protocol_id, charge, multiplicity, medium,
                  parameter_filename=LOCAL_PARAMETERS, observed_runtime=None):
    """Construct (do not execute) a task binding from the recorded runtime."""
    runtime = load_json(runtime_binding)
    require(runtime.get("schema") == "native_runtime_binding.v1", "unsupported runtime binding schema")
    executable = runtime["driver"]
    runtime_files = {str(Path(p["path"]).resolve()): p for p in runtime["files"]}
    if observed_runtime is not None:
        observed = load_json(observed_runtime)
        require(observed["source_binding"]["sha256"] == sha256(runtime_binding), "observed runtime source pin differs")
        for p in observed["observed_files"]:
            key = str(Path(p["path"]).resolve())
            require(key not in runtime_files or runtime_files[key]["sha256"] == p["sha256"], "observed runtime hash conflicts")
            runtime_files[key] = p
    return {"schema_version": 1, "protocol_id": protocol_id,
            "parameter_filename": parameter_filename, "parameter_sha256": sha256(parameters),
            "input": pin(input_path), "xyz": pin(xyz), "elements": xyz_elements(xyz),
            "source_inputs": [pin(p) for p in source_inputs], "executable": executable,
            "runtime_files": [dict(p, role=Path(p["path"]).name) for p in runtime_files.values()
                              if Path(p["path"]).resolve() != Path(executable["path"]).resolve()],
            "environment": dict(runtime["environment"], **runtime["thread_environment"]),
            "state_convention": STATE_CONVENTION,
            "state": {"charge": charge, "multiplicity": multiplicity},
            "model": {"method": METHOD, "medium": medium}}


def validate_binding(binding, variant, *, check_environment=True):
    require(binding.get("schema_version") == 1 and binding.get("protocol_id"), "binding needs schema_version/protocol_id")
    require(binding.get("state_convention") == STATE_CONVENTION, "explicit original-GFN2 state convention required")
    require(type(binding["state"]["charge"]) is int and type(binding["state"]["multiplicity"]) is int
            and binding["state"]["multiplicity"] > 0, "invalid charge/multiplicity")
    check_pin(binding["executable"])
    require(isinstance(binding.get("runtime_files"), list), "runtime_files must be explicit, including empty diagnostic lists")
    roles = [p.get("role") for p in binding["runtime_files"]]
    require(all(roles) and len(roles) == len(set(roles)), "runtime file roles missing/duplicated")
    require(isinstance(binding.get("source_inputs"), list) and binding["source_inputs"], "source input pins required")
    for p in binding["runtime_files"] + binding["source_inputs"]:
        check_pin(p)
    require(isinstance(binding.get("environment"), dict) and binding["environment"], "explicit environment pins required")
    for key, value in binding["environment"].items():
        require(isinstance(key, str) and (value is None or isinstance(value, str)), "invalid environment pin")
        if check_environment:
            require(os.environ.get(key) == value, f"environment changed: {key}")
    path = validate_input(binding, variant)
    elements = xyz_elements(check_pin(binding["xyz"]))
    require(binding.get("elements") == elements, "declared elements differ from actual XYZ")
    return path, elements


def validate_approval(approval_record, binding, audit):
    """Check reviewer-issued runtime registry entry and its real control receipts.

    The registry binds runtime/helper/library/environment identity. It is not an
    independent proof that every loaded dependency has been enumerated. Reviewers
    must obtain these pins from the actual qualification runtime.
    """
    approval_path = check_pin(approval_record)
    approval = load_json(approval_path)
    require(approval.get("schema_version") == 1 and approval.get("status") == "qualified_for_research",
            "runtime approval unavailable")
    require(binding["runtime_files"], "research requires actual helper/library pins")
    require(approval.get("runtime_identity") == runtime_identity(binding, audit["parameter_sha256"]),
            "runtime approval identity differs")
    runtime_pin = approval.get("runtime_binding")
    require(runtime_pin is not None, "current runtime binding receipt required")
    runtime = load_json(check_pin(runtime_pin))
    require(runtime.get("schema") == "native_runtime_binding.v1", "unsupported runtime binding schema")
    require(runtime.get("driver", {}).get("sha256") == binding["executable"]["sha256"], "runtime driver pin differs")
    expected_files = {(str(Path(p["path"]).resolve()), p["sha256"]) for p in runtime["files"]}
    require(approval.get("observed_runtime") is not None, "observed runtime helper/library pins required")
    observed = load_json(check_pin(approval["observed_runtime"]))
    require(observed["source_binding"]["sha256"] == runtime_pin["sha256"], "observed runtime source differs")
    expected_files.update((str(Path(p["path"]).resolve()), p["sha256"]) for p in observed["observed_files"])
    actual_files = {(str(Path(p["path"]).resolve()), p["sha256"])
                    for p in [binding["executable"], *binding["runtime_files"]]}
    require(actual_files == expected_files, "runtime helper/library file set differs")
    expected_environment = dict(runtime["environment"], **runtime["thread_environment"])
    require(binding["environment"] == expected_environment, "runtime environment scope differs")
    require(binding["model"]["medium"] in approval.get("qualified_media", []), "medium lacks runtime qualification")
    require(set(binding["elements"]) <= set(approval.get("qualified_elements", [])), "element lacks runtime qualification")
    require(binding["state"]["multiplicity"] in approval.get("qualified_multiplicities", []),
            "multiplicity bookkeeping lacks runtime qualification")
    receipts = approval.get("receipts")
    require(isinstance(receipts, list) and receipts, "actual qualification receipts required")
    observed_elements, observed_media, covered_states, observed_multiplicities = set(), set(), {}, set()
    for record in receipts:
        receipt_path = check_pin(record)
        receipt = load_json(receipt_path)
        require(receipt.get("runtime_binding", {}).get("sha256") == runtime_pin["sha256"],
                "control receipt lacks matching current runtime binding")
        require(receipt.get("passed") is True and receipt.get("inputs_after_match") is True, "failed/unaccepted runtime receipt")
        require(receipt.get("binary", {}).get("sha256") == binding["executable"]["sha256"], "receipt executable differs")
        process, output = receipt["process"], receipt["output"]
        require(process.get("exit_code") == 0 and not process.get("error") and not process.get("stop_reason"), "receipt execution failed")
        checks = output.get("checks", {})
        require(output.get("passed") is True and output.get("normal_termination") is True and checks
                and all(v is True for v in checks.values()), "runtime output checks failed")
        require(output.get("accepted") is True and output.get("property_checks")
                and all(v is True for v in output["property_checks"].values()), "current state/population checks failed")
        job = receipt["job"]
        source = job.get("parameter_source")
        # Some real Ca import receipts retain a null historical source descriptor.
        # Actual pinned input bytes + post-run hash + runtime import text below
        # are authoritative; a present descriptor must additionally agree.
        if source is not None:
            require(source.get("sha256") == audit["parameter_sha256"], "receipt imported parameter identity differs")
        files = job["files"]
        xyz_names = [name for name in files if name.endswith(".xyz")]
        require(len(xyz_names) == 1, "receipt requires one pinned XYZ")
        for name in files:
            require(Path(name).name == name, "receipt filename must be local")
            file_path = receipt_path.parent / name
            check_pin({"path": str(file_path), "sha256": files[name]["sha256"]})
        present_elements = set(xyz_elements(receipt_path.parent / xyz_names[0]))
        observed_elements.update(present_elements)
        medium = {"water": "ALPB(Water)", "vacuum": "vacuum"}.get(job.get("medium"))
        require(medium is not None, "unsupported receipt medium")
        observed_media.add(medium)
        multiplicity = output.get("multiplicity")
        require(type(multiplicity) is int and multiplicity > 0, "receipt multiplicity absent")
        observed_multiplicities.add(multiplicity)
        covered_states.setdefault((medium, multiplicity), set()).update(present_elements)
        parameter_names = [name for name, entry in files.items()
                           if entry["sha256"] == audit["parameter_sha256"]]
        require(len(parameter_names) == 1, "receipt lacks exact imported full parameter file")
        parameter_name = parameter_names[0]
        require(receipt["artifacts"].get(parameter_name, {}).get("sha256") == audit["parameter_sha256"],
                "receipt parameter bytes changed after execution")
        output_names = [name for name in receipt["artifacts"] if name.endswith(".out") or name == "stdout.txt"]
        require(len(output_names) == 1, "receipt lacks unambiguous stdout artifact")
        output_name = output_names[0]
        require(Path(output_name).name == output_name, "receipt stdout must be local")
        output_path = receipt_path.parent / output_name
        check_pin({"path": str(output_path), "sha256": receipt["artifacts"][output_name]["sha256"]})
        output_text = output_path.read_text()
        require("ORCA TERMINATED NORMALLY" in output_text, "qualification stdout termination missing")
        runtime_reads = re.findall(r"^\s*Reading XTB parameter file\s+([^\r\n]+?)\s*$", output_text, re.M)
        require(runtime_reads and set(runtime_reads) == {parameter_name}, "qualification runtime import differs")
        lines = output.get("parameter_lines", [])
        require(f"Reading XTB parameter file {parameter_name}" in lines, "receipt lacks runtime import evidence")
    require(set(approval["qualified_elements"]) <= observed_elements, "approval claims unobserved elements")
    require(set(approval["qualified_media"]) <= observed_media, "approval claims unobserved media")
    require(set(approval["qualified_multiplicities"]) <= observed_multiplicities, "approval claims unobserved multiplicities")
    state_key = (binding["model"]["medium"], binding["state"]["multiplicity"])
    require(set(binding["elements"]) <= covered_states.get(state_key, set()),
            "task element/medium/multiplicity combination lacks runtime receipt")
    return {"approval": approval_record, "receipt_count": len(receipts),
            "scope": "qualified parameter/runtime controls; no affinity or structural-transfer claim"}


def stage(binding, parameters, policy, *, approval=None, mode="research", variant="corrected"):
    """Stage byte-exact full parameters beside the already rendered task input.

    Caller controls preparation and execution. This function never starts ORCA.
    Both returned manifest and every source pin must be retained with receipts.
    """
    require(mode in ("research", "diagnostic"), "unknown guard mode")
    require(variant in ("corrected", "original", "native-default"), "unknown parameter variant")
    require(mode == "diagnostic" or variant == "corrected", "historical variants require explicit diagnostic mode")
    binding = load_json(binding) if not isinstance(binding, dict) else binding
    input_path, elements = validate_binding(binding, variant)
    audit = audit_parameters(parameters, policy)
    if variant in ("original", "native-default"):
        require(audit["parameter_sha256"] == ORIGINAL_SHA and not audit["corrections"], "historical variant must pin unchanged original")
    else:
        require(audit["corrections"], "corrected variant requires explicit named corrections")
    require(binding.get("parameter_sha256") == audit["parameter_sha256"], "binding parameter pin differs")
    require(all(e in audit["elements"] for e in elements), "XYZ contains an unknown element")
    qualification = None
    if mode == "research":
        require(all(audit["elements"][e]["parameter_consistent"] for e in elements), "used element has inconsistent shell parameters")
        require(approval is not None, "runtime approval required for research")
        qualification = validate_approval(approval, binding, audit)
    task_dir = input_path.parent
    target, manifest_path = task_dir / _parameter_name(binding), task_dir / LOCAL_MANIFEST
    require(not target.exists() and not manifest_path.exists(), "refuse to overwrite staged parameters/manifest")
    payload = {"schema_version": 1, "mode": mode, "variant": variant, "binding": binding,
               "policy": load_json(policy) if not isinstance(policy, dict) else policy,
               "parameter_audit": audit, "qualification": qualification,
               "guard_implementation": pin(__file__)}
    payload["cache_identity"] = digest(payload)
    if variant != "native-default":
        with target.open("xb") as stream:
            stream.write(Path(parameters).read_bytes())
        require(sha256(target) == audit["parameter_sha256"], "staged parameter digest differs")
    payload["local_parameters"] = {"path": str(target), "sha256": audit["parameter_sha256"]}
    with manifest_path.open("x") as stream:
        json.dump(payload, stream, indent=2, allow_nan=False); stream.write("\n")
    return {"manifest": pin(manifest_path), "cache_identity": payload["cache_identity"],
            "runtime_qualified": mode == "research", "mode": mode, "variant": variant}


def _manifest(path):
    manifest = load_json(path)
    require(manifest.get("schema_version") == 1 and manifest.get("mode") in ("research", "diagnostic"), "unknown guard manifest")
    require(manifest.get("variant") in ("corrected", "original", "native-default"), "unknown parameter variant")
    require(manifest["mode"] == "diagnostic" or manifest["variant"] == "corrected", "historical variant outside diagnostic mode")
    expected = manifest["cache_identity"]
    identity = {k: v for k, v in manifest.items() if k not in ("cache_identity", "local_parameters")}
    require(digest(identity) == expected, "guard manifest/cache identity changed")
    binding = manifest["binding"]
    require(Path(manifest["local_parameters"]["path"]) == Path(binding["input"]["path"]).parent / _parameter_name(binding),
            "local parameter path differs from task input")
    require(manifest["local_parameters"]["sha256"] == binding["parameter_sha256"], "local parameter pin differs")
    check_pin(manifest["guard_implementation"])
    return manifest


def preflight(manifest_path):
    manifest = _manifest(manifest_path)
    binding, variant = manifest["binding"], manifest["variant"]
    validate_binding(binding, variant)
    if variant == "native-default":
        require(manifest["mode"] == "diagnostic", "native defaults are diagnostic only")
        require(not Path(manifest["local_parameters"]["path"]).exists(), "default diagnostic requires absent export before execution")
    else:
        audit = audit_parameters(check_pin(manifest["local_parameters"]), manifest["policy"])
        if manifest["mode"] == "research":
            validate_approval(manifest["qualification"]["approval"], binding, audit)
    return {"parameter_preflight_passed": True, "runtime_qualified": manifest["mode"] == "research",
            "cache_identity": manifest["cache_identity"], "manifest": pin(manifest_path)}


def postflight(manifest_path, stdout, exit_code):
    """Check immutable inputs/import evidence; DOES NOT accept energy or SCF.

    Caller's scientific collector must independently validate convergence,
    gradients, charges, termination, and the complete scoring protocol.
    """
    manifest = _manifest(manifest_path)
    binding, variant = manifest["binding"], manifest["variant"]
    validate_binding(binding, variant)
    check_pin(manifest["local_parameters"])
    audit = audit_parameters(manifest["local_parameters"]["path"], manifest["policy"])
    if manifest["mode"] == "research":
        validate_approval(manifest["qualification"]["approval"], binding, audit)
    require(type(exit_code) is int and exit_code == 0, "execution did not exit successfully")
    text = Path(stdout).read_text()
    require("ORCA TERMINATED NORMALLY" in text, "normal termination absent")
    # Anchored runtime text, excluding the echoed ORCA input (| 12> ...).
    imports = re.findall(r"^\s*Reading XTB parameter file\s+([^\r\n]+?)\s*$", text, re.M)
    if variant != "native-default":
        require(imports and set(imports) == {_parameter_name(binding)}, "actual exact-file runtime import evidence absent/different")
    else:
        require(not imports, "native-default diagnostic unexpectedly imported parameters")
    return {"parameter_postflight_passed": True, "runtime_qualified": manifest["mode"] == "research",
            "scientific_result_accepted": None,
            "scientific_result_note": "caller must separately admit SCF/state/energy/forces",
            "cache_identity": manifest["cache_identity"], "manifest": pin(manifest_path),
            "stdout": pin(stdout), "parameter_file": pin(manifest["local_parameters"]["path"])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    a = sub.add_parser("audit"); a.add_argument("--parameters", required=True); a.add_argument("--policy", required=True)
    sub.add_parser("directives")
    a = sub.add_parser("bind")
    for name in ("input", "xyz", "parameters", "runtime-binding", "protocol-id", "medium"):
        a.add_argument("--"+name, required=True)
    a.add_argument("--source", action="append", required=True)
    a.add_argument("--charge", type=int, required=True)
    a.add_argument("--multiplicity", type=int, required=True)
    a.add_argument("--parameter-filename", default=LOCAL_PARAMETERS)
    a.add_argument("--observed-runtime", help="pinned actual process/library observation supplement")
    a = sub.add_parser("stage")
    for name in ("binding", "parameters", "policy"):
        a.add_argument("--"+name, required=True)
    a.add_argument("--approval", help="runtime registry JSON (hashed automatically)")
    a.add_argument("--mode", choices=("research", "diagnostic"), default="research")
    a.add_argument("--variant", choices=("corrected", "original", "native-default"), default="corrected")
    a = sub.add_parser("preflight"); a.add_argument("--manifest", required=True)
    a = sub.add_parser("postflight"); a.add_argument("--manifest", required=True); a.add_argument("--stdout", required=True); a.add_argument("--exit-code", type=int, required=True)
    args = parser.parse_args()
    try:
        if args.command == "directives":
            print(import_directives(), end=""); return
        if args.command == "audit":
            result = audit_parameters(args.parameters, args.policy)
        elif args.command == "bind":
            result = build_binding(input_path=args.input, xyz=args.xyz, parameters=args.parameters,
                                   runtime_binding=args.runtime_binding, source_inputs=args.source,
                                   protocol_id=args.protocol_id, charge=args.charge, multiplicity=args.multiplicity,
                                   medium=args.medium, parameter_filename=args.parameter_filename,
                                   observed_runtime=args.observed_runtime)
        elif args.command == "stage":
            result = stage(args.binding, args.parameters, args.policy,
                           approval=pin(args.approval) if args.approval else None, mode=args.mode, variant=args.variant)
        elif args.command == "preflight":
            result = preflight(args.manifest)
        else:
            result = postflight(args.manifest, args.stdout, args.exit_code)
        print(json.dumps(result, indent=2, allow_nan=False))
    except (GuardError, KeyError, TypeError, ValueError, OSError) as exc:
        print(json.dumps({"status": "parameter_contract_failed", "reason": str(exc)}), file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
