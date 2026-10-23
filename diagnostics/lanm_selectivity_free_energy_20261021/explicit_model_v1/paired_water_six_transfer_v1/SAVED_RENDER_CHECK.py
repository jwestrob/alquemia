"""Actual consumed-input rendering checks; zero native calls and no candidate."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time


def pin(p):
    return {"path": str(p.resolve()), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


def module(p, name):
    spec = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    start, cpu = time.monotonic(), time.process_time()
    assert all(os.environ.get(k) == "1" for k in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"))
    here = Path(__file__).resolve().parent
    root = here.parents[3]
    workspace = root / "workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1"
    old_path = here.parent / "paired_water_continuation_transfer_v1/RUN.py"
    assert pin(old_path)["sha256"] == "6380231a2a060add2fbb9812bcb7b50f626e9b2ac6499b5bd8cecf042a50f0fc"
    old = module(old_path, "saved_six_render_old")
    new = module(here / "RENDER.py", "saved_six_render_new")
    data_pin = old.PINS["pins"]["prepared32"]
    data = old.loaded(old.verify(data_pin))
    epsilon_path = workspace / "paired_water_epsilon_diagnostic_v1/execution_v1/run_v1/execution/RESULT.json"
    assert pin(epsilon_path)["sha256"] == "a93f1502a82ee41f34cd9ff40999b37100de514ff3c19d07b547dfcf70396d80"
    epsilon = json.loads(epsilon_path.read_text())
    manifest = old.read(epsilon["manifest"])
    values = manifest["baseline"]["parameters"]
    eps = [old.d.e.EPSILON[m] for m in ("La", "Dy")]
    checks = []

    def check(name, condition):
        checks.append({"name": name, "passed": bool(condition)})
        assert condition, name

    for gid, row in data["rows"].items():
        old_key, old_render = old.d.render(data["config"], row, values)
        key, rendered = new.render(old.d, data["config"], row, values + eps)
        check("Exact original-epsilon key " + gid, key == old_key)
        check("Exact original pair/recipient " + gid, all(
            rendered[k] == old_render[k] for k in (
                "water_recipient", "original_key_reused", "resolved_water_pair")))
    for probe in manifest["probes"]:
        v = values + eps
        v[4 if probe["metal"] == "La" else 5] = probe["epsilon"]
        for gid, inp in probe["inputs"].items():
            key, _ = new.render(old.d, data["config"], data["rows"][gid], v)
            check("Exact measured epsilon input " + probe["id"] + " " + gid,
                  key == old.verify(inp["files"]["key"]).read_text())
        for gid in old.d.DRY:
            row = data["rows"][gid]
            key, rendered = new.render(old.d, data["config"], row, v)
            check("Dry input unchanged " + probe["id"] + " " + gid,
                  key == old.verify(row["input"]["key"]).read_text()
                  and rendered["original_key_reused"])
    for name, bad in (("zero epsilon", values + [0., eps[1]]),
                      ("nonfinite epsilon", values + [float("nan"), eps[1]]),
                      ("above declared bound", values + [1.51 * eps[0], eps[1]])):
        rejected = False
        try:
            new.render(old.d, data["config"], data["rows"][old.d.CENTERS[0]], bad)
        except ValueError:
            rejected = True
        check(name, rejected)
    out = workspace / "paired_water_six_transfer_v1/render_controls_v1"
    out.mkdir(parents=True, exist_ok=False)
    result = {"status": "passed_saved_render_checks", "implementation": pin(here / "RENDER.py"),
              "control_source": pin(Path(__file__)), "old_adapter": pin(old_path),
              "prepared32": data_pin, "actual_epsilon_result": pin(epsilon_path),
              "checks": checks, "new_native_calls": 0, "candidate_bound": False,
              "elapsed_seconds": time.monotonic() - start,
              "process_cpu_seconds": time.process_time() - cpu}
    (out / "RESULT.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"result": pin(out / "RESULT.json"), "checks": len(checks),
                      "elapsed_seconds": result["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
