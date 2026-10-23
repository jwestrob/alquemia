#!/usr/bin/env python3
"""Apply only the reviewed export patch; no build or native invocation."""
import hashlib
import json
from pathlib import Path
import shutil
import time


def pin(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def main():
    start, cpu = time.monotonic(), time.process_time()
    d = Path(__file__).resolve().parent
    root = Path.cwd()
    w = root / d.relative_to(root).as_posix().replace("diagnostics/", "workspaces/", 1)
    inputs = json.loads((d / "PATCH_INPUTS.json").read_text())
    for item in inputs.values():
        if isinstance(item, dict) and set(item) == {"path", "sha256"}:
            assert pin(item["path"]) == item, item
    review = pin(d / "ROOT_REVIEW.json")
    assert review["sha256"] == "ce7cc1cb3f106f4ca953b0bb1cc96b2163f009130fff6e087d36a804bebf64a7"
    parent = json.loads(Path(inputs["parent_manifest"]["path"]).read_text())
    base = Path(parent["source_root"])
    out = w / "source_v2"
    out.mkdir(parents=True, exist_ok=False)
    for item in parent["source_outputs"]:
        assert pin(item["path"]) == item, item
        dst = out / Path(item["path"]).relative_to(base)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(item["path"], dst)
    source = out / "source/induce.f"
    old, anchor = source.read_text(), inputs["hook_anchor"]
    assert old.count(anchor) == 1
    source.write_text(old.replace(anchor, anchor.replace("      call fftback", "      call pmenormexport\n      call fftback")) + Path(inputs["helper"]["path"]).read_text())
    assert pin(source)["sha256"] == inputs["prospective_patched_source_sha256"]
    shutil.copy2(inputs["patch"]["path"], out / "PME_NORM_EXPORT.patch")
    outputs = [pin(out / Path(v["path"]).relative_to(base)) for v in parent["source_outputs"]]
    changed = [Path(a["path"]).name for a, b in zip(outputs, parent["source_outputs"]) if a["sha256"] != b["sha256"]]
    assert changed == ["induce.f"], changed
    assert pin(out / "cmake/CMakeLists.txt") in outputs
    manifest = dict(status="reviewed_PME_export_patch_applied_source_only", source_root=str(out),
                    parent_manifest=inputs["parent_manifest"], export_ready=pin(d / "READY.json"),
                    root_source_review=review, patch_inputs=pin(d / "PATCH_INPUTS.json"),
                    preparation_source=pin(__file__), patch=pin(out / "PME_NORM_EXPORT.patch"),
                    changed_files=changed, source_outputs=outputs, new_native_calls=0,
                    periodic_response_guards_remain_fatal=True,
                    cost=dict(wall_seconds=time.monotonic()-start, CPU_seconds=time.process_time()-cpu, threads=1))
    (out / "MANIFEST.json").write_text(json.dumps(manifest, indent=2)+"\n")
    # This existing recipe is already single-worker and covers the actual manifest schema.
    builder_src = Path(inputs["parent_ready"]["path"]).parent / "BUILD_LOGIN.py"
    builder = d / "BUILD_LOGIN_ONE.py"
    assert not builder.exists()
    shutil.copy2(builder_src, builder)
    ready = dict(status="ready_for_one_authorized_single_worker_build", source_manifest=pin(out / "MANIFEST.json"),
                 builder=pin(builder), parent_builder=pin(builder_src), root_source_review=review,
                 build_output=str(w / "build_v1/execution"), build_count_cap=1, compiler_workers=1,
                 native_calls=0, qualification_execution_authorized=False,
                 argv=["/home/jwestrob/.pyenv/versions/3.12.0/bin/python3", str(builder), str(out / "MANIFEST.json"),
                       pin(out / "MANIFEST.json")["sha256"], str(w / "build_v1/execution")])
    (d / "BUILD_READY.json").write_text(json.dumps(ready, indent=2)+"\n")
    print(json.dumps(ready, indent=2))


if __name__ == "__main__":
    main()
