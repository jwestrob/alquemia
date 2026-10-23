#!/usr/bin/env python3
"""Report a terminal donor fit from saved results; no evaluations or selection."""
import argparse
import hashlib
import json
from pathlib import Path
import time


def pin(path):
    p = Path(path).resolve()
    return {"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}


def summarize(score):
    energy = {
        name: {"count": len(values), "failed": sum(abs(v) > 1 for v in values.values()),
               "max_abs_kcal_mol": max(abs(v) for v in values.values())}
        for name, values in score["energy_errors_kcal_mol"].items()
    }
    force = {
        kind: {"count": len(values), "failed": sum(not v["passed"] for v in values.values()),
               "max_fixture_RMS_kcal_mol_A": max(v["rigid_RMS_kcal_mol_A"] for v in values.values()),
               "max_component_kcal_mol_A": max(v["maximum_body_component_error_kcal_mol_A"] for v in values.values())}
        for kind, values in score["per_fixture"].items()
    }
    return {"objective": score["objective"], "blocks": score["blocks"],
            "energy": energy, "forces": force,
            "all_development_checks_passed": score["development_checks_passed"]}


def main(result_path, output):
    started, cpu = time.monotonic(), time.process_time()
    result_path = Path(result_path).resolve()
    result = json.loads(result_path.read_text())
    if result["status"] != "bounded_paired_donor_complete" or result["error"] is not None:
        raise ValueError("A complete, error-free terminal donor result is required")
    if not result["candidate_frozen"] or not result["original_results_unchanged"]:
        raise ValueError("Unfrozen selection or modified inputs")
    before, after = result["initial"], result["selected"]
    initial, selected = summarize(before["score"]), summarize(after["score"])
    parameters = []
    for definition, a, b in zip(result["family"]["parameters"], before["parameters"], after["parameters"]):
        lo, hi = definition["lower"], definition["upper"]
        parameters.append({"definition": definition, "initial": a, "selected": b,
                           "fraction_of_interval": (b - lo) / (hi - lo),
                           "distance_to_nearest_bound_fraction": min(b - lo, hi - b) / (hi - lo)})
    iterations = result["iterations"]
    rank = None
    if iterations and "linear_directions" in iterations[-1]:
        rank = {"origin_key": iterations[-1]["origin_key"],
                "is_selected_origin": result.get("current_key") == iterations[-1]["origin_key"],
                "diagnostic": iterations[-1]["linear_directions"],
                "condition": iterations[-1].get("origin_condition")}
    report = {"status": "saved_terminal_summary", "source": pin(result_path),
              "reporter": pin(__file__), "initial": initial, "selected": selected,
              "stopping_reason": result["stopping_reason"], "parameters": parameters,
              "completed_outer_steps": result["completed_outer_steps"],
              "new_parameter_vectors": result["new_parameter_vectors"],
              "new_native_program_invocations": result["new_native_program_invocations"],
              "last_derivative_origin": rank,
              "frozen_water_invariants": result["frozen_water_invariants"],
              "new_molecular_calls_for_report": 0,
              "interpretation": "Consumed-data development; neither independent validation nor binding free energy."}
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    labels = {"balanced_exchange": "Aquo-referenced exchange",
              "water_addition": "Wet minus dry exchange",
              "individual_radial_work": "Individual radial work",
              "paired_radial_work": "Paired radial work (control only)"}
    lines = ["# Joint donor fit: saved result", "",
             f"**All development controls pass: {selected['all_development_checks_passed']}.** "
             f"Stop: `{result['stopping_reason']}`.", "",
             f"Objective {initial['objective']:.6g} → {selected['objective']:.6g}; "
             f"{result['completed_outer_steps']} completed steps, {result['new_parameter_vectors']} new vectors, "
             f"{result['new_native_program_invocations']} native program invocations.", "",
             "## Energy controls", "", "| Control | Failed, initial → selected | Maximum absolute error, initial → selected (kcal/mol) |",
             "|---|---:|---:|"]
    for key, a in initial["energy"].items():
        b = selected["energy"][key]
        lines.append(f"| {labels[key]} | {a['failed']}/{a['count']} → {b['failed']}/{b['count']} | {a['max_abs_kcal_mol']:.4f} → {b['max_abs_kcal_mol']:.4f} |")
    lines.extend(["", "## Rigid-body force controls", "", "| Block | Failed fixtures, initial → selected | Worst fixture RMS, initial → selected (kcal/mol/Å) |",
                  "|---|---:|---:|"])
    for key, a in initial["forces"].items():
        b = selected["forces"][key]
        lines.append(f"| {key} | {a['failed']}/{a['count']} → {b['failed']}/{b['count']} | {a['max_fixture_RMS_kcal_mol_A']:.4f} → {b['max_fixture_RMS_kcal_mol_A']:.4f} |")
    lines.extend(["", "Every energy error must be ≤1 kcal/mol; every fixture requires RMS ≤2 and maximum component ≤5 kcal/mol/Å. Water parameters and results remain frozen. An average improvement does not establish these controls.", "",
                  "This uses consumed development geometries and provisional PBE0-D4 targets. It does not establish independent physical accuracy, a global optimum, or a protein exchange free energy. Any reported Jacobian applies to its recorded origin, which may precede the selected point.", "",
                  f"Source: `{result_path}` (SHA256 `{report['source']['sha256']}`).", ""])
    (out / "REPORT.md").write_text("\n".join(lines))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    fig, axes = plt.subplots(1, 3, figsize=(16, 6.7), gridspec_kw={"width_ratios": [1.25, 1, .9]})
    energies = after["score"]["energy_errors_kcal_mol"]["balanced_exchange"]
    ids = list(energies)
    y = np.arange(len(ids))
    short = [k.replace("hydrated_", "wet ").replace("matched_dry_", "dry ").replace("acetate_", "acetate ").replace("_central", " center").replace("_minus", " inward").replace("_plus", " outward") for k in ids]
    ax = axes[0]
    ax.axvspan(-1, 1, color="#ddebd9", label="±1 criterion")
    ax.scatter([before["score"]["energy_errors_kcal_mol"]["balanced_exchange"][k] for k in ids], y, color="#777777", marker="x", label="Initial")
    ax.scatter([energies[k] for k in ids], y, color="#1665ad", label="Selected")
    ax.set(yticks=y, yticklabels=short, xlabel="Exchange error (kcal/mol)", title="Each donor environment")
    ax.invert_yaxis(); ax.axvline(0, color="black", lw=.5); ax.legend(fontsize=8)
    ax = axes[1]
    grid = np.array([[after["score"]["per_fixture"][kind][k]["rigid_RMS_kcal_mol_A"] for kind in ("La", "Dy", "paired")] for k in ids])
    im = ax.imshow(grid, aspect="auto", cmap="YlOrRd", vmin=0)
    for i in range(len(ids)):
        for j in range(3):
            rec = after["score"]["per_fixture"][("La", "Dy", "paired")[j]][ids[i]]
            ax.text(j, i, f"{grid[i,j]:.1f}" + ("*" if not rec["passed"] else ""), ha="center", va="center", fontsize=8)
    ax.set(xticks=range(3), xticklabels=["La", "Dy", "Dy−La"], yticks=y, yticklabels=[], title="Selected rigid-force RMS")
    fig.colorbar(im, ax=ax, shrink=.75, label="kcal/mol/Å (* = failed fixture)")
    ax = axes[2]
    losses = [before["score"]["objective"]] + [result["vectors"][row["selected_key"]]["score"]["objective"] for row in iterations if row["status"] == "complete"]
    ax.plot(range(len(losses)), losses, "o-", color="#1665ad", ms=3)
    ax.set(xlabel="Completed step", ylabel="Declared weighted objective", title="Actual selected objective")
    ax.grid(alpha=.2)
    fig.suptitle("Shared donor parameters with the passing water model fixed", fontsize=14)
    fig.text(.5, .015, "Consumed PBE0-D4 development data • one shared carboxylate pair covers mono- and bidentate acetate • no affinity claim", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .04, 1, .94))
    for extension in ("png", "pdf"):
        fig.savefig(out / ("donor_response." + extension), dpi=180)
    plt.close(fig)
    report["report_wall_seconds"] = time.monotonic() - started
    report["report_process_CPU_seconds"] = time.process_time() - cpu
    (out / "SUMMARY.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"summary": pin(out / "SUMMARY.json"), "report": pin(out / "REPORT.md")}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    main(args.result, args.output)
