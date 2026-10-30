#!/usr/bin/env python3
"""Read one completed, fresh CP2K-2026.2 ALMO FULL_SCF quartet output.

Usage (postprocessing only; this program never launches CP2K):
  python PARSE.py --output native.out --geometry-id aquo8_central_La \
      --exit-code 0 --result parsed.json

ENERGY_LEDGER.json must be adjacent and match LEDGER_SHA256. Result creation is
exclusive: choose a new result filename for each actual attempt. CLI exits 0
only when all gates pass, otherwise 1 after saving the failure and raw evidence.
Missing/ambiguous native records fail closed. Input/build/geometry identity and
root runtime admission remain separate: output parsing cannot certify them.

Source anchors in the exact downloaded 2026.2 release:
  mscfg_methods.F:376-414: fragment delimiter, charge/multiplicity, atom list.
  qs_scf_loop_utils.F:595-616: native fragment convergence/failure messages.
  qs_scf_output.F:167-174,571-581,873: counts, iteration row, total energy.
  almo_scf.F:790-793,970-1184: R/F, settings and domain table.
  almo_scf_optimizer.F:242-250,404-426: DIIS row and fatal nonconvergence.
  almo_scf_optimizer.F:1492,1855-1858,1973: PCG row/fatal nonconvergence.
  force_env_methods.F:386-391: final ENERGY| record (15 decimal places).
The FULL_SCF branch does not require perturbative CORRECTED ENERGY output.
Native convergence is checked separately from the telescoping energy identity.
"""

import argparse
import datetime
from decimal import Decimal, InvalidOperation, localcontext
import hashlib
import json
from pathlib import Path
import re
import sys


LEDGER_SHA256 = "24c699cc97c3c55857c56a5ee07dda6213eeaf12e59af69395f23bbd05b732c5"
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
AGREEMENT = Decimal("1e-4")
CLOSURE = Decimal("1e-8")
ORDER = Decimal("1e-8")
FRAGMENT_SCF = Decimal("1e-10")
ALMO_SCF = Decimal("1e-9")
# Exact converted spherical orbital basis counts; DECK_REVIEW.md records the
# atom/kind counts. These are independent expected values, not self-comparisons
# between two native reports of a potentially incorrect basis.
EXPECTED_AOS = {
    "aquo8_central_La": (47, 344),
    "aquo8_central_Dy": (47, 344),
    "hydrated_acetate_mono_central_La": (47, 443),
    "hydrated_acetate_mono_central_Dy": (47, 443),
}


def decimal(token):
    value = Decimal(token.replace("D", "E").replace("d", "e"))
    if not value.is_finite():
        raise ValueError("Nonfinite numeric token")
    return value


def pin(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def half_unit(token):
    """Printing resolution only; never an SCF/error estimate."""
    return Decimal(5).scaleb(decimal(token).as_tuple().exponent - 1)


class ParseFailure(ValueError):
    pass


class Output:
    def __init__(self, text, result):
        self.lines = text.splitlines()
        self.result = result

    def matches(self, pattern, lo=0, hi=None):
        regex = re.compile(pattern)
        return [{"line": i + 1, "raw": self.lines[i], "tokens": list(m.groups())}
                for i in range(lo, len(self.lines) if hi is None else hi)
                if (m := regex.fullmatch(self.lines[i]))]

    def one(self, name, pattern, lo=0, hi=None):
        hits = self.matches(pattern, lo, hi)
        self.result["evidence"][name] = hits
        if len(hits) != 1:
            raise ParseFailure(f"{name}: expected one record, observed {len(hits)}")
        return hits[0]

    def gate(self, name, passed, **evidence):
        self.result["gates"][name] = {"passed": bool(passed), **evidence}


def native_rows(out, name, lo, hi, allowed):
    pattern = (rf"\s*ALMO SCF (DIIS|CG|LS)\s+(\d+)\s+({NUMBER})\s+"
               rf"({NUMBER})\s+({NUMBER})\s+({NUMBER})\s*")
    rows = out.matches(pattern, lo, hi)
    out.result["evidence"][name] = rows
    candidates = [i + 1 for i in range(lo, hi)
                  if re.match(r"\s*ALMO SCF\s", out.lines[i])]
    if not rows or [r["line"] for r in rows] != candidates:
        raise ParseFailure(f"{name}: absent or unrecognized native iteration row")
    if any(r["tokens"][0] not in allowed for r in rows):
        raise ParseFailure(f"{name}: unexpected optimizer")
    for row in rows:
        for token in row["tokens"][2:]:
            decimal(token)
    steps = [int(r["tokens"][1]) for r in rows]
    out.gate(name + "_iteration_sequence", steps[0] in (0, 1) and
             all(b == a + 1 for a, b in zip(steps, steps[1:])) and steps[-1] <= 300,
             first=steps[0], last=steps[-1], row_count=len(rows))
    final = rows[-1]
    residual = decimal(final["tokens"][4])
    out.gate(name + "_printed_residual", Decimal(0) <= residual <= ALMO_SCF,
             value=final["tokens"][4], threshold=str(ALMO_SCF),
             rounding_half_unit=str(half_unit(final["tokens"][4])),
             note="Printed F14.9 error alone is insufficient: normal completion, no early-stop/failure, and correct native branch are also required.")
    return final


def parse(output_path, geometry_id, exit_code):
    result = {"status": "failed_cp2k_almo_output_checks", "geometry_id": geometry_id,
              "exit_code": exit_code, "implementation": pin(__file__),
              "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "gates": {}, "evidence": {}, "error": None,
              "limits": ["Output checks do not replace exact input, basis/ECP, build and geometry admission.",
                         "FRZ is combined frozen interaction; electrostatic/Pauli split is unavailable.",
                         "Printed precision and arithmetic closure are not SCF or physical uncertainty.",
                         "A single system does not establish a quartet double difference."]}
    try:
        ledger_path = Path(__file__).resolve().with_name("ENERGY_LEDGER.json")
        result["ledger"] = pin(ledger_path)
        if result["ledger"]["sha256"] != LEDGER_SHA256:
            raise ParseFailure("Adjacent energy ledger hash mismatch")
        ledger = json.loads(ledger_path.read_text())
        if geometry_id not in ledger["systems"]:
            raise ParseFailure("Geometry is outside the exact archived quartet")
        archived = ledger["systems"][geometry_id]
        state = archived["geometry_state"]
        result["archive"] = archived["provenance"]
        xyz_pin = archived["provenance"]["xyz"]
        if pin(xyz_pin["path"]) != xyz_pin:
            raise ParseFailure("Archived atom-identity source hash mismatch")
        xyz = Path(xyz_pin["path"]).read_text().splitlines()
        symbols = [line.split()[0] for line in xyz[2:]]
        if int(xyz[0]) != state["atoms"] or len(symbols) != state["atoms"]:
            raise ParseFailure("Archived XYZ record count mismatch")
        raw = Path(output_path).read_bytes()
        result["output"] = {"path": str(Path(output_path).resolve()),
                            "sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}
        out = Output(raw.decode("utf-8", errors="strict"), result)
        out.gate("native_exit_code_zero", exit_code == 0)
        start = out.one("program_start", r".*PROGRAM STARTED AT\s+(.+)")
        end = out.one("program_end", r".*PROGRAM ENDED AT\s+(.+)")
        failures = [{"line": i + 1, "raw": line} for i, line in enumerate(out.lines)
                    if re.search(r"\bABORT\b|SCF run NOT converged|not converged!|FAILED to converge|"
                                 r"Leaving inner SCF loop|EPS_ERROR_EARLY is on|MAX_ITER_EARLY is on|"
                                 r"Read restart ALMO|ASPC order:", line, re.I)]
        result["evidence"]["native_failures_or_unapproved_stopping"] = failures
        out.gate("native_no_failure_or_early_stop", not failures)

        R = out.one("R", rf"\s*Single-molecule energy:\s+({NUMBER})\s*")
        F = out.one("F", rf"\s*Energy of the initial guess:\s+({NUMBER})\s*")
        block = out.one("block_heading", r"\s*-+ Optimization of block-diagonal ALMOs -+\s*")
        full = out.one("full_heading", r"\s*-+ Optimization of fully delocalized MOs -+\s*")
        U = out.one("U", rf"\s*ENERGY\| Total FORCE_EVAL \(\s*QS\s*\) energy \[a\.u\.\]\s+({NUMBER})\s*")
        order = [start["line"], R["line"], F["line"], block["line"], full["line"], U["line"], end["line"]]
        if not all(b > a for a, b in zip(order, order[1:])):
            raise ParseFailure("Native energy/stage records are out of order")
        out.gate("native_stage_order", True, lines=order)
        deloc = out.one("delocalization", r"\s*delocalization:\s+(\S+)\s*")
        radius = out.one("delocalization_radius", r"\s*delocalization cutoff radius:\s+(\S+)\s*")
        out.gate("native_full_scf", deloc["tokens"] == ["FULL_SCF"] and radius["tokens"] == ["infinite"])
        domains = out.one("total_fragments", r"\s*Total fragments:\s+(\d+)\s*")
        out.gate("two_domains", domains["tokens"] == ["2"])
        bsettings = out.one("block_settings", r"\s*optimization of block-diagonal ALMOs:\s*")
        fsettings = out.one("full_settings", r"\s*optimization of extended orbitals:\s*")
        for name, low, high, method in [("block", bsettings["line"], deloc["line"] - 1, "DIIS"),
                                        ("full", fsettings["line"], domains["line"] - 1, "PCG")]:
            opt = out.one(name + "_optimizer", r"\s*optimizer type:\s+(\S+)\s*", low, high)
            target = out.one(name + "_target", rf"\s*target error:\s+({NUMBER})\s*", low, high)
            maximum = out.one(name + "_maximum", r"\s*maximum iterations:\s+(\d+)\s*", low, high)
            out.gate(name + "_native_settings", opt["tokens"] == [method] and
                     decimal(target["tokens"][0]) == ALMO_SCF and maximum["tokens"] == ["300"])
            if name == "full":
                outer = out.one("full_outer_maximum", r"\s*maximum outer loop iterations:\s+(\d+)\s*", low, high)
                out.gate("no_outer_iteration_extension", outer["tokens"] == ["0"])

        table = out.one("domain_table_heading", r"\s*Fragment\s+Basis Set\s+Occupied\s+Virtual\s+Charge\s+Deloc Neig\s*")
        table_rows = out.matches(r"\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(-?\d+)\s+(ALL|NONE|N/A|\d+)\s*",
                                 table["line"], table["line"] + 5)
        result["evidence"]["domain_table"] = table_rows
        if len(table_rows) != 2 or [int(t["tokens"][0]) for t in table_rows] != [1, 2]:
            raise ParseFailure("Missing/unknown two-fragment domain table")
        expected_aos = EXPECTED_AOS[geometry_id]
        for i, row in enumerate(table_rows):
            _, ao, occ, virt, charge, neighbor = row["tokens"]
            expected = state["fragments"][i]
            out.gate(f"domain_{i+1}_state", int(ao) == expected_aos[i] and int(occ)*2 == expected["explicit_electrons"] and
                     int(virt) == int(ao)-int(occ) and int(charge) == expected["charge"] and neighbor == "ALL",
                     native_AOs=int(ao), expected_AOs=expected_aos[i], native_occupied=int(occ), native_virtual=int(virt), native_charge=int(charge))
        out.gate("parent_AO_count_from_domains", sum(int(t["tokens"][1]) for t in table_rows) == sum(expected_aos),
                 native=sum(int(t["tokens"][1]) for t in table_rows), expected=sum(expected_aos))

        fragment_headers = out.matches(r"\s*-\s*MOLECULAR GUESS:\s*FRAGMENT\s+(\d+)\s+OUT OF\s+(\d+)\s*-\s*")
        result["evidence"]["fragment_headers"] = fragment_headers
        if [h["tokens"] for h in fragment_headers] != [["1", "2"], ["2", "2"]]:
            raise ParseFailure("Expected exactly two ordered native molecular-guess fragments")
        fragment_energies = []
        for i, header in enumerate(fragment_headers):
            lo = header["line"]
            hi = fragment_headers[i+1]["line"]-1 if i == 0 else R["line"]-1
            expected = state["fragments"][i]
            charge = out.one(f"fragment_{i+1}_charge_multiplicity",
                             r"\s*-\s*CHARGE\s+(-?\d+)\s+MULTIPLICITY\s+(\d+)\s*-\s*", lo, hi)
            out.gate(f"fragment_{i+1}_state", list(map(int, charge["tokens"])) ==
                     [expected["charge"], expected["effective_multiplicity"]])
            atoms = out.matches(r"\s*-\s+(\d+)\s+([A-Za-z][A-Za-z0-9_]*)\s*-\s*", lo, hi)
            result["evidence"][f"fragment_{i+1}_atoms"] = atoms
            out.gate(f"fragment_{i+1}_membership", [[int(a["tokens"][0]), a["tokens"][1]] for a in atoms] ==
                     [[n, symbols[n-1]] for n in expected["atom_indices_one_based"]])
            converged = out.one(f"fragment_{i+1}_converged",
                                r"\s*\*\*\* SCF run converged in\s+(\d+) steps \*\*\*\s*", lo, hi)
            energy = out.one(f"fragment_{i+1}_energy", rf"\s*Total energy:\s+({NUMBER})\s*", converged["line"], hi)
            fragment_energies.append(decimal(energy["tokens"][0]))
            # qs_scf_output.F: the five numeric suffixes are parameter, time,
            # convergence, total energy, energy change. Method may contain spaces.
            scf_rows = out.matches(rf"\s*(\d+)\s+(OT\s+(?:CG|LS|SD|DIIS))\s+({NUMBER})\s+({NUMBER})\s+"
                                   rf"({NUMBER})\s+({NUMBER})\s+({NUMBER})\s*", lo, converged["line"]-1)
            if not scf_rows:
                raise ParseFailure(f"fragment {i+1}: no recognized final ordinary OT SCF row")
            last = scf_rows[-1]
            result["evidence"][f"fragment_{i+1}_last_scf_row"] = last
            out.gate(f"fragment_{i+1}_native_convergence", last["tokens"][0] == converged["tokens"][0] and
                     1 <= int(converged["tokens"][0]) <= 300 and
                     Decimal(0) <= decimal(last["tokens"][4]) <= FRAGMENT_SCF,
                     residual_token=last["tokens"][4], threshold=str(FRAGMENT_SCF),
                     native_message_required=True, residual_rounding_half_unit=str(half_unit(last["tokens"][4])))
            for label, expected_count in [("electrons", expected["explicit_electrons"]),
                                           ("occupied orbitals", expected["explicit_electrons"]//2),
                                           ("orbital functions", expected_aos[i])]:
                counts = out.matches(rf"\s*Number of {label}:\s+(\d+)\s*", lo, converged["line"]-1)
                result["evidence"][f"fragment_{i+1}_{label}"] = counts
                # Domain table and native fragment state are mandatory even when
                # ordinary QS count printing is disabled; no invented count.
                out.gate(f"fragment_{i+1}_{label}_if_emitted", not counts or
                         (len(counts) == 1 and int(counts[0]["tokens"][0]) == expected_count),
                         available=bool(counts), expected=expected_count)

        p = native_rows(out, "block_diis", block["line"], full["line"]-1, {"DIIS"})
        utrace = native_rows(out, "full_pcg", full["line"], U["line"]-1, {"CG", "LS"})
        result["convergence_interpretation"] = (
            "Two explicit ordinary-SCF convergence messages plus their final residuals are required. "
            "ALMO has no separate success sentence: the final DIIS/PCG residuals, declared target/iteration settings, "
            "absence of early stopping/abort, and normal program end together establish its native termination. "
            "The inspected source aborts on nonconvergence. None of these is inferred from FRZ+POL+CT arithmetic.")
        values = {"R": decimal(R["tokens"][0]), "F": decimal(F["tokens"][0]),
                  "P": decimal(p["tokens"][2]), "U": decimal(U["tokens"][0])}
        components = {"FRZ": values["F"]-values["R"], "POL": values["P"]-values["F"],
                      "CT": values["U"]-values["P"]}
        interaction = values["U"]-values["R"]
        closure = sum(components.values())-interaction
        ref_residual = sum(fragment_energies)-values["R"]
        full_trace_residual = decimal(utrace["tokens"][2])-values["U"]
        out.gate("fragment_reference_sum", abs(ref_residual) <= CLOSURE, residual_hartree=str(ref_residual), tolerance_hartree=str(CLOSURE))
        out.gate("full_trace_matches_final_energy", abs(full_trace_residual) <= CLOSURE,
                 residual_hartree=str(full_trace_residual), tolerance_hartree=str(CLOSURE))
        out.gate("component_arithmetic_closure", abs(closure) <= CLOSURE,
                 residual_hartree=str(closure), tolerance_hartree=str(CLOSURE),
                 note="Telescoping bookkeeping identity, not independent convergence evidence.")
        out.gate("variational_block_not_above_frozen", values["P"]-values["F"] <= ORDER,
                 difference_hartree=str(values["P"]-values["F"]), tolerance_hartree=str(ORDER))
        out.gate("variational_full_not_above_block", values["U"]-values["P"] <= ORDER,
                 difference_hartree=str(values["U"]-values["P"]), tolerance_hartree=str(ORDER))
        archive_energy = archived["energy_hartree"]
        archived_el = decimal(archive_energy["electronic_after_final_cosx_integration"]["token"])
        archived_full = decimal(archive_energy["full_pbe0_d4"]["token"])
        d4 = decimal(archive_energy["d4_printed"]["token"])
        electronic_error = values["U"]-archived_el
        total_error = values["U"]+d4-archived_full
        out.gate("archived_electronic_agreement", abs(electronic_error) <= AGREEMENT,
                 difference_hartree=str(electronic_error), tolerance_hartree=str(AGREEMENT))
        out.gate("archived_total_agreement", abs(total_error) <= AGREEMENT,
                 difference_hartree=str(total_error), tolerance_hartree=str(AGREEMENT),
                 note="Same archived directly printed D4 added exactly once; no tolerance padding.")
        factor = decimal(ledger["units"]["hartree_to_kcal_mol_literal_from_original_collector"])
        result["energies_hartree"] = {k: str(v) for k, v in values.items()}
        result["components_hartree"] = {k: str(v) for k, v in components.items()}
        result["components_kcal_mol"] = {k: str(v*factor) for k, v in components.items()}
        result["interaction_electronic_hartree"] = str(interaction)
        result["dispersion_bookkeeping"] = {
            "archived_D4_hartree": str(d4), "full_PBE0_D4_hartree": str(values["U"]+d4),
            "CT_suppressed_PBE0_D4_hartree": str(values["P"]+d4),
            "interaction_relative_to_electronic_fragments_hartree": str(interaction+d4),
            "true_D4_corrected_isolated_fragment_interaction_hartree": None,
            "limitation": "Isolated whole-shell D4 reference unavailable; cancels in matched quartet DD under ledger conditions. POL/CT remain electronic."}
        result["printing"] = {"R_half_unit_hartree": str(half_unit(R["tokens"][0])),
                              "F_half_unit_hartree": str(half_unit(F["tokens"][0])),
                              "P_half_unit_hartree": str(half_unit(p["tokens"][2])),
                              "U_half_unit_hartree": str(half_unit(U["tokens"][0])),
                              "note": "No estimate of actual SCF, discretization, basis or physical error."}
        result["status"] = ("passed_cp2k_almo_output_checks" if all(g["passed"] for g in result["gates"].values())
                            else "failed_cp2k_almo_output_checks")
    except (OSError, UnicodeError, ValueError, KeyError, TypeError, InvalidOperation) as error:
        result["error"] = f"{type(error).__name__}: {error}"
        result["gates"]["complete_recognized_evidence"] = {"passed": False}
    result["passed"] = result["status"] == "passed_cp2k_almo_output_checks"
    result["new_molecular_calls"] = 0
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--geometry-id", required=True)
    ap.add_argument("--exit-code", required=True, type=int)
    ap.add_argument("--result", required=True, type=Path)
    args = ap.parse_args()
    if args.result.exists():
        ap.error("Result already exists; preserve the prior receipt and choose a fresh filename")
    with localcontext() as ctx:
        ctx.prec = 50
        result = parse(args.output, args.geometry_id, args.exit_code)
    args.result.parent.mkdir(parents=True, exist_ok=True)
    with args.result.open("x") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"status": result["status"], "result": pin(args.result), "error": result["error"]}))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
