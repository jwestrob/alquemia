"""Stage A inputs: Dy_ssaq_w21_s{1,2} = metal swap La->Dy in the final geometry of the finished La_ssaq_w21_s{1,2}.

job.inp = the La job.inp with only the metal %basis block replaced by the Dy block the project already
generated for the cancelled run_ss/Dy_ssaq_w21_s{1,2}, and `nprocs 24` written as `nprocs @NPROCS@`
(orca_job.sh fills it from $SLURM_CPUS_ON_NODE). Checks are asserted; BUILD.json records hashes and CN.
Login node, single thread, a few files. Run once; refuses to overwrite an existing run directory.
"""
import difflib, hashlib, json, re
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "workspaces/lady_relaxed_exchange_20261024/run_ss"
OUT = REPO / "workspaces/lady_second_shell_20261109/runs"
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def split(txt):
    """header (up to %basis), basis block (%basis .. line before '* xyzfile'), tail (from '* xyzfile')."""
    i = txt.index("%basis"); j = txt.index("* xyzfile")
    return txt[:i], txt[i:j], txt[j:]


def cn_report(xyz_lines):
    el = [l.split()[0] for l in xyz_lines]; X = np.array([[float(v) for v in l.split()[1:4]] for l in xyz_lines])
    d = np.linalg.norm(X - X[0], axis=1)
    o = sorted(float(x) for e, x in zip(el, d) if e == "O")
    return {"CN_lt2p95": sum(x < 2.95 for x in o), "O_2p95_3p10": [round(x, 3) for x in o if 2.95 <= x < 3.10],
            "nearest_O": [round(x, 3) for x in o[:11]]}


def main():
    meta = {}
    for s in (1, 2):
        la, dy_old = SRC / f"La_ssaq_w21_s{s}", SRC / f"Dy_ssaq_w21_s{s}"
        t = (la / "job.out").read_text(errors="replace")
        assert "ORCA TERMINATED NORMALLY" in t and "Final Gibbs free energy" in t, la
        la_inp = (la / "job.inp").read_text(); dy_inp_old = (dy_old / "job.inp").read_text()
        h_la, b_la, t_la = split(la_inp); h_dy, b_dy, t_dy = split(dy_inp_old)
        norm = lambda x: re.sub(r"nprocs \S+ end", "nprocs @NPROCS@ end", x)
        assert norm(h_la) == norm(h_dy) and t_la == t_dy, "old Dy input differs from La outside the basis block"
        assert "NewECP Dy" in b_dy and "La" not in b_dy and "NewECP La" in b_la
        new_inp = norm(h_la) + b_dy + t_la
        assert new_inp == norm(dy_inp_old), "rebuilt Dy input != project's earlier Dy_ssaq_w21 input"
        assert t_la.startswith("* xyzfile 3 1 start.xyz")
        L = (la / "job.xyz").read_text().splitlines()
        n = int(L[0]); atoms = [l for l in L[2:] if l.strip()]
        assert n == 64 == len(atoms) and atoms[0].split()[0] == "La" and all(a.split()[0] != "La" for a in atoms[1:])
        new_atoms = [atoms[0].replace("La", "Dy", 1)] + atoms[1:]
        name = f"Dy_ssaq_w21_s{s}"; d = OUT / name
        d.mkdir(parents=True, exist_ok=False)
        (d / "job.inp").write_text(new_inp)
        (d / "start.xyz").write_text(f"{n}\nDy swap of {la.name}/job.xyz (La final Opt geometry), coordinates unchanged\n"
                                     + "\n".join(new_atoms) + "\n")
        diff = [l for l in difflib.unified_diff(la_inp.splitlines(), new_inp.splitlines(), lineterm="", n=0)
                if l[:1] in "+-" and not l.startswith(("+++", "---"))]
        meta[name] = {"la_source": str(la), "la_job_inp_sha256": sha(la / "job.inp"), "la_job_xyz_sha256": sha(la / "job.xyz"),
                      "la_job_out_sha256": sha(la / "job.out"), "dy_basis_source": str(dy_old / "job.inp"),
                      "dy_basis_source_sha256": sha(dy_old / "job.inp"), "new_job_inp_sha256": sha(d / "job.inp"),
                      "new_start_xyz_sha256": sha(d / "start.xyz"), "la_final_geometry": cn_report(atoms),
                      "diff_vs_la_lines_minus": sum(l.startswith("-") for l in diff),
                      "diff_vs_la_lines_plus": sum(l.startswith("+") for l in diff)}
        print(name, json.dumps(meta[name]["la_final_geometry"]), "diff -/+", meta[name]["diff_vs_la_lines_minus"], meta[name]["diff_vs_la_lines_plus"])
    (OUT.parent / "BUILD.json").write_text(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
