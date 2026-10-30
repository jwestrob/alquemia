"""CPCM metal-radius sensitivity: ORCA fell back to 2.40 A for both La and Dy.
CPCM single points at the conformer-best geometries with the metal radius set to
  r_ionic: La 2.40, Dy 2.14 (2.40 x Shannon CN9 ratio 1.083/1.216)
  r_small: both 2.00     r_large: both 2.80
Everything else is the production model (large-core PBE0-D4/CPCM)."""
import make_tests as M
RAD = {"r_ionic": {"La": 2.40, "Dy": 2.14}, "r_small": {"La": 2.0, "Dy": 2.0}, "r_large": {"La": 2.8, "Dy": 2.8}}
Z = {"La": 57, "Dy": 66}
runs = []
for t, rad in RAD.items():
    M.HEAD[t] = M.HEAD["lc_cpcm_ion"]
    for s in M.STATES + ["La_ion", "Dy_ion"]:
        m = s.split("_")[0]
        x = f"1\nbare {m}3+\n{m} 0.0 0.0 0.0\n" if s.endswith("_ion") else M.geom(s).read_text()
        q = 3 if ("_aq" in s or s.endswith("_ion")) else 2
        r = M.write(s, t, x, q, m)
        p = M.OUT / r / "job.inp"
        p.write_text(p.read_text().replace("* xyzfile", f"%cpcm\n radius[{Z[m]}] {rad[m]:.2f}\nend\n* xyzfile"))
        runs.append(r)
(M.OUT / "RUNS_radius.txt").write_text("\n".join(runs) + "\n")
print(len(runs))
