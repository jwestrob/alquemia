"""Water H-bonds and metal-donor distances in the selected EDTA/CDTA/TMDTA states (final job.xyz)."""
import sys
from pathlib import Path
import numpy as np
B = Path("/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces")
OLD, NEW = B / "lady_relaxed_exchange_20261024/runs", B / "lady_panel_20261028/runs"
def load(n):
    p = (OLD if (OLD / n).exists() else NEW) / n / "job.xyz"
    L = p.read_text().split("\n")[2:]
    a = [l.split() for l in L if l.strip()]
    return [x[0] for x in a], np.array([[float(v) for v in x[1:4]] for x in a])
def report(n):
    el, X = load(n)
    D = np.linalg.norm(X[:, None] - X[None], axis=-1)
    m = 0
    bondH = lambda i: [j for j in range(len(el)) if el[j] == "H" and D[i, j] < 1.2]
    bondC = lambda i: [j for j in range(len(el)) if el[j] == "C" and D[i, j] < 1.6]
    wat = [i for i in range(len(el)) if el[i] == "O" and len(bondH(i)) == 2 and not bondC(i)]
    ligO = [i for i in range(len(el)) if el[i] == "O" and bondC(i)]
    don = sorted([(D[m, i], el[i] + ("w" if i in wat else "")) for i in range(1, len(el)) if el[i] in "NO" and D[m, i] < 3.1])
    hb = []
    for w in wat:
        if D[m, w] > 3.1: continue
        for h in bondH(w):
            acc = [(D[h, j], "Olig" if j in ligO else "Ow") for j in ligO + wat if j != w and D[h, j] < 2.3]
            hb.append(min(acc)[1] + f"{min(acc)[0]:.2f}" if acc else "free")
    s = lambda t: np.mean([d for d, e in don if e == t]) if any(e == t for _, e in don) else float("nan")
    print(f"{n:13s} CN {len(don):2d}  <M-N> {s('N'):.3f}  <M-Olig> {s('O'):.3f}  <M-Ow> {s('Ow'):.3f}   water H-bonds: {' '.join(hb)}")
for n in sys.argv[1:]: report(n)
