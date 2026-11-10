"""Fingerprint of a single rare-earth pattern (PLAN.md): additive step between Gd and Tb on a cubic polynomial in
Shannon nine-coordinate radius, the step in element units, and the yttrium offset. Lead agent, Claude Fable."""
import numpy as np
LN = ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
R9 = dict(zip(LN, [1.216, 1.196, 1.179, 1.163, 1.132, 1.120, 1.107, 1.095, 1.083, 1.072, 1.062, 1.052, 1.042, 1.032]))
R9_Y = 1.075
HEAVY = {"Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"}
VARIANTS = {"no_Ce": {"Ce"}, "no_Ce_La": {"Ce", "La"}, "no_Ce_La_Eu": {"Ce", "La", "Eu"}}


def fingerprint(vals, y_val=None, exclude=("Ce",), min_points=9):
    """vals: dict element -> log10 value (higher = more bound or more concentrated). Returns dict or None."""
    el = [e for e in LN if e not in exclude and e in vals and vals[e] is not None and np.isfinite(vals[e])]
    if len(el) < min_points or not ({"Gd", "Tb"} <= set(el)):
        return None
    r = np.array([R9[e] for e in el]); v = np.array([vals[e] for e in el], float)
    z = (r - 1.124) / 0.06
    X = np.c_[np.ones(len(el)), z, z ** 2, z ** 3, [1.0 if e in HEAVY else 0.0 for e in el]]
    b, *_ = np.linalg.lstsq(X, v, rcond=None)
    res = v - X @ b; dof = len(el) - 5
    s2 = float(res @ res / dof) if dof > 0 else np.nan
    cov = s2 * np.linalg.inv(X.T @ X) if dof > 0 else np.full((5, 5), np.nan)
    h = float(b[4]); h_se = float(np.sqrt(cov[4, 4])) if dof > 0 else np.nan
    zm = (0.5 * (R9["Gd"] + R9["Tb"]) - 1.124) / 0.06
    dfdr = (b[1] + 2 * b[2] * zm + 3 * b[3] * zm ** 2) / 0.06
    s = float(dfdr * (R9["Tb"] - R9["Gd"]))                  # change of the smooth curve over one Gd->Tb spacing
    out = dict(n=len(el), h=h, h_se=h_se, slope_step=s, u=(h / s if abs(s) >= 0.02 else np.nan), rms=float(np.sqrt(s2)) if dof > 0 else np.nan)
    if y_val is not None and np.isfinite(y_val):
        zy = (R9_Y - 1.124) / 0.06
        out["y"] = float(y_val - (b[0] + b[1] * zy + b[2] * zy ** 2 + b[3] * zy ** 3 + b[4]))
        for tw in ("Dy", "Ho", "Er", "Tm"):
            if tw in vals and vals[tw] is not None and np.isfinite(vals[tw]): out[f"Y_minus_{tw}"] = float(y_val - vals[tw])
    return out
