# Primary documentation and convergence distinction

The original native output identifies **ORCA 6.1.1**. Its complete output SHA is `6249c519168d9fdbda327a9d28eb21bfcf8cf587d56d6cf72792f80ee6315114`; exact input and recipe pins are in ASSESSMENT.json.

The [official ORCA 6.1 manual, SCF §2.6.7.4](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html#input-parameters), accessed 21 October 2026, documents the nested TRAH control:

> QuadRegionStart 1.e-4 — start Newton-Raphson if || G || < QuadRegionStart

The native 6.1.1 output independently prints that same default NR-entry threshold, 1.000e−4 (line 628), separately from its TRAH convergence threshold and TolG, both 2.000e−6 (lines 625 and 668; see actual file for the exact neighboring labels). The proposal changes only the former to 1e−6. It keeps VeryTightSCF and all printed convergence tolerances unchanged: TolG/TolX 2e−6, TolE 1e−9 Eh, one-electron-energy 1e−6 Eh and TolErr 1e−8. Active solver criteria must still be read from the actual recovery output; stale inactive-mode residuals cannot substitute for them.

No LooserSCF, changed TolG, increased MaxIter, smearing or modified electronic state is proposed. The manual separately distinguishes an orbital-rotation local minimum from the global electronic minimum; even successful recovery would not certify global state, spin–orbit physics or reference accuracy. The proposal is supported by the documented control and the observed oscillation, not by a guarantee that changing it will converge.
