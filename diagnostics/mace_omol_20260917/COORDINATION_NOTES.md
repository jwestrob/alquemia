# Interpretation notes for the matched coordination candidate

Primary-source check during the running frozen benchmark, 2026-09-17. This
does not change inputs, method, thresholds or execution scope.

The OMol25 dataset paper reports 2–350 atoms (mean50), charges -10 to+10,
and multiplicities1–11. It explicitly identifies charge/spin localization as
a limitation of a global embedding. Its distance-scaling evaluation distinguishes
connected and disconnected neighbor graphs at6Å; metal-complex scans extend a
single metal–ligand separation by up to4Å. These evaluations do not establish
accuracy for our fully separated metal references. The paper's original MACE
baseline used neutral training only and is **not** the later charge-conditioned
MACE-OMOL-0 checkpoint used here. See [dataset profile, model discussion and
Appendix H.6](https://arxiv.org/html/2505.08762v1).

The later checkpoint's own [training configuration](https://raw.githubusercontent.com/ACEsuit/mace-foundations/main/mace_omol/mace-omol.sh)
and our loaded-model/actual-batch receipts establish the charge/spin inputs for
this experiment. Neither publication-level average errors nor finite-distance
scans are validation of these specific La/Ca protein fragments.

Our cancellation of geometry-independent terms is exact model algebra. Our
repeat/rotation/distance checks establish numerical behavior. Neither certifies
the electronic state of each separated fragment. A successful direction test
would therefore support a useful empirical coordination descriptor, subject to
independent testing, rather than a computed solution binding free energy. A
failure must not be rescued by moving the reference closer or changing charges
after inspecting its predictions. Any such alternative is a separate model.
