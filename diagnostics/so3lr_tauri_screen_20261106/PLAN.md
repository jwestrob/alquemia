# Bounded Tauri screen — frozen 6 November 2026 before Tauri evaluation

Question: does the pinned Tauri L potential improve the local water/carboxylate
balance enough to earn one solvent feasibility study? Primary L-2606; M-2606 is
fixed capacity sensitivity, never a replacement chosen for agreement or speed.
HEAD f1181b4, unrelated dirty files and dnv1_fold45665 preserved. No prior Tauri
execution found. Scientific numerics use Slurm standard_expand, root submissions,
shared allocations, no scheduler memory/time overrides. New work only here and in
workspaces/so3lr_tauri_screen_20261106. Existing environments remain unchanged.

The exact 24 configurations: four complete W8/W9 La/Dy JSON geometries, two dry
acetate centers plus the sixteen archived radial/angular endpoints, two hydrated
AcW7 committed XYZs. Preserve atom order, decimal coordinates, source hashes and
operators. No relaxation, training or new biological/holdout data. Multiplicity
La1/Dy6, Q=3 aquo, Q=2 acetate. Omitted explicit spin defaults are inadmissible.

Source70707f5bb0c333ecec707e5a12b1efdf1c83603a; GLP2294f4f59e940be64c45bf5e4d1bea5cb90809f8;
other dependencies frozen in environment-lock.txt before model loading. New Python
3.12 environment. Float64 model/state, jax_enable_x64 before loading, highest
matmul precision, head5, isolated no PBC; short cutoff from model=5 A; primary LR
1000 A, dispersion damping2 A. Store raw device arrays in NPZ before any formatted
report. Native standard loader/graph/batch/make_so3lr for all frames; ASE as an
intact-state implementation check. Do not fake an applied field or add -mu.E.

Primary observables: radial metal–donor force components; collective water-shell
force (rigid OHH translations away from metal); dry-acetate radial and angular
energy derivatives in the original frames, full forces secondary. DeltaF=Dy-La
only for exact matching coordinates, atom order, composition and Q apart from
metal, with physical multiplicities separately explicit. Balanced total-energy
contrasts CW=(DyW9-DyW8)-(LaW9-LaW8), CAc=(DyAcW7-DyW8)-(LaAcW7-LaW8), in kcal/mol.
These are electronic configuration contrasts, not solution free energies.

Implementation gates on real fixtures: all24 Q/spin/head graph fields asserted;
LaW8 native/ASE energy and force differences <=1e-6 eV and eV/A; energy/force
finite difference along LaW8 water-shell direction, step1e-4 A, <=1e-4 eV/A.
One translated REAL dry-acetate frame with metal(+3)/acetate(-1) fragment metadata,
asserted one-based selections and native residue fields; compare constrained and
unconstrained only as a diagnostic. Never constrain bound periodic fragments.
CW at LR100,1000,and no-cutoff where supported, independently for both models.
One real water spectator translated50 and100 A from LaW8, Q3,m1; report force and
nonadditive-energy sensitivity, without a fitted correction. No imaginary atoms
or manufactured primary perturbations. Warm timing on LaAcW7, ten evaluations;
compilation separate; small-cluster timing is not bulk-solvent MD throughput.

Before comparison audit old gas/sc_gas inputs, outputs and exact geometry for
DyW8/W9/ion and La partners: native ECP/electrons, occupations, S2, residuals,
SCF/root/stability evidence, energy closure. Reproduce the reported 27.55 kcal/mol
core discrepancy from raw outputs, but do not treat normal termination as state
qualification. Their optimized ingredient geometries/methods differ from this
panel and no matched physical-Dy gradients are supplied.

Conditional reference release: if audit confirms no compatible reference, exactly
six fresh fixed energy/gradient calls: La/Dy W8,W9,AcW7 at these archived coordinates.
wB97M-V/def2-TZVPD, full pinned external basis from FAIRchem d4e10d0f4e87593922fda1d1ffb7dd28e504cadd;
RIJCOSX def2/J NoUseSym DIIS NOSOSCF NormalConv DEFGRID3 ALLPOP NONBO NONPA;
Tight SCF/maxiter300, THRESH1e-12/TCUT1e-13; EnGrad; no appended D4/SCNL. ORCA6.1.1
port (original6.0). PModel guess declared before runs for all six because the
archived La PAtom startup is known unsupported; this changes initialization, not
the Hamiltonian. UKS for all, original20-degree beta HOMO/LUMO rotation for La
singlets only. Request STABPerform true, STABRestartUHFifUnstable false: inspect
orbital stability without automatic root search. No MO reuse, optimizer, retries
or additional SPs. Twelve MPI/job; allocation-derived MaxCore with25% headroom.
Capture actual automatic NL/EnGrad behavior, ECP and basis rather than assuming
from a dictionary. Actual La46/Dy28 core/electron counts, finite gradients and SCF
criteria/occupations/S2/stability are mandatory reference records. Diagnostic
qualification requires all printed Tight residuals passing, spin deviations <=.5
from S(S+1), no negative stability eigenvalues below -1e-4 Eh; if unavailable,
failed or ambiguous label unresolved. Those are provisional screening gates, not
proof of a unique physical root or an exact identity with mixed Tauri training.

Screen rule frozen: for qualified references, >1 kcal/mol error on either balanced
contrast ends the route before MD. A wrong local radial-force sign with reference
magnitude >1 kcal/mol/A also stops. Old effective-singlet Dy data are descriptive
only; missing/ambiguous roots mean unresolved, never model failure/success.
Numerical interface errors must be resolved separately; preserve failed attempts.
No acceptance thresholds or model/feature changes after inspecting values.

Experimental admission is independent: test ONLY599/600 at1M (including trans's
alternate row), then fallback576/597 at.1M. Verify raw native annotations, species,
protonation, temperature, ionic strength, explicit electrolyte and source. K+ or
Na+ salt without the anion does not establish shared medium. If source evidence
cannot resolve it, neither pair is admitted and no substitute pair is screened.
The conditional D values are development references, not pristine holdouts.

Only a passed physical screen AND a clean experimental pair can release a solvent
feasibility pilot. No automatic MD here. Before any such pilot, freeze bound
association domains, periodic electrostatics and sensitivity, exact electrolyte,
charge/spin, independent starts and one-H200-day aggregate budget. The first
solvent report would concern stability, mixing/overlap and measured cost, never
endpoint means or a trapped BAR estimate as the requested D. No peptide or full
protein work is released. M1 baseline and reserved SpyCI-LAMBS/PLM outcomes remain
untouched. End with reject, unresolved, or one justified next calculation/cost.
