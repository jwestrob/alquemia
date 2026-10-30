# Bounded ALMO implementation tests

Status at preparation: **unexecuted**. Root owns all numerical Slurm execution
under `--qos=standard_expand`; this package performs no installation or submission.
Passing these tests qualifies sampled implementation behavior, not the La/Dy
quartet's scientific result or its archived Hamiltonian conversion.

## Fixed systems and method separation

`tests/FIXTURE.json` contains the six literal coordinate rows from CP2K 2026.2
`tests/QS/regtest-almo-strong/atomic-water.inp`. The first and second OHH groups
are two neutral singlet water fragments, each with ten electrons. Coordinates
are unchanged. The original regression is periodic GTH/BLYP with a different
fragment setup; **its published regression energy is not our numerical oracle**.

The new molecular benchmark explicitly uses a 16 Å nonperiodic cell, full
nuclear Coulomb interactions, no ECP, and literal H/O spherical 6-31G* entries
copied from this same CP2K release's `data/BASIS_SET`. Those exact decimal
exponents/contractions are also supplied to PySCF, bypassing named-library
rounding differences. Each water has 18 AOs; the complex has 36 AOs/20 electrons.
The accompanying CP2K `ALL` potentials are the full-nuclear GAPW representation,
not a GTH valence Hamiltonian. Source hashes and the literal files are retained.

Two distinct checks use this one fixed fixture:

- **Core tests:** restricted PBE0, full Coulomb J/K, no density fitting or
  dispersion; fixed unpruned PySCF 99-radial/590-angular grid. This exercises
  the hybrid exchange factor and numerical XC derivative used by the proposed
  quartet implementation.
- **Independent native comparison:** restricted PBE in PySCF and CP2K's own
  ALMO implementation, with no HF, ADMM, ECP or dispersion. PySCF uses an
  unpruned 150/974 grid. Two CP2K jobs use GAPW/MT, the same 16 Å cell, 400/600 Ry,
  relative cutoff 80 Ry, and 200/590 atomic integration grids. This tests an
  independent constrained electronic optimizer. It does not change the quartet
  PBE0/ECP protocol.

The installed ALMO regression inputs were inspected rather than blindly reused:
`regtest-almo-eda/almo-eda-x.inp` uses GTH/BLYP, periodic ice, XALMO_X and loose
fragment convergence. Its scalar reference cannot validate this molecular
constraint. A separately executed matched native comparison is therefore
required; there is no claimed literature R/F/P/U number for this exact fixture.

## Independent correctness checks

`tests/run_tests.py --mode core` builds the backend directly, without the archive
adapter. It runs two isolated fragments, constrained ALMO, ordinary full SCF,
one-fragment ALMO, and the public constraint-release path. There is no geometry
optimization, finite nuclear displacement, adaptive rerun or altered threshold.

1. Reconstruct density independently by Cholesky whitening of S, thin QR of
   the occupied vectors, and back transformation. This avoids the production
   occupied-overlap inverse expression. Compare density, Tr(DS)=20,
   DSD=2D, and exact forbidden coefficient support at frozen, nonstationary
   perturbed and optimized states. Density/idempotency error limit is 2e-9;
   forbidden coefficient limit is 1e-13.
2. Compare analytic dE/dC against **backend energies** reconstructed by that QR
   route, at four deterministic supported directions and two steps, 1e-4 and
   5e-5. Both fragment-only and combined directions are present. All eight
   slopes must pass `2e-7 + 1e-6*abs(analytic)` Hartree per coefficient step.
   Require nontrivial signal and rejection of a deliberately doubled gradient.
   This is 16 displaced density evaluations, not SCF or nuclear FD calculations.
3. Apply nonsingular, nonorthogonal mixing within each fragment's occupied
   columns. Energy and density must remain unchanged (2e-10 Hartree and 2e-9
   maximum density entry). This tests occupied-subspace gauge invariance.
4. Require rejection of an explicitly rank-deficient occupied matrix and a
   coefficient outside its fragment's AO support.
5. Require constrained convergence at the declared settings and
   `E_full <= E_ALMO <= E_frozen + allowance`, with 2e-8 Hartree allowance for
   each ordering comparison. Compare optimized energy against the independent
   backend-density route at 2e-10 Hartree.
6. Removing the partition to one fragment must recover ordinary SCF within
   2e-7 Hartree, starting from the genuine two-water frozen occupied space.
   Public constraint release must also converge to ordinary SCF within 2e-8.

Every optimizer is capped at 300 iterations/1,500 evaluations; ordinary SCF
is capped at 150 cycles with 1e-12 energy/1e-9 orbital-gradient targets.
ALMO targets are 1e-12 energy/1e-8 gradient. Failures remain failures. The
script records all prerequisite outcomes, nonstationary FD values, energies,
optimizer history, backend log, versions, source hashes and process costs.

## Native frozen-state definition and admission

In the pinned CP2K source, `src/almo_scf.F:716–793` normalizes/orthogonalizes
within fragments using `retain_locality=.TRUE.`, then calls
`almo_scf_t_to_proj` with `orthog_orbs=.FALSE.`, full S and occupied-overlap
inverse. It doubles the projector for the restricted density and evaluates
the actual KS energy. `src/almo_scf_methods.F:1460–1555,1582–1735` supplies
these routines. Thus native F is the antisymmetrized nonorthogonal occupied
determinant, not an unantisymmetrized sum of fragment densities. Within-block
normalization preserves the occupied space. The new solver uses the same
physical frozen-state definition.

Each CP2K run must have two explicit neutral-singlet fragment convergence
messages and final residuals <=1e-10, the expected atoms/18+18 AOs/5+5 occupied
orbitals, block-DIIS and FULL_SCF-PCG final residuals <=1e-8, no early stop or
restart, exit zero and a normal program end. Maximum native iterations remain
300, with no outer-loop extension. Missing or ambiguous records fail closed.
R and F come from their native records; P is the final block energy; U is the
final FORCE_EVAL energy independently joined to the full-SCF trace.

`compare_cp2k.py` checks the independently printed fragment-energy sum against
R and the full trace against U, plus native ordering at 1e-8 Hartree. It reports
FRZ=F−R, POL=P−F, CT=U−P and interaction=U−R. Telescoping closure is labelled
bookkeeping, never convergence evidence. Both native outputs and the PySCF PBE
benchmark must be complete; no partial passing benchmark is produced.

For each of FRZ/POL/CT/interaction, require both:

- Fine CP2K minus PySCF absolute difference <=1e-5 Hartree.
- Native 600-minus-400 Ry absolute difference <=1e-5 Hartree.

No tolerance is enlarged using the observed grid difference. That difference
is a measured resolution sensitivity, not a rigorous error bound. Raw R/F/P/U
and both fragment absolute-energy differences are also reported to expose
backend offsets and cancellation; exact absolute-energy equality across
different integration schemes is not assumed. Failure means that independent
agreement is not established, without uniquely blaming the optimizer or grid.

## Exact execution interface

Root snapshots all files, including extensionless BASIS/POTENTIAL assets, and
uses the installed PySCF 2.14.0 environment. Core and PBE benchmark jobs each
use one numerical thread; the scripts require actual Slurm/standard_expand.

```sh
<python> <snapshot>/tests/run_tests.py --almo-source <snapshot>/almo.py --mode core --output <fresh-core>/RESULT.json
<python> <snapshot>/tests/run_tests.py --almo-source <snapshot>/almo.py --mode benchmark --output <fresh-pbe>/RESULT.json
<python> <snapshot>/tests/run_cp2k.py --cutoff 400 --task <fresh-cp2k-400> --seconds <root-declared-cap>
<python> <snapshot>/tests/run_cp2k.py --cutoff 600 --task <fresh-cp2k-600> --seconds <root-declared-cap>
```

The two CP2K jobs may run concurrently. Each requires exactly 8 allocated CPUs
and invokes the **container's own MPICH** with one rank and eight OpenMP threads;
BLAS threads stay one. Build/source bind read-only and each fresh task binds
read-write. The runner verifies installed binary/image and deck/asset hashes,
records intent before its single parent process, keeps raw stdout/output, and
limits only its own process group. Root declares memory (32–64 GiB proposed)
and a finite native cap <=3,600 seconds before submission. No fake allocation,
automatic retry, job cancellation or queued-job modification is included.
Actual Slurm parent accounting is root's responsibility.

After **both** native attempts terminate, root supplies their actual exit codes:

```sh
<python> <snapshot>/tests/compare_cp2k.py --pyscf <fresh-pbe>/RESULT.json --coarse <fresh-cp2k-400>/endpoint.out --coarse-input <fresh-cp2k-400>/input.inp --coarse-exit-code <actual-code> --fine <fresh-cp2k-600>/endpoint.out --fine-input <fresh-cp2k-600>/input.inp --fine-exit-code <actual-code> --output <fresh-comparison>/RESULT.json
```

This is one fixed implementation panel. No quartet run, physical CT attribution,
charged-fragment qualification or new scientific campaign follows automatically.

## Method references

The original ALMO constraint is orbital AO-block support, with polarization
defined by relaxation within that constraint. See the
[original method paper](https://doi.org/10.1021/jp073685z) and the
[Q-Chem first-generation ALMO documentation](https://manual.q-chem.com/5.0/sect-EDA-CTA.html).
The executable independent authority is the installed, source-pinned CP2K
2026.2 ALMO implementation; the matched benchmark result remains pending.
