# Frozen-charge transfer challenger v1

**Baseline remains default.** This opt-in model is a perturbative electrostatic
descriptor. Geometry, direct donors, waters, PQQ(3-), native r2SCAN-3c/ECP,
CPCM(water) and singlet charge policy are retained within every pair.
No threshold is inherited by this challenger.

## Energy accounting

Let q be one endpoint's MBIS atomic charges, Q the fixed protein charges after
the documented covalent-boundary exclusion, T the physical protein cavity, R
the isolated capped-core cavity, and H a homogeneous dielectric of 1. APBS
returns a charging energy F including grid self energies. Define RF_D(x) =
F_D(x)-F_H(x), always on identical grids within each subtraction.

    DeltaU_M = Coulomb(q_M,Q)/epsilon_in
               + RF_T(q_M+Q) - RF_T(Q) - RF_R(q_M)
    R = E_Ca - E_La
    S = R - (E_Ca,aquo - E_La,aquo)
    S_env = S + DeltaU_Ca - DeltaU_La

The explicit Coulomb term includes all core/environment pairs, including cap
charges; it is not just the potential at the metal. Environment-only terms
cancel in the common T cavity. Homogeneous intramolecular core terms cancel;
grid self terms cancel within each RF difference. The homogeneous grid cross
term is retained as an independent numerical check against analytic Coulomb.
Reaction-field core self and core/environment cross effects remain in DeltaU.
No new nonpolar cavity energy, dispersion, mechanical relaxation or entropy is
added. Polarization represented by fixed ff19SB charges and solvent dielectric
is an effective approximation, not self-consistent protein electronic response.

The original energy already contains quantum CPCM. The subtracted RF_R is a
matched *point-charge APBS* reference, not the exact CPCM functional. Thus the
transfer approximates changing the environment of the CPCM-derived density;
it does not claim exact removal/replacement of quantum CPCM solvation.

## Frozen state and numerical settings

- Development structures: already consumed 1H4I fixed core, qm33 and qm36.
  The qm33/qm36 pair alone is the partition test; fixed-core protonation is
  not silently equated to the older boundary source.
- Physical assembly: full deposited catalytic chain A for each input. This
  is explicitly a chain-level test, not the biological alpha2beta2 assembly.
- MBIS monopoles from each metal's own CPCM calculation. Require reported
  charge sum within 1e-4 e; never infer a 46-electron ECP charge shift or
  normalize charges. Check the actual quantum ESP on deterministic points
  outside atomic exclusion spheres before using charges in APBS. Quality
  gate: relative RMS <=10% or absolute RMS <=0.005 atomic units. This is a
  development approximation check, not a validated universal tolerance.
- Protein ff19SB protonation/coordinates stay fixed. Reuse the archived
  local boundary map where applicable: remove represented source atoms and
  MM1 CA charge; distribute the local closure correction equally to the
  bonded backbone N/C. No global neutralization. The partition test must
  judge this approximate representation; charge closure alone does not.
- Bondi H/C/N/O/S radii 1.20/1.70/1.55/1.52/1.80 A; fixed 1.80 A radius for
  both Ca and La. This common metal cavity is an explicit approximation.
  Full source protein atoms remain in T with zero charge if represented in QM.
  Synthetic cap spheres must be wholly contained in the physical protein
  cavity. PQQ and metal are QM-only, never assigned generic force-field charges.
- APBS 3.4.1 official Linux release, isolated under the task workspace.
  Linear Poisson equation, zero salt, dielectric 1/78.54, 298.15 K, molecular
  surface, 1.4 A probe, spl2 charges, mdh boundary. Primary grid spacing 0.5 A;
  refinement 0.4 A at identical physical box extent; nominal padding >=20 A
  and prescribed box extension >=30 A. Grid sizes obey APBS multigrid rules.
- Identity changes BOTH cavity and charges to R/q (Q absent), not just Q=0.
  Require |DeltaU| <=0.01 kcal/mol. Repeat, rigid transform, refinement and
  box-extension movements must be <=0.5 kcal/mol. Require partition-score
  movement <=2 kcal/mol. Failures disqualify this model; do not retune.
- Any unsupported cofactor, missing nearby atom, charge overlap, inconsistent
  coordinate mapping or failed solver creates an unavailable score.

## Cost and execution

Historical job 1198050 used 154 s on 344 allocated CPUs, 52,976 allocated
core-seconds for four endpoints. Individual endpoints used 16 MPI ranks and
117.10–152.80 s, concurrently; summing their wall times does not give job time.
Its batch MaxRSS was 12,401,728 KiB. These are archive measurements, not a new
performance result. Standard/memory hardware may differ and will be recorded.

Initial pilot: six new MBIS endpoints, at most two retries, no gradients or
displacement calculations. High-level admission budget 158,928 allocated
core-seconds; low-level/solver budget 52,976; aggregate 211,904. These one-time
development ceilings derive from 3+1 multiples of the measured campaign.
No ordinary-score cost claim follows. Production feasibility separately
requires <=2 times a matched baseline endpoint-pair cost and Jacob's promotion.
Finite task admission stops after budget exhaustion; local policy prohibits
subprocess timeouts, so a launched task can overrun its estimate. Overruns
are recorded and block further admission rather than being hidden.

## Primary implementation references

- [APBS matched self-energy subtraction](https://apbs.readthedocs.io/en/stable/using/examples/solvation-energies.html)
- [APBS Coulomb/reaction-field components](https://apbs.readthedocs.io/en/stable/using/examples/binding-energies.html)
- [ORCA MBIS](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/population.html)
- [ORCA multiscale coupling](https://www.faccts.de/docs/orca/6.1/manual/contents/multiscalesimulations/qmmm-molecules.html)

The installed ORCA 6.1.1 and actual outputs determine compatibility. The manual
documents CPCM/B as passing large-system surface charges to small-system
calculations; that does not establish fully self-consistent high-level protein
polarization. No integrated hybrid backend is substituted in this pilot.
