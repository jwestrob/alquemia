# Finite local Mex GFN2 field scout

Declared 16 October 2026 before new molecular calls. Root owns execution.

## Question and fixed source

Can the qualified standalone GFN2 implementation return converged local energies,
populations and analytic QM/external-site gradients for the actual repaired Mex8FNS
EF3 origin, with and without its existing shifted MM field? This is four calls:
La/Dy × isolated/embedded. It is not a whole-protein SCC calculation, a MACE
calculation, solvent calculation, optimization, affinity score or classifier.

Use configuration A of
`workspaces/metal_environment_response_20260926/lanm_ef3_CboundH_repaired_v2/Mex8FNS/INPUTS.json`.
Preserve its exact 209-atom XYZs: 205 real atoms including the target, four caps,
local charge0. Preserve `environment.pc` and `environment_atoms.json` order,
charges and positions: 1840 external physical sites, field charge+2, including
three formal Nd3+ spectators. Four zero-field boundary CA atoms remain part of
the 2049-atom physical system but are absent from the electronic field file.
The actual source has no virtual charge sites. No protonation, water, geometry,
cap, charge redistribution, spectator occupancy or source selection changes.

Both evaluations use original published spin-independent f-in-core GFN2,
effective UHF0. Physical local LaIII singlet and DyIII sextet metadata remain
separate; this does not qualify Dy open-shell response or whole-protein MACE
spin coupling with Nd spectators.

## External field and energy reference

Use the exact standalone xTB6.7.1 executable, author parameter file and runtime
already admitted by Dy qualification12488 and whole-LBT solvent12510. Pin the
installed API header/library and archived native source as interface evidence;
execution uses the established CLI and process helper, not a new engine/binding.

Embedded control is `$embedding` with `interface=orca`, `input=external.pc` and
`gradient=pcgrad`. The first four point-charge tokens are the original q[e] and
x/y/z[angstrom] tokens. Append the actual source element's `GAM` read from the
pinned author parameter file as an explicit numeric fifth token. In this parser
numeric means hardness, not atomic number. This gives H0.405771, C0.538015,
N0.461493, O0.451896 and Nd0.658531 in native atomic units. Do not use the silent
generic-N default. Preserve each site's element, ID, shifted charge, position
and selected hardness in metadata. `interface=orca` converts positions to Bohr.

This is native GFN2 hardness-damped shell electrostatics, not ORCA bare1/r
embedding. `electro_pcem` adds the QM–field interaction once; it adds no PC–PC
self-energy. The external potential uses its finite-space zero at infinite
separation, without offsets/background, periodic electrostatics or fitted gauge.
Embedded minus isolated includes local electronic response and field interaction.
It cannot be interpreted as binding energy.

A later Hamiltonian could be
`MACE_vac(Q) + GFN2_emb(Q;MM) - GFN2_vac(Q) + existing_MM_ledger`,
plus one separately defined consistent solvent term. Its field interaction is
already present in the difference; adding another QM–MM Coulomb term would
double-count it. This scout evaluates only the two GFN2 pieces. Existing MM
mechanics, MACE and solvent are not evaluated or assembled into a score.

## Settings, validation and outputs

All preparation and scientific evaluation run under shared Slurm. The separate
`prepare` stage performs zero molecular calls and validates all four states,
source pins, La/Dy paired coordinates, full physical/cap/field mappings, actual
hardness/field serialization and charge accounting before writing its manifest.
Root reviews that result before `run`. Inputs and preparation remain immutable.

Use exactly `--gfn 2 --chrg 0 --uhf 0 --etemp 300 --acc 0.0001
--iterations 500 --norestart --json --pop --grad`; no solvent. Four independent
single-thread workers at most. No retries, alternative SCF settings, seeded
restart, native-ORCA fallback or screening by energy sign.

Use existing local standalone gates: charge closure≤1e−5e; JSON/stdout energy
agreement≤1e−7Eh; engrad/stdout energy≤1e−8Eh; JSON/charge-file populations≤1e−7e;
engrad/XYZ positions≤1e−6Å. Require normal process termination, SCC convergence
evidence, correct printed method/version/charge/electrons/UHF/temperature/accuracy/
iteration limit/restart state, and correct PC on/off state. Embedded additionally
requires printed1840 point charges and sumPC+2. No silent solver failure or
unsupported numerical gradient. All arrays must be finite and correctly ordered.

Save209 QM gradients and1840 external gradients as dE/dR in Eh/Bohr and
kcal/mol/Å. Isolated external gradients are explicit exact zeros because that
Hamiltonian has no field, not fabricated backend output. `pcgrad` is the native
headerless N×3 output with eight decimals; do not use ORCA's count-prefixed parser.
Record its5e−9Eh/Bohr per-component rounding bound. Save total translation
residuals descriptively, without a tighter gate than serialization warrants.
Map cap derivatives through both existing anchor Jacobians and external-site
gradients by physical source IDs. Report embedded−isolated energies, populations
and mapped gradients only when both cells pass. Retain every failed cell in the
four-cell denominator. These returned analytic gradients are not yet independently
finite-difference qualified for this external-field protocol.

Pin raw outputs, implementation/config/plan, source mapping, preparation, parent
qualification and runtime. Record attempted/accepted counts, wall time, observed
CPU usage/RSS and application-allocated CPU-seconds; root records scheduler cost
and completion wake separately. No production or shared-default changes.
