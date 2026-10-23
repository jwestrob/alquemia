# Five-state cutoff12 term observer

Source preparation only. Authority: ROOT_DECLARATION.json (7b8af01e). No molecular execution, compilation, new geometry, initialization, factorization or solver is released by this package.

## Question and evidence

All nine inactive cutoff12 gradients were admitted. Whole-box Cartesian directional checks passed; all three strain1 checks failed. Strain residuals were 0.82748646 and 0.81429096 kcal/mol/Å. The two-step change was 0.013195495, of which 0.013192967 was nonpolar and 0.000002528 was polarization. That decomposition attributes the step change, not the absolute residual: saved outputs have total and polarization energies but only the total virial.

Native empole1/emrecip1 includes explicit reciprocal-volume/lattice virial, multipole-coordinate terms and moving-frame torque contributions. Native ehal1 includes the taper derivative. Polarization eprecip1 subtracts stored permanent reciprocal virial before its own accumulation. These facts preclude attributing the failure to a missing cell derivative from source inspection alone. All routines and their order must stay intact.

## Exact observer change

One override object contains the original gradient routine plus marked insertion-only calls and one new private-named helper. Removing every marked insertion recovers the pinned original gradient.f byte for byte. Every original kernel, accumulator and operation is retained. The helper writes a separate term_trace.dat and never assigns any imported native quantity. It records the initial zero virial, each of 28 ordered terms (including disabled terms), and the final virial after watfour2. Each term records its scalar energy, enabled flag and full cumulative virial. A final list of all 28 scalar energies must match their earlier snapshots. The optional VDW correction flag is also reported and required false for these exact baseline settings.

One observer object is linked before the unchanged libtinker.a, with the three existing main/gate/state-writer objects. Only the observer object is compiled. No library rebuild, kernel edit or response activation. GNU archive symbol inventory verifies gradient.f.o supplies only gradient_ as a global defined symbol; the explicit override resolves it without pulling that archive member. Actual link and execution remain separately released root actions.

## Required observation checks

For all five states, the existing main writes the original stage2_gradient.dat. Total energy, all 35,874 gradient components and all nine virial component tokens must match that state's admitted cutoff12 result exactly. Polarization energy and unit-conversion constants also match exactly. Ticket hashes legitimately differ; role/state must match. Any difference stops admission. No retry is automatic.

Term energy summation uses the existing arithmetic roundoff allowance only; virial increments are decimal differences of printed cumulative snapshots and telescope exactly. The final post-term virial increment must be zero. Every term stays in the report, including zero/disabled terms. The permanent/polarization split means their net increments in the original shared call order, not independently re-evaluated kernels.

## Strain units and analysis

Frozen DIRECTIONS specifies q in Å and L0=1 Å. Coordinates and cell transform with dimensionless epsilon=q/L0. Therefore dE/dq=(S:virial)/L0 is kcal/mol/Å, as the inherited helper labels it. The equivalent derivative dE/depsilon=S:virial is kcal/mol; the numerical values coincide because L0=1 Å. Nothing is rescaled or retuned. Use the original strain1 yy matrix and q steps ±1e−4 and ±5e−5 Å. The three original total gates retain 0.001 kcal/mol/Å plus literal-print allowances. Per-term slopes/errors are diagnostic decompositions, not new physical acceptance tests.

## Execution and limits

Five exact existing cutoff12 states: central and strain1 ±h/±h/2. One CPU per job, standard partition / standard_expand QoS / banfield account, default scheduler memory/time, no requeue. Existing child limits: compile/link each180s and3GiB; each gradient3600s and3GiB, all numerical thread variables1. New build/result/release pins remain null until root binds them. Native counts:5 gradients maximum;0 energy-only/initializer/factor programs. Build cost and five native allocations are reported separately. No complete force-field, active-response, periodic stability or trajectory qualification follows from this observer.
