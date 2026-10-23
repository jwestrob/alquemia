# Common-basis water-shell response

## Question and scope

Does unequal electronic basis flexibility between isolated fragments and full
water shells distort the inferred La/Dy collective-response error? Execute the
finite design in `explicit_model_v1/paired_water_coordination_joint_v1/mechanics_next_v1/AFTER_PAIR_RESULT.md`.
The goal was reactivated by the user. Existing authorization covers this test
and subsequent evidence-led mechanics development.

Use only the four consumed W8/W9 La/Dy configurations and their literal atomic
coordinates. Reuse four admitted full complexes and eight admitted metal/shell
counterpoise fragments. Calculate 34 metal-plus-water fragments and 34 single
waters, each with every parent basis center retained. Inactive centers are
ghosts; their basis functions remain, but their nuclear charges, electrons and
ECP potentials do not. A real metal retains its original orbital/auxiliary
basis and 46-electron La or 55-electron Dy core. No optimization, force call,
parameter fit, new model evaluation or reserved validation outcome is included.

## Frozen calculation and arithmetic

Use the qualified CP8 PBE0/RIJCOSX/DefGrid3 electronic recipe without D4, exact
source orbital and auxiliary bases, original SCF settings and effective singlet
state. Metal-plus-water has charge +3 and 18 explicit electrons; water-only has
charge zero and 10 electrons. W8/W9 retain 25/28 centers and 391/434 orbital
basis functions. Record actual native atom/basis/electron/ECP tables and SCF
stopping criteria. Energy-check admission does not imply that every printed
density/DIIS criterion, stability or physical electronic state is qualified.

For each parent, C = E_full - E_ion - E_shell,
P = sum_i(E_metal_water_i - E_ion - E_water_i), and B = C - P.
The ion in P is subtracted N times. Form the paired difference
Dy9 - La9 - Dy8 + La8 only after computing the complete per-parent quantities.
Ghost-basis water, shell and ion energies cannot be assumed to cancel across
different parents. Require all 68 new and 12 reused cells; missing terms make
the complete decomposition unavailable. C must replay the completed CP8
electronic result. Preserve original D4 C/P/B contributions separately and add
each once when displaying the electronic-plus-original-D4 sensitivity.

## Execution and interpretation

Root owns submissions. First run the W9 Dy water-01 metal-plus-water and
water-only cells, one invocation each. Independently inspect their native
contracts, then submit the remaining 66 unchanged cells if they pass. These
two scouts belong to the total of 68 and must not be repeated. No automatic
retry or replacement state is authorized by this finite execution record.

Each cell uses four single-thread MPI ranks, `banfield/standard` with
`standard_expand`, scheduler-default memory/time and allocation-derived
MaxCore. Live policy on resumption permits 24 CPU/job, 500 CPU/user, 20 running
and 200 submitted jobs. Use independent parallel jobs within those limits;
retain exact commands, source pins, completion wakes and actual costs.
Preparation/review uses one local thread per agent/root, within the aggregate
16-thread login allowance. Saved-result collection is a separate one-CPU task.

Compare B's basis sensitivity with the existing 5.602 kcal/mol discrepancy and
the declared 1 kcal/mol development energy resolution. Repartitioning cannot
erase the already measured full-shell discrepancy. A persistent B error
motivates a physically distinct collective-response test; it does not uniquely
identify polarization, charge transfer or penetration. Both references and
sampling models still require qualification before exchange production.
