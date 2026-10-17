# DFT supports the large local structural penalty

The four-state PBE0 reference supports the strong La penalty between the two
consumed LBT3 source structures seen with MACE. It does not support substituting
GFN2 simply because that model reduces source dependence. The next question is
which physically supported structural motions can accommodate this strain.
No preference classifier or production default has changed.

## Matched local comparison

All methods use the exact same four 48-atom cores, charge −1, from La-conditioned
7CCO and **Lu-conditioned** 7CCN, evaluated with La and Dy. 7CCN is not a Dy
structure. Define S(M)=E_M(7CCN)−E_M(7CCO), with same-metal differences formed
in native units before conversion. This is a local vacuum source response,
not a binding free energy or a measured structural population.

| Method | S(La) | S(Dy) | S(Dy)−S(La), kcal/mol |
|---|---:|---:|---:|
| MACE OMOL | +77.640037 | +3.950156 | −73.689881 |
| Standalone GFN2 | −14.563185 | −22.508206 | −7.945021 |
| PBE0-D4 | +69.911132 | +3.584798 | −66.326335 |

MACE minus PBE0 is +7.728904 for La, +0.365358 for Dy and −7.363546 for
the differential. GFN2 differs by −84.474317, −26.093003 and +58.381314.
Thus the broad MACE response has independent electronic-method support, with
a material quantitative La discrepancy. The prior 1 kcal/mol interpretation
scale is unchanged; it is not a statistical uncertainty or a fitted threshold.

This result follows the first branch of LOCAL_REFERENCE_INTERPRETATION.md.
A flatter GFN2 result is not a supported accuracy improvement here. Conversely,
agreement with an approximate DFT treatment does not prove biological accuracy.
DFT remains eligible for subsequent La/Dy work where justified; these results
do not require replacing MACE or conducting a functional sweep.

## Reference and numerical limitations

PBE0-D4/def2-TZVP/def2-J RIJCOSX uses DefGrid3, VeryTightSCF, native ORCA6.1.1
EnGrad and metal lcecp1-TZVP/AutoAux-J. La ECP46 and Dy ECP55 each leave 200
explicit electrons and an effective singlet. Physical La1/Dy6 metadata remain
separate: the Dy4f9 shell is frozen in the core, with no spin–orbit treatment.
MACE receives physical multiplicity1/6; GFN2 is an effective closed-valence model.

All four endpoints passed the original native convergence/state/gradient
admission. All four all_printed_residuals_pass flags remain false; detailed
native density/orbital residuals and small-gap messages remain in SCOUT_REPORT.md,
REMAINING_REPORT.md and the result. Native convergence is not electronic stability
or a numerical derivative test. Forces and physical Dy-state accuracy are not
qualified by this energy comparison.

The 48-atom region is a capped diagnostic, not the whole peptide. Existing
whole-peptide free-terminal chemistry remains a modeling assumption; primary
experimental terminal specifications were not recovered. See ../LBT_TERMINAL_CHEMISTRY.md.
No source was removed, no cap/charge/water inventory was changed and no reserved
SpyCI-LAMBS outcomes were opened.

## Execution and costs

Four molecular jobs12826/12844/12845/12846 each used12MPI slots, lasting
839/1068/866/995 seconds respectively: **45,216 allocated CPU-seconds** in total.
Each used its actual48,000MiB grant to derive MaxCore3100MB/rank. Scheduler
TotalCPU and MaxRSS are retained separately in the endpoint reports.
Three earlier pre-call serialization failures cost72CPU-seconds; the saved
JSON diagnosis12840 cost2CPU-seconds. The completed saved comparison12936 used
7seconds on1CPU, with zero molecular calls. No endpoint was repeated.

The complementary four-force MACE job12937 completed all four states in11seconds
on oneH200/12CPUs (132CPU-seconds and11GPU-seconds allocated). Every energy
exactly reproduced its saved12471 value; model parameters were unchanged.
These are four new force observables on existing coordinates. Its original
wrapper was first rejected for requesting a GPU in standard, then job12935
failed before model calls because the launcher selected an environment without
PyTorch. Both failures remain preserved. Launcher-only recovery used the actual
12500 GPU partition/virtualenv; no scientific code, state or gate changed.
Job12935 cost12CPU-seconds and1GPU-second allocated. The separate saved force comparison12938 subsequently completed all four cases
in10seconds on1CPU. Target vector differences are0.100–0.193eV/Å; all-atom RMS
differences0.154–0.182eV/Å, with no force-accuracy threshold or qualification.
See ../lbt_local_force_comparison/REPORT.md for the complete vectors and limits.

## Reproduction and decision

RESULT.json is the exact output of saved comparator12936, SHA256
95a988fd90e9af1bbc16f293b24b48bcc7dfdeea49ad349c520ac6e1bedc9f7d.
It contains all four owning reference pins, source identities, native differences,
continuous deviations and numerical flags. Exact source/config/commands are in
EXECUTION.json and the immutable lbt_local_reference_comparison_v2 workspace.
The three-result V2 adapter repairs only JSON tuple/list representation;
original scout and failed attempts stay preserved.

Next: use the completed same-coordinate force comparison to help inspect
whether the earlier torsion-only freedom left relevant bond or coordination
strain fixed. Declare any new mechanics test before execution. Do not automatically
widen bounds, repeat a completed search, fit expected preferences or promote
GFN2 on the basis of its smaller source spread. The source-robust LanM profiler
goal remains unmet.
