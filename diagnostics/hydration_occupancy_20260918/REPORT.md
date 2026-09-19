# Joint water arrangements: complete electronic table, occupancy model unfinished

**All12 arrangements /24 metal endpoints are complete.** Twenty new native DFT
energy/analytic-gradient calculations succeeded; four compatible full-water
endpoints were reused. All16 possible one-water additions shift the electronic
contrast toward La in this fixed-context experiment. This is a useful mechanism
result, not16 independent biological successes or an equilibrium occupancy model.

Production/default and released PQQ inputs, energies and bands remain unchanged.
The two alpha-lactalbumin structures are one consumed biological group. This
experiment establishes no additional classification accuracy beyond the preceding
water-orientation preparation result.

## What the calculations show

`R = E_Ca − E_La`, converted once with the existing627.509474kcal/mol/Hartree
constant. Positive changes are more La-like. No aquo offset or old decision band
is applied to this new preparation protocol.

| Structure | Variable water | Change in R on addition, kcal/mol |
|---|---|---:|
|1F6S|A211|+1.290 to+2.259|
|1F6S|A212|+5.104 to+6.074|
|6IP9|A310|+4.780 to+6.221|
|6IP9|A322|+4.756 to+6.089|
|6IP9|A326|+3.292 to+4.963|

Ranges span every arrangement of the other variable waters; they are not
statistical confidence intervals. Adding the entire variable inventory shifts
R by+7.363304kcal/mol in1F6S and+14.988150kcal/mol in6IP9 relative to retaining
none of those variable waters. One/two outer waters remain present throughout.

**Water identity matters even at the same water count.** With one variable water,
1F6S's lowest prepared DFT state usesA211 forCa butA212 forLa. The alternatives
are1.0028 and2.8114kcal/mol higher, respectively. In6IP9, the lowest two-water
prepared states also differ (Ca:A310+A326; La:A322+A326), although the competing
choices are only0.4709/0.2242kcal/mol higher. These are conditional electronic
rankings over prepared states, not free-energy populations or proven DFT minima.

**Single-water contributions cannot simply be added.** For6IP9A322+A326,
`E(pair)+E(empty)−E(A322)−E(A326)` is+7.499821kcal/mol forCa and+6.656096 forLa.
These values include all changes in the electronic model and state-specific
hydrogen reorientation; they do not isolate a unique physical cause.

The previously reported negative effects of poorly oriented waters cannot be
used as the occupancy model for these repaired preparations. The expanded
hydrogen-bond context also differs from the earlier small-core deletion trial;
the difference between those trials is not attributed solely to orientation.

## What remains unknown about occupancy

Water exchange uses `DeltaG_add = G(site+water)−G(site)−mu_liquid`.
The computed bulk reference is reused unchanged. The table reports exact
thresholds for the still-uncomputed bound-state contribution; no missing term
is filled with zero and no Boltzmann occupancy is supplied.

For a concrete example, addingA322 to6IP9'sA310+A326 arrangement has the following
**conditional diagnostic**, assuming bound internal water ZPE equals gas internal
ZPE and leaving the remaining correctionC explicitly unknown:

- Ca: `DeltaG_add = +3.175 + C_Ca` kcal/mol.
- La: `DeltaG_add = −2.300 + C_La` kcal/mol.

Thus a water can be costly forCa and favorable forLa before the remaining terms
are determined. The signs can change with realistic-sized missing contributions;
this is not permission to assume `C=0`, equal metal corrections, or an actual
occupancy. Bound internal vibration shifts, translational/rotational confinement,
multiple basins, and non-electrostatic solvent terms remain unresolved. The
protein/oxygen framework is fixed, so oxygen/scaffold relaxation is also absent.

## Cheap preparation and physical motion

The existing native unmasked MACE-OMOL model prepared rigid water orientations;
DFT supplies all energies compared across water counts.32 searches across16
nonempty, nonfull metal states converged (407 energy/force objective evaluations).
The empty states need no MACE calculation.

The two seeds reach different minima for6IP9Ca patterns001 and101 (MACE energy
gaps1.9906/2.1503kcal/mol). The declared lower-MACE-energy geometry was selected
within each fixed composition. Both alternatives and their receipts remain.
Their DFT ordering has not been tested; no global-minimum or complete-basin claim.

Saved forces/analytic gradients yield physical water translation/rotation
projections at32 retained-water endpoints without additional energy calls.
Translation-direction cosine has median0.959109, but minimum0.061715 forA322 in
6IP9La pattern110. DFT translation-gradient norms range3.407–30.650kcal/mol/A.
Strong average agreement does not validate all water-motion directions or a
curvature/entropy model. The next local check should include the problematic
partially occupied state as well as a full state; see[NEXT_MOTION.md](NEXT_MOTION.md).

## Missing-site proposals

A35-protein-heavy-atom fit (RMSD0.596439A) maps the two deposited water sites in
1F6S to two sites in6IP9. Transferring6IP9A322 produces an unmatched candidate
in1F6S,2.953A from the metal. It is also1.912A from existing waterA211 and2.235A
from a protein oxygen. It is **not ready to insert**: augmented site inventories
need water-water exclusion/geometry handling. No transferred water was added to
any quantum calculation. All proposed positions and actual target neighbors are
retained in[SITE_PROPOSALS.json](SITE_PROPOSALS.json).

## Execution and validation

- MACE1201908:291s, oneGPU/16CPUs,4656allocated CPU-core-seconds.
- DFT1201910:1795s on64CPUs, four16-rank workers,114880allocated core-seconds.
  Twenty successful endpoints, no retries. Endpoint wall times313.407–392.178s;
  median354.632s. The executor measured1792.226s separately from scheduler wall.
- Total new allocation:119536CPU-core-seconds and291GPU-seconds. Preparation,
  analysis and fixture-test CPU were not individually metered. Prior full-state
  and bulk-reference costs are reused historical work, not new computations.
- MACE worker peak host RSS2,077,280KiB; peak CUDA allocation1.876838GiB.
  Slurm batch MaxRSS reports1,162,716KiB(MACE) and12,031,020KiB(DFT), a separate
  accounting measure. Requested host RAM64474MiB/256GiB, respectively.
-49 real-fixture tests pass, zero skips:26water/preparation/collection/algebra/
  alignment checks and23 archived baseline/development tests. Old full-water
  MACE collection also reproduces exactly. Parser/algebra tests launch no quantum
  executable; the actual scientific integrations are the two completed jobs.
- Native DFT orientation comparators1201824/1201825 from the preceding experiment
  remain running. They were neither duplicated nor interrupted and their cost
  is not charged to this arrangement experiment.

These are development costs for complete enumeration, not the cost of a proposed
routine scanner score. No production occupancy scorer is enabled.

New protocol IDs:
`mace_omol_rigid_water_occupancy_proposals_v1` and
`native_r2scan3c_cpcm_mace_prepared_water_arrangements_v1`.
No new aquo reference, fitted threshold or default change.

## Artifacts and next action

- [Unrounded results](RESULT.json),[motion diagnostic](MOTION.json),
  [allocation receipts](ACCOUNTING.txt),[commands](COMMANDS.md).
- [Readable state table](export_v1/TABLE.md),[state CSV](export_v1/states.csv),
  [all addition edges](export_v1/water_additions.csv),[figure](export_v1/water_arrangements.svg).
- Workspaces: `workspaces/hydration_occupancy_20260918/{states_v1,mace_v1,dft_v1}`.

**Recommendation: pursue the occupancy model.** This experiment identifies
metal-dependent water-position rankings, substantial coupled effects, and
specific transitions sensitive to the missing free-energy terms. Next test
local water motion against DFT, including the weak MACE translation case;
retain the existing discriminator while this separate model is developed.
