# Native ALPB shielding: a bounded component diagnostic is feasible

**Recommendation:** test the electrostatic reaction-field change from adding the
matched protein's omitted, uncharged cavity centers to the existing capped core.
The native API supports this without SCC, new radii or new electronic energies.
Call it a **capped-union exclusion sensitivity**, not a source-only physical
cavity or a classifier correction. The synthetic caps remain extra cavity
centers; a finite cap-charge and cap-volume control is needed to interpret it.

This is a source/code feasibility audit on 14 October 2026. No scientific
calculation, numerical test, installation, submission or correction was performed.
The preceding midpoint experiment is complete and closed; this is a separate
proposed phase. [Midpoint result](../../exaa_metal_midpoint_20261014/REPORT.md).

## 1. Exact native energy and the proposed change

Inspected xTB 6.7.1 commit `26b28010e805f7d1aeeef39813feb473e69cc4be`, using the
actual archived source under
`workspaces/pede_solvent_response_20261013/fixed_charge_replay_v3/execution_v1/sources/xtb/`.
Relevant implementation: `src/solv/gbsa.f90` (`initBorn`, `update`,
`addBornMatrix`, `getEnergyParts`, `getADet`, `compute_fhb`),
`src/solv/model.f90` (`paramToModel`, `newBornModel`) and `src/solv/born.f90`.
The already-qualified thin caller is
`diagnostics/pede_solvent_response_20261013/frozen_distribution_candidate/implementation/REPLAY.f90`.

For a cavity atom list C and fixed atomic charges q, the native components are:

- `Gelec = 1/2 q^T B_e(C) q`, where `B_e = bornMat − diag(2 hbw)`.
- `Ghb = sum_i hbw_i(C) q_i^2`; `hbw` depends on solvent-accessible surface.
- `Gsasa = sum_i gamma_i SASA_i(C)`; it includes every cavity atom, charged or not.
- `Gshift` is the model/standard-state shift, not an atomwise environment energy.

The electrostatic matrix includes the native ALPB shape term
`Gshape = 1/2 keps * alpbet / aDet(C) * Q^2`. Thus an added zero-charge atom can
change core Born radii, local solvent access **and the global shape contribution**.
It does not become irrelevant merely because its charge is zero or its sphere is
contained within another sphere. The pairwise descreening and radius-weighted
shape construction are not a Boolean union of physical spheres.

Let L be the unchanged 192-atom local cluster and U its deduplicated full-source
union, retaining the 14 synthetic caps. The primary diagnostic is exactly

`S_M = Gelec(U, [q_M,vac, 0_environment]) − Gelec(L, q_M,vac)`.

Report `S_Ca − S_La` per source and the sample3-minus-sample0 difference. These
are component sensitivities, not affinities or class calls. Keep `Gshape`
separate from the remaining P16 electrostatic response. Record the other native
components for accounting, but **do not add full-protein Gsasa/Gshift or a second
complete solvation energy to the local composite**. `Ghb` also changes and is
outside this first electrostatic-only question. No permanent protein Coulomb
term, protein polarization, density relaxation, optimization or MM energy is
being represented. Zero environment charges are an algebraic isolation, not a
claim that this protein or its spectator Ca is chemically neutral.

Preserve the qualified factory: internal GFN2-xTB/ALPB, water, P16, `gsolv`,
298.15 K, 230 angular points, zero salt and native elemental radii/descreening.
The existing vacuum populations came from the separately recorded 300 K SCC
recipe. Do not change either temperature or use a common invented metal radius.
No charge normalization is permitted: each local-to-union comparison retains
its exact saved Q (Ca−3/La−2 nominal states). The native potential convention
and standard-state convention remain unchanged.

## 2. Physical source and boundary mapping

An actual matched, standard-residue-protonated source exists for each ExaA
two-Ca sample; missing full-source hydrogen preparation is not the blocker.
The authoritative inventory is
`workspaces/exaa_donor_response_20261014/inventory_v1/INVENTORY.json`, SHA256
`757a13523549697fd52cf9e94b3ba2ffdd01c463f3969fd8a5d01c1459d28343`.
Both source rows pin `source_preparation`, maps, state, full source and every
original native charge/energy receipt. Their context preparation's `source`
points to the exact `source_preparation/protonated.pdb`:

| Source | Protonated source SHA256 |
|---|---|
| sample0, control | `6b1c2f1c57842f1f6354fabc47d6a70ac78e927f0dd1ad9eaf063ebc01edbb41` |
| sample3, target | `4905b139bde56fdd14113092d0ce1f8382966e5404d45a197dc59fc289b9d536` |

The saved protonation receipt reports no detected missing residues, no repaired
missing atoms and both native Ca ions retained. Standard protein H are already
present. PQQ H were deliberately supplied by the pinned local carver, not
PDBFixer: retain the three existing `H_N1`, `H_C3`, `H_C8` coordinates and their
chemical mapping. This is the existing modeled protonation/oxidized-PQQ3− state,
not independent experimental state validation.

Construct U by starting with the exact local origin XYZ, then appending only
source atoms absent from its mapped physical atom identities. Match residue,
chain, insertion code, atom name and element; do not deduplicate by proximity.
Verify all mapped source coordinates against the pinned origin. The exchanged
source metal must be represented once by local atom0; the spectator remains
once at its source coordinates with its native Ca radius and diagnostic q=0.
Sample3's active ion is chain D; sample0's is chain B. Preserve every physical
source atom and all three existing PQQ H, rather than appending the whole PDB
and duplicating the local atoms. The proposed scout uses origins only, avoiding
a new moving-core/fixed-exterior coordinate construction.

**Exact limitation:** TBorn has one shared `nat/num/xyz` list for charge sites,
Born descreening, SASA/HB and shape. It has no independent charge-site versus
physical-cavity list. Keeping the charged link H therefore keeps their native H
cavity centers. Removing them loses charged sites; radius-zeroing or charge
redistribution would create another model. Restoring their omitted source bond
partners also produces expected short cap/source separations. These must be
mapped and reported, not passed off as an intact molecular geometry or removed
by an overlap filter. Accidental duplicate centers and failed mappings remain
preparation failures. A rigorous source-only cavity with all original cap charges
is **unsupported by this unchanged API**.

## 3. A finite cap control without arbitrary parameters

Use the native electrostatic matrix to report the physical/physical,
cap/physical and cap/cap quadratic blocks, with shape separated. Their sum must
reproduce the full native Gelec. This exposes charged-cap participation, not an
error estimate.

For a second, explicitly algebraic control, use the same real-atom charges and
set only the cap charges to exact zero. Evaluate this vector in (a) U and
(b) U with those now-uncharged cap centers removed. All real coordinates,
charges and radii are identical between the two controls, so their common
`Q_real` is conserved. Their difference isolates the cap centers' effect on the
real-charge electrostatic response, including shape. This removes no physical
source atom and changes no radius or fitted scale. It is not a chemically
assigned electronic state, must not be scored as one, and does not establish
the error of the full charged-cap calculation. A material cap-dependent effect
would leave the protein-shielding interpretation unresolved.

## 4. Proposed scout and numerical gates — unexecuted

First use **sample0 Ca origin**, chosen as the already-supported, fully mapped
control before looking at cavity energies. Reuse its accepted xTB/tblite
component-v2 provenance and fixed-charge native diagonal. Source/charge pins are
already joined in
`diagnostics/exaa_metal_position_20261014/solvent_input_audit/EVIDENCE.json`.

1. Replay its local cavity with the unchanged factory; compare against the
   archived native Gelec with the existing `1e−9 Eh` allowance and propagated
   eight-decimal charge precision. Verify all state/radius/coordinate pins.
2. Evaluate U and the cap-removed control. Repeated unchanged charge-vector
   calls must agree; native energy and quadratic/block reconstructions must
   close within `1e−9 Eh`. Preserve Q and report shape, cap and real-atom terms.
3. Repeat the complete U with one joint translation, one proper rotation and
   one atom permutation with charges/maps permuted identically. Apply the same
   `1e−9 Eh` numerical tolerance to Gelec; report residuals. SASA/HB rotational
   differences are separate diagnostics, not an implicit electrostatic pass.
4. Only after those gates pass, add sample0 La and sample3 Ca/La origins using
   the identical construction and controls. No accommodated structures,
   classifier calibration or expanded source set belongs to this scout.

This is four source/metal primary union energies and four cap-removed controls,
plus one local replay and three transformed/reordered union checks: **12 native
cavity updates**, with further q-cap-zero/block evaluations reusing those
matrices. Four old local energies are reused; there are zero new SCC, MACE or
DFT calls. Freeze actual transformation values, inputs and interpreter settings
before execution. A failed numerical or mapping gate stops expansion; do not
loosen it. Require the propagated charge-precision bound to remain informative
under the existing `1e−6 Eh` ceiling. No classifier success threshold is set.

## 5. Cost and implementation risks

The native model allocates `bornMat(N,N)`, `brdr(3,N,N)` and `dsdrt(3,N,N)` plus
all-pair lists and neighbor arrays even for zero environment charges. It has
quadratic memory/work components, not sparse scaling in the 192 charged atoms.
Its factory keeps the existing 35 Å Born-neighbor cutoff; global `aDet` still
uses the entire declared cavity. Do not replace the whole source with a new
arbitrary distance shell or change this cutoff in the first test.

The old wrapper also copies the full N×N matrix and writes every entry as text.
That output policy is inappropriate to inherit blindly for a protein. A thin
diagnostic caller can retain the needed charged-core matrix block, per-atom
cavity metadata and native scalar components without changing native energies;
this reduces output/copy overhead, not TBorn's internal quadratic allocations.
Run one fixture at a time initially and measure actual memory/time before
expanding; no throughput claim follows from static inspection.

The old wrapper assumes a `5e−9 e` rounding error on every input charge. Newly
appended zeros are exact, not rounded electronic populations. The correct
uncertainty bound must restrict the uncertain vector and matrix block to saved
charged-core populations (and to real-core populations in the cap-zero control).
Applying the old uniform error to thousands of exact zeros would artificially
inflate the bound. This changes uncertainty bookkeeping, not the native energy.

## 6. What the previous global-PB failure does and does not imply

The [global electrostatic challenger](../../global_electrostatic_20260916/REPORT.md)
failed its frozen surface/rotation and partition gates. Its approximately
18 kcal/mol partition shift was dominated by incomplete cancellation of vacuum
QM and direct core/environment terms; the reaction-field contribution to that
shift was about −0.63 kcal/mol. Those results do not prove that cavity shielding
is useless, nor validate this ALPB alternative. That model introduced permanent
charges, a different electronic partition and a numerical PB surface together;
the proposed fixed-q component test introduces none of those couplings or
surface solves. It still needs its own invariance, mapping and cap controls.

**Decision:** a small capped-union electrostatic sensitivity is executable and
potentially informative. A physically source-only, charge-preserving native
cavity correction is currently blocked by the shared charge/cavity atom list.
Do not implement or calibrate a correction from this feasibility report. The
next result should establish whether the finite exclusion effect is sizeable,
metal/source dependent and distinguishable from cap/shape sensitivity before
spending on a broader representation change.
