# 1F6S normalized 52-atom response fixture — mapped and ready

**The archived 52-atom amide representation and its exact protein field are
coherent enough to prepare this explicitly conditional response diagnostic.**
Four real-fixture checks, including the existing six-cell dry-run, pass. No
molecular calculation or job was launched by this task.

Manifest:
`workspaces/metal_environment_response_20260926/transfer_preparation/1F6S_normalized52_v1/INPUTS.json`.
This is a new named use of the existing normalized NMA representation, not the
original 40-atom repair. The 40-atom core is neither modified nor combined with
the 52-atom field. Root decides execution and scientific interpretation.

## What the mapping review verified

All 52 core atoms are accounted for by the archived physical source or explicit
synthetic caps. Every mapped source coordinate agrees with the **normalized**
physical preparation, including hydrogens. No raw/normalized H mixture is made.
Eleven caps reproduce their coordinates from actual retained/omitted source
anchors to a maximum 4.93e−11 Å error. Both physical Jacobians are saved for every
cap; finite-difference checks give at most 6.56e−10 absolute error against the
predeclared 1e−8 mapping tolerance.

Both backbone donors contain the real peptide chemistry and actual bonded
neighbors: Lys79 C/O with Phe80 N/H, and Asp84 C/O with Leu85 N/H. The C–O,
C–N and N–H connections were verified against the archived source bond graph;
they are not inferred merely from sequential residue numbers. Asp82, Asp87 and
Asp88 carboxylate states and the original proton inventory remain unchanged.
Both coordinating waters211/212 retain all three atoms in the electronic region.
No cofactor parameter or additional water was invented.

The field exactly regenerates from its archived physical indices and charge
weights. A is byte-identical to the earlier native point-charge file, and both
metal XYZ files are reused byte-for-byte. Its **1,880 charges sum to −4 e**.
The mapped QM projection support and environment identities together cover the
1,932-atom conditional physical preparation. Support includes physical link
anchors as well as directly retained QM atoms; removed direct field charges are
not physical atom deletions.

This is the archived ff19SB background, not AMOEBA. Its general peptide-boundary
ledger removes the actual mapped support and restores each retained residue's
prescribed charge over recorded bonded exterior recipients. Recipients can
include H atoms. We preserve this existing policy and its exact increments,
without pretending it is the PQQ-specific CA-removal/equal-N,C prescription.
Charge closure alone does not qualify that electrostatic boundary approximation.

Ca uses core charge −1, multiplicity1, 214 explicit/all electrons. La uses charge0,
multiplicity1, 250 all electrons and 204 explicit electrons with the native
46-electron ECP. Both endpoints have identical nonmetal atom inventories and
coordinates. Combined core/field charges are Ca−5 and La−4. Runtime electronic
state checks remain required for any later calculation.

## Same deterministic A/B rule

The nearest complete exterior hydroxyl is **Thr86**, oxygen-to-core-heavy distance
4.64009 Å. B rotates its actual HG1 by +10° about CB→OG1. Only point-charge index
1292 changes coordinates, by0.15509 Å. O–H length0.96 Å and C–O–H angle are
preserved; the closest nonbonded H-heavy distance after rotation is1.99772 Å,
above the1.2 Å screen. Charges, core geometry, protonation, waters and all other
physical coordinates remain fixed.

This is a more remote perturbation than the two PQQ cases and may produce a
small signal. A near-zero response would not demonstrate missing metal sensitivity
at direct donors. We will not strengthen or replace the perturbation because its
energy response is inconvenient. B is a modeled orientation, not a crystallographic
observation or equilibrium population.

## Evidence and limits

The existing evidence stratum is `qualified_qualitative_strong_site_affinity`,
group `bovine_alpha`, consumed method-development material. That biological
information motivated a relevant fixture, not a target for this environmental
response sign. No new classification, affinity, Kd, occupancy or generalization
claim follows from these endpoints. A conditional protein electronic response is
not a binding observable.

The source's archived protein radial-H and water-H normalization are retained
without new optimization. The field remains fixed-charge; the native electronic
region can respond. Numerical qualification remains unestablished after prior
rigid/SCF issues. Full hybrid forces/energies need their own interaction ledger,
parameters and boundary treatment; this preparation supplies none by default.

## Reproduce and resume

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  diagnostics/metal_environment_response_20260926/transfer_preparation/1F6S/prepare.py \
  --repository "$PWD" \
  --output workspaces/metal_environment_response_20260926/transfer_preparation/1F6S_normalized52_reproduction
```

`check.py` takes `--repository`, `--inputs` and a new JSON `--output`. Four checks
passed: paired geometry/electron parity, unchanged charges with only the selected
H moving, all11 physical cap Jacobians, and six-cell preparation/dry-run through
the existing runner at6×57MPI. These are real-artifact geometry/parser tests,
not scientific integration tests. No new force-field System was built.

Current input files and generating implementation are frozen. New preparations
must use a new output directory. Root owns any submitted jobs and collector;
this task used zero molecular calls and zero allocated compute. Local preparation
and test process CPU/RSS were not formally measured. No production, library,
remote push or external email changes occurred.
