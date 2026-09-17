# Locate the error before constructing a DFT/MACE correction

Approved under Jacob's active MACE discriminator goal and standing pilot
authorization. Declared before the new local-core outputs. This responds to
both whole-protein solvent descriptors failing the five-structure development
screen. Do not calibrate that failed direct score on a larger panel yet.

## Question and inputs

Is the failed ordering already present in the local MACE/GB description, or
introduced by the whole-protein contribution? Does the independent hydrogen
preparation materially change that conclusion? Answer using the same five
consumed structures, labels, groupings and frozen metal/PQQ/water identities as
mace_global_benchmark_20260916/PLAN.md. Baseline and old inputs remain intact.

For each case prepare two explicit local geometries:

1. `archived`: exact existing baseline core XYZs from the pinned accuracy
   inventory. Reuse its actual DFT/CPCM endpoint results when verified. No
   quantum recomputation or inferred/missing energy is allowed.
2. `global_H`: same core atoms, source-to-core mapping, charge, caps and heavy
   positions; replace only physical source H coordinates with the already
   declared whole-protein preparation. Keep synthetic cap coordinates exact.
   Keep PQQ coordinates exact. Use the corrected retained-water H coordinates.
   This matches the new full system at every physical core atom. No matching
   DFT exists for this changed geometry; report it unavailable, not as baseline.

Use explicit generic v3 graph maps where available. For canonical PQQ use
its recorded fragment roles, atom inventories and source residue identities;
map source atoms against the archived unmodified protein coordinates before H
projection. The last documented fragment H is the covalent cap; no cap maps to
a physical source atom. Require a unique matching atom within0.002A, preserving
coordinates at their recorded precision. Mapping ambiguity or unsupported
chemistry fails. Do not recarve by distance or change core membership/charge.
Record each source H move, exact paired coordinates and all nonphysical caps.

## Fixed models and comparisons

Use primary medium analytic MACE-POLAR, its exact already validated checkpoint
and software. Large is not rerun in this diagnosis: the primary model's failure
is sufficient to investigate. Evaluate both local states for every endpoint.
Add the same native OBC-II frozen-monopole correction used globally, with all
parameters/limitations unchanged. Save component energies, densities and forces.
Neither local state is promoted or given an inherited aquo reference/class.

For each state report local Rvac and Rsolv. For global_H additionally report

    C_environment = R_full,solv - R_core,solv(global_H).

This is a matched low-level full-minus-core contribution. It changes cavity,
permanent interactions and learned charge response together; it does not
uniquely identify one of those as the cause. Retain its separate vacuum and
GB components. Report each of the three predeclared site contrasts before and
after this contribution. Record the archived-to-global_H local geometry effect
separately, rather than confusing it with an environment effect.

Compare archived local Rsolv with the actual archived DFT/CPCM R on the same
coordinates. Different electronic and solvent approximations mean this is a
combined local-method discrepancy, not a pure electronic-error diagnosis.
Do not add that discrepancy to the global score and call the result a valid
hybrid: its DFT geometry would not match. A coherent future candidate would be

    E_hybrid = E_DFTcore,CPCM(global_H)
             + E_MACEfull,vac + G_GBfull
             - E_MACEcore,vac(global_H) - G_GBcore(global_H).

That candidate needs matching new DFT endpoints, explicit CPCM/GB mismatch
qualification and appropriate partition checks. None are silently supplied in
this experiment. This experiment decides what subsequent calculation is useful;
it does not complete the active goal, even if a retrospective contrast improves.

## Finite inventory and physical checks

Five cases × two local geometries × La/Ca =20 new MACE energy/force calls,
followed by20 GB calls. Zero new DFT endpoints. No geometry search, gradient
finite differences, threshold fit or biological relabeling. Reuse the passed
analytic MACE kernel and GB numerical checks; all prepared pairs must preserve
physical source mapping and formal/electron parity accounting. Charge tolerance
1e-5e. Missing/failed results stay unavailable. Real-fixture parser/algebra,
coordinate and cache/recovery checks accompany the implementation.

One A5000/16CPU/64474MiB, existing manifested runner. Small cores are expected
to take seconds per endpoint plus initialization, so a few GPU-minutes total;
record the actual costs. No artificial project CPU/time budget. Keep artifacts
under workspaces/mace_local_correction_20260916 and compact records here.
