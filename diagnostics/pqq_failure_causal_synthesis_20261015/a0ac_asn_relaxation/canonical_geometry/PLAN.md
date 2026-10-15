# Canonical25 fixed terminal-amide geometry check

## Question and frozen population

Can the same intact 180° terminal-amide alternative pass existing geometric gates across the exact 25 consumed canonical sources? Use the completed Slurm4847 graph inventory, in its original canonical order. Every source already declares a neutral ASN through the fixed-core `anchor_asparagine` role. There is no proximity selection, label selection or remote-amide substitution.

For each paired Ca/La origin, rotate only OD1/ND2/HD21/HD22 around its mapped CB→CG axis. Keep CG fixed, with all other XYZ lines unchanged. This declares 25 source pairs / 50 possible geometry outputs. Preserve the exact original state, charge, water, cap and cofactor inventories using each source's actual context size. Preserve original four-mode IDs and report whether terminal χ2 was selected; do not reselect modes.

## Existing gates and minimal generalization

Reuse the executed A0AC preparation's source parser and original physical `Kinematics`. A focused copy of its fixed-flip geometry/contact checks substitutes the pinned role, context count, cap identities and state for A0AC-specific constants. Original helpers remain immutable. The contact helper changes only its omitted-source index boundary from 204 to the actual local count and its explanatory label.

Keep exact paired geometry and physical identities, source coordinate tolerance 1e−12 Å, physical-map replay 1e−8 Å, and covalent/rigid-group tolerance 1e−10 Å. Keep the existing relative severe-contact rule: fail a newly introduced distance strictly below 1.0 Å for heavy pairs or 0.55 Å when either atom is H, if its origin distance was at least that threshold. Include the full deduplicated protonated source, all physical H, fixed synthetic caps and unchanged carver cofactor H; report inherited contacts and omitted-source contacts separately. No boundary atoms are dropped.

This is a discrete endpoint check. As in the original fixed flip, it does not apply the old ±0.8-radian local domain or claim that a continuous rotation path is admitted. It does not build a rebased optimizer or evaluate Jacobian/force quality for a new search.

## Real regression and execution

Canonical A0AC La-conditioned seed1/sample3 has an exact existing positive: maps, state, source preparation, protonated source and both origin XYZ pins match the executed fixed-flip preparation. Run that source first in the same one-CPU Slurm job. Require byte-identical paired flipped XYZ and identical numerical geometric/contact evidence; only descriptive generic labels may differ. Failure stops further preparation and preserves all 25 result rows. There are no synthetic fixtures.

Then prepare the remaining sources, retaining each paired failure and reason in the full 25-row result. Failure of a source cannot create a partial metal pair. Metadata freeze only reads, hashes and copies inputs. All coordinate parsing, rotations, contact checks and the real regression run only inside actual Slurm, using the pinned numerical interpreter, one allocation CPU and one thread. No node, memory or walltime override is prescribed.

## Outputs and limits

`execution/PREPARATION.json` retains all 25 rows; successful rows contain paired `matrix[metal].amide_flip` XYZ pins, original maps/state/source pins, geometry admission and original mode membership. Each source preserves geometry/contact evidence and contact identities. `A0AC_REGRESSION.json` and `EXECUTION.json` record actual qualification and costs, including failures.

No optimization, energy/model call, score, class, calibration or production change is performed. A passing endpoint proves neither an equilibrium conformer nor the applicability of a subsequent search. This stage ends at geometry preparation; geometry admission and the regression are UNRUN until the actual job finishes.
