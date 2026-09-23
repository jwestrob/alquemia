# Structure-informed alternative starts — frozen mapping rule

Approved scope: the parent's [common-eight plan](../nikasha_parallel_pilots_20260922/PLAN.md)
and hash-pinned INPUTS.json. No new structures, chemistry, labels, DFT or default.
This is developmental transplantation; homology-independent transfer is not claimed.

Before energies, derive two starts per target from exact original source geometry:
Ca template1H4I, La template4MAE. If target and template share the declared biological
protein group, replace Ca with canonicalQ9Z4J7, or La with canonicalA0A3F2YLY8.
No template search by energy/outcome or alternate fallback after geometry rejection.

Use each target's already frozen four selected angular coordinates. Map by exact
canonical role, residue chemistry and chi number. Transfer the absolute source
N–CA–CB–CG, CA–CB–CG–CD/OD1, or CB–CG–CD–OE1 dihedral into the target's own
source graph. Require mapped quartet bonds and the existing rotation axis.
Use the shortest periodic delta. Terminal deprotonated Asp/Glu carboxylates permit
O-name exchange (periodpi) only when both archived fragment ledgers declare−1,
zero carboxyl H, and both oxygen atoms are present. Otherwise period2pi.
No foreign coordinates, atom identities, charges or proton states are copied.
A chemically absent role (notably Ca-template extraAsp) stays at targetq0 and is
explicitly recorded as unmapped; do not borrow a different Asp or invent it.

Reproduce every mapped source dihedral under actual target Kinematics, allowing
only the declared O equivalence. Reject unsupported mapping, zero displacement,
any exact start beyond the existing±0.8radian /0.8Å source domain, and chemical/
overlap guard failures. Never clip or rescale a template into the domain.
The same admitted starts and selected angular space are used for both metals.

Each admitted start gets one native float64OMOL SLSQP search per metal, with the
existing hartree objective scale,200iterations, fixed source-origin domain and
chemical guards. Preserve all failures. Compare final native energy to its own
actual starting energy; origin and old candidates remain separately available.
No assertion of a minimum or favorable biological score follows energy lowering.

At most32 searches; up to32 final geometries cross-scored for both metals means
at most128 nativeGFN2 calls and32 additional cross-MACE evaluations before exact
deduplication. Reuse only compatible actual outputs. Root's finite-candidate
adapter supplies the shared old+new pool and unchanged composite scoring;
there is no second scoring workflow. Both metals see every admitted geometry.
Report individual works, solvent/native contributions, rawR and old-band
transitions against static/adaptive, including all8sources and unsupported starts.
No new threshold or pilot-based classifier calibration.

First gate is source-derived geometry preparation only. Send actual admitted
counts/domain failures to root before any submission. No expanded domain or
replacement template is authorized by this plan.
