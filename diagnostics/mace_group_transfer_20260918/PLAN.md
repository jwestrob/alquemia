# Declared: typed-contact whole-protein group transfer

Status: authorized under Jacob's active discretionary MACE goal; declared before
new model evaluations. Previous goal turn made progress: completed the24+24
charge-group study, found a qualified PQQ gain but1/7 total directions, and
identified the vacuum component as a new lead. No old result is reclassified.

## Model and physical inventory

Protocol:intact_POLAR_medium_typed31_group_vacuum_compatibility_v1.
An empirical fixed-geometry electronic compatibility descriptor, not a binding
free energy or gas-phase substitute for solvation in an affinity calculation.
Same medium POLAR checkpoint, full atom inventory, learned weights, formal
protonation and analytic kernel. Same task-local grouped Fukui normalization
at all3 density-restoration stages. No OBC-II/CPCM term, aquo offset, threshold
fit, geometry optimization, training or new DFT. Vacuum selection is motivated
by already inspected results, so the original7 structures are development.

Replace carve-dependent group selection with one uniform whole-protein rule.
Use existing typed protein/PQQ/water direct-donor definitions and3.1A cutoff.
A protonated histidine N with a bonded H cannot join as an N donor. Group the
entire donor residue with each contacted metal. For a coordinating backbone
carbonyl, also join the actual amide N's complete residue using recorded C-N
bonds, never residue-number arithmetic. Join disulfide partners. Merge all
overlaps through this source graph, including metal sites sharing residues.
PQQ remains one complete -3 group; waters are complete neutral groups. A
retained water not contacting a metal remains a separate neutral group.
Unselected protein residues stay complete individual formal-charge groups.
All metal sites are processed together, independently of the selected endpoint.
Background metals remain Ca2+; the selected endpoint is Ca2+ orLa3+.
Replay ff19SB template residue charge sums, including actual termini/acetyl.
No force-field atomic charge is inserted in the learned energy. Unsupported
chemistry, missing source mapping, duplicate atoms or nonclosure fails openly.
Preserve all heavy and hydrogen coordinates, assembly and existing water union.
Canonical atom order is lexical source ID within a protein; no atom is removed.

This is a new group policy. Recompute all7 old structures to measure the policy
change; do not mix new-family scores with the old carve-defined scores as if
identical. Formal group charges remain unvalidated density constraints. The
old vacuum grouping failures remain recorded, irrespective of the new result.

## Fixed scope and cost

Sources:old7 (PQQ1H4I/4MAE,alpha1F6S/6IP9,GGR1GLG/2FW0/2FVY), all22Khoury
sites (A0A7six,HEW5eight,RTXeight), parvalbuminCD/EF and aequorinEF1/EF3/EF4.
These are9 biological groups,8 with directional/supporting labels. All sources
and labels have been inspected previously under other models; none is called
blind. Khoury ITCLaKd versus CaCDfolding thresholds remains supporting evidence;
Ca threshold is not a CaKd and does not assign individual site affinities.

55newMACEcalls:14old-structure endpoints;one all-Ca reference for each of5new
proteins plus27single-La substitutions(32);6GGR endpoints with the old complete
connected-core residue set additionally merged;RTXfirst-siteCa/La rotations
and first-siteLa reversed-atom permutation(3). Source-equal all-Ca cases share
one evaluated reference with explicit atom/group equality checks. No other
endpoint reuse. Approximately20–30GPUminutes expected from previous actual
9–59s medium endpoints; oneA5000/16CPUs/64474MiB. Count failures and overhead.
No project time/CPU budget; this finite manifest defines scope, not a cost cap.

## Frozen comparisons and checks

R=E(Ca)-E(La), exactly one eV-to-model-kcal conversion. No absolute classifier.
All7old relative directions retained. Each Khoury arithmetic mean of all sites
must exceed each of3GGR scores:9supporting comparisons, not9independent groups.
Both parvalbumin sites compared separately with all3GGR:6supporting comparisons.
Total22 directional comparisons, reported by stratum. Aequorin remains the
ordered3-site vector with no sitewise pass/fail label. Never select a favorable
site or structure. No threshold refitting or label/protonation changes.

Numeric:energy rotation/permutation<=.01modelkcal; force<=.001eV/A; group/total
charge<=1e-5e; Fukui abs(sum)/sum(abs)>1e-10 at every group/stage. Paired atom
and coordinate equality;all-Ca source/group equivalence across every site.
All3GGR grouping contrast changes retain the2modelkcal tolerance. Raw accuracy
is always reported even if representation qualification fails. A positive
result requires additional compatible calibration and external testing before
any promotion; this experiment cannot by itself complete the broad goal.
