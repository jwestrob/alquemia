# Four conditional W9 radial reference cells

Implement the parent declaration only. Use Planck's pinned four-geometry manifest,
`aquo9_{minus,plus}_{La,Dy}`, with rigid ±0.10 Å translations of every intact
water along its original metal–oxygen direction and a fixed metal. Preparation
reads metadata and the two previously admitted central W9 references; it opens no
new validation outcome or fit trajectory and makes zero native calls.

Copy the original central native cells into a four-cell panel with only geometry
and identifier changes. Preserve the original `REFERENCE.py` adapter, external
ECP/basis assets, PBE0-D4/def2-TZVP/def2-J/RIJCOSX, DefGrid3, VeryTightSCF,
PModel/NoAutostart, MaxIter300 and EnGrad. Each cluster has charge+3,98explicit
electrons and effective multiplicity1; physical La/Dy multiplicities remain1/6
with core sizes46/55. The template deck is exactly the matching central deck;
the geometry file carries the displacement.

Prepare four templates and a finite root-release template with candidate/review
fields null. Native release requires terminal job17676 with every joint control
passing, its frozen six selected parameters and the shared independent
`selected_ten_raw_review_passed` contract. The release binds exact sources,
geometry/native manifests and four unique output directories. Each native cell
uses the existing four-MPI launcher, one thread/rank, standard_expand and default
memory/time; MaxCore remains floor(0.75×actual granted MiB/4). Exclusive task
receipts prohibit an automatic rerun in the released directories. Maximum four
new native invocations; no central reruns or fallback.

The saved collector retains all four slots and both central anchors. Missing or
failed results remain unavailable. It replays the original native admission,
full gradients, native recipe and source/input identities. It retains numeric
SCF lines, warnings and final COSX corrections, and flags changes in the central
energy-stage signature. Such metadata does not establish analytic derivative
accuracy; this finite radial transfer test supplies no infinitesimal derivative
qualification. References stay on the same provisional uncorrected ECP-PBE0-D4
surface. No model scoring or parameter fitting occurs in this adapter.
