# Matched water removal: six fixed central controls

Declared 21 October 2026 before preparation/execution. Root owns review and
submission. This task prepares inputs and runner bindings only; no molecular
calls, scheduler submissions, geometry optimization, or parameter fitting.

## Question and finite inventory

The existing hydrated central fixtures differ from the original dry dimers in
both water content and donor placement. Separate these effects by deleting
only water atoms from the **six existing hydrated central fixtures**:
monodentate acetate, bidentate acetate, and acetamide, each with La and Dy.

Retain the metal and every ligand atom, including hydrogens. Preserve surviving
symbols, coordinate text, order and chemical state exactly. No rotation,
translation, hydrogen adjustment, relaxation or geometry rebuilding is allowed.
Acetate complexes have charge +2; acetamide complexes have charge +3. Physical
multiplicities remain La 1 / Dy 6. Removing seven waters from mono-acetate or
acetamide leaves 8 or 10 atoms respectively; removing six waters from
bi-acetate leaves 8 atoms. Every frozen-f ECP cell has 40 explicit electrons
and effective multiplicity 1, distinct from physical Dy sextet metadata.

Exactly six native ECP-PBE0 energy/gradient cells and six MACE energy/force
evaluations are proposed. No new aquo8 or hydrated calculation is in this
inventory. Reuse the eight admitted existing ECP endpoints and their MACE
counterparts (six hydrated central plus two aquo8) through exact result pins.

## Unchanged methods

Bind the admitted `reference_chemistry_v1/REFERENCE.py` executor and original
ECP-PBE0-D4 recipe: def2-TZVP/def2-J, archived La ECP46 / Dy ECP55 orbital and
auxiliary assets, RIJCOSX, DefGrid3, VeryTightSCF, PModel, MaxIter300, EnGrad,
NoAutostart. Do not change native admission gates or reinterpret successful
termination as physical-state or reference-accuracy qualification.

Use root's unchanged `mace_chemistry_v1/RUN.py`, admitted helper, runtime,
checkpoint and float64 model configuration. Pass physical La1/Dy6 inputs.
Do not switch checkpoints or equate MACE's state input with the frozen-f
reference Hamiltonian.

Suggested execution matches the admitted controls: native singleton jobs with
four MPI ranks/one thread each, `standard_expand`, scheduler default memory and
walltime, MaxCore derived from the actual memory grant; one finite MACE job
with one H200 and the existing 12-CPU launch wrapper. Resource choice and
submission remain root-owned. No native call is made in preparation.

## Balanced comparisons

For each method and donor define, using the **same saved aquo8 reference**:

`Cwet = E(Dy,wet) + E(La,aquo8) − E(La,wet) − E(Dy,aquo8)`

`Cdry = E(Dy,matched dry) + E(La,aquo8) − E(La,matched dry) − E(Dy,aquo8)`

Each expression is an element- and charge-balanced La/Dy exchange. Then:

`Cwet − Cdry = [E(Dy,wet) − E(Dy,matched dry)]
              − [E(La,wet) − E(La,matched dry)]`.

This isolates the metal-dependent effect of adding the specified six or seven
waters to the fixed metal–ligand geometry. Explicit free-water energies cancel
between metals because the same number of chemically identical waters is added.
The equivalent balanced reaction is
`La(wet) + Dy(matched dry) -> Dy(wet) + La(matched dry)`.

It is a many-body total-energy response: water–water, metal–water,
ligand–water, polarization and dispersion interactions can contribute. It is
not a unique ligand–water energy, equilibrium hydration free energy, or protein
affinity. It does not include entropic or chemical-potential population terms.

Donor-placement effects can later be evaluated against available original dry
endpoints with the same reference convention. Do not subtract a historical
water-monomer-referenced contrast directly from an aquo8-referenced contrast;
reassemble the common cycle from admitted endpoint energies first. No missing
historical endpoint is authorized as a new call by this six-cell preparation.

## Saved validation and stop rules

Check exact parent pins, central-only six-cell denominator, paired La/Dy
coordinates, surviving line equality, deleted-water composition/count, total
charge, physical/effective multiplicities and explicit electrons. Validate
rendered decks through the frozen executor's saved-only `validate` function.
Pin all source, model configuration, basis and geometry inputs. Preserve
full force arrays; the old metal-force projection is not the derivative of a
rigid ligand translation. That derivative is minus the sum over all ligand
forces dotted with the saved unit translation axis.

Root must review the prepared bundle before submission. Future admission uses
the unchanged runners, preserves each failed attempt and missing row, and never
substitutes a different solver/state/result. The specialized original56-cell
native collector cannot collect this six-cell experiment; use a separately
declared saved collection with complete four-term denominators.
