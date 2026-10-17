# Fixed-coordinate MM force attribution at the four LBT starts

Declared 17 October 2026 UTC under the active coupled-mechanics plan, before
new component evaluations. Root owns execution. This is a contained mechanical
attribution diagnostic alongside the three existing DFT scouts.

## Question and fixed denominator

The admitted saved force audit13032 finds the largest whole/hybrid force change
at Phe1 CA (physical ID `A/1//CA`, index3) in all four response starts. Its
archived MM group vectors identify retained bonded interactions as dominant,
but that group combines bond stretching, angle bending, torsions and CMAP.
Which unchanged terms produce this modeled restoring force?

Use every original baseline-region start, in order:7CCO/La,7CCO/Dy,7CCN/La,
7CCN/Dy. Bind the actual four-row force audit and each POINT, archived MM result,
execution/config, prepared XML, positions, physical map and interaction ledger.
No new coordinates, source choice, optimizer, chemical state or Hamiltonian.
The choice of this atom is explicitly diagnostic and based on the old audit;
it is not a confirmatory accuracy test or a new motion-selection rule.

## Finite calculation and gates

Run only the existing OpenMM8.5.1 Reference platform in its pinned Python
interpreter, one Slurm CPU with all numerical thread pools set to one.
Use the exact prepared native XML and all256 physical coordinates. No GPU,
DFT, MACE, GFN2, ALPB or new force-field parameter calls are requested.
All OpenMM energy/force queries are counted as new MM diagnostic work.

1. Verify hash-linked provenance, all four source/metal identities, exact atom
   order/positions, archived successful MM execution, zero constraints, declared
   force classes/groups and unchanged native/retained term parameters. Confirm
   the target is wholly MM. Retain the native source indices and unit conventions.
2. Replay the archived four MM groups and total full256-atom force. Require
   absolute energy residual at most1e-8 kcal/mol and maximum Cartesian force
   residual at most1e-5 kcal/mol/angstrom. These are numerical closure gates,
   not physical-accuracy tolerances; they are frozen before this execution.
3. Separate retained bonded interactions into harmonic bonds, harmonic angles,
   periodic torsions and CMAP by changing only in-memory force-group labels.
   Preserve all native parameters, particles and interactions. Require each
   regrouped sum to close to the original group/total at the same tolerances.
4. Attribute every retained bonded term touching Phe1 CA. Static lineage finds
   four bonds, fourteen angles and thirty-four periodic-torsion terms per state,
   with no CA-touching CMAP:52 terms/state,208 terms total. Preserve every term,
   including small, cancelling and zero contributions. Evaluate batches of at
   most31 target terms, reserving group31 for all unchanged remainder terms.
   OpenMM groups apply to Force containers: split those containers into exact
   single-term copies plus remainder containers, retaining every original term
   exactly once. Verify full-system closure for each batch; change no force
   equation, parameter, atom index or physical coordinate.
   Require summed CA vectors to close to each corresponding bonded class.
5. Report each term's native index, full atom IDs, parameter values with units,
   energy and target force. Include current bond lengths/angles with their
   reference parameters as descriptions of the modeled force, not assertions
   that the supplied structure is defective. Keep full-class force arrays and
   all four states. Sorted display does not change the denominator.

The finite upper bound is four original/relabelled-system pairs and two
local-term batches per state: at most4 contexts/state (16 total) and80
state queries/state (320 total), with actual counts recorded. No retries,
finite differences, sampling, optimization or new decomposition family is
launched automatically. A failed input, replay or closure keeps that state
unavailable and preserves its partial output. No survivor-only conclusion.

## Interpretation and delivery

This can identify the terms and parameter assignments responsible for the MM
force, and detect a concrete copying/counting/closure defect. It cannot decide
whether MM or MACE is more accurate, validate the molecular state, establish a
minimum, or predict metal preference. Large force or geometric displacement
from an MM reference alone is not proof of an incorrect parameter.

The independent whole-peptide DFT scout13051 remains the physical comparison;
Mex scouts12977/12978 remain unchanged. No parameter fitting, broader response,
new SCC recovery or production change follows automatically. Preserve original
failures and settings. Deliver executable source, pinned inputs, actual numerical
closure, complete term tables, measured cost and a concise mechanics conclusion.
