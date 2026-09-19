# A reference reproducibility issue, without a reversed class

The completed comparison has one DFT score outside the predeclared 0.01
kcal/mol reproduction tolerance: Q9L935 changes from 9.6498206631937 to
9.774587501023811, a +0.12476683783011211 kcal/mol difference. Both classify
as Ca under the frozen bands. MACE reproduces every archived score exactly.
The combined qualification gate remains **false**; no tolerance was widened.

## Actual endpoint evidence

| Endpoint | Archived energy / hartree | Fresh energy / hartree |
|---|---:|---:|
| Ca | -2953.523479200126 | -2953.523479202328 |
| La | -2307.461311336739 | -2307.461510167553 |

The change is almost entirely the La endpoint. Source and fresh XYZ hashes
are identical. ORCA is 6.1.1 in both outputs. The echoed inputs differ only
in MPI ranks (15 archived, 16 fresh); method, grid, CPCM, charge, singlet,
NoAutostart and memory setting are unchanged. Both terminate normally and
print SCF convergence, after 9 and 23 cycles respectively for La.

Both La outputs report final density and DIIS changes above their printed
tolerances. This alone does not contradict ORCA's convergence stamp:
the default convergence check uses energy criteria rather than requiring every
printed criterion. See the [ORCA 6.1 SCF manual](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html).
Numerical convergence sensitivity is a plausible explanation, not a proven
unique cause or evidence that either wavefunction is the stable ground state.
No tighter-SCF calculation, stability calculation or rescue rerun was launched.

Separately, A0A3F2YLY8 shifts by -0.00024073708348737455 kcal/mol and crosses
its exact La calibration boundary into inconclusive. It remains within the
reproduction tolerance. The fresh DFT panel therefore has 24 correct calls,
one inconclusive, zero opposite-class calls, and 24/25 reproduced scores.
These are different issues; neither is a MACE prediction failure.

The completed MACE runs retain 25/25 reference calls, and the timing advantage
is measured. The strict combined benchmark is not fully qualified because
of the DFT reproduction failure. Do not erase these useful findings or claim
an accuracy improvement from them. Retain the production default and frozen
records; discuss a versioned convergence check before changing the protocol.

Exact source hashes, input differences and convergence tables are in
`workspaces/mace_pqq_utility_20260918/DFT_reproducibility_v1.json`.
