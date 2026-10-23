# Baseline operator export: fixed consumed inventory

`FIXTURE_INVENTORY.json` identifies **54 distinct zero-field systems**, without
symmetry replication or selection by an observed error:

| Fixtures | Count | Existing authority |
| --- | ---: | --- |
| Original water development | 10 | Selected17676 compact DATA: energies and forces for the original25 checks |
| Literal extracted W8/W9 metal–water dimers | 34 | Saved78-moment result: energies, components and per-atom induced dipoles |
| Near/far metal–two-water trimers | 8 | Actual18343 trimer model result: energies, components and induced dipoles |
| Bare La/Dy ions | 2 | Saved78-moment result: actual zero energies and dipoles; retained as export controls |

The four conditional-susceptibility zero states reuse W9 water1 dimers and W9
far(1,3) trimers, for both metals. Their16 signed-field results and prior zero-Dy
identity are cached checks, not additional baseline geometries. A field changes
the right-hand sides of the unchanged linear model; it does not require another
zero-field matrix once the applied-field coupling and units are verified.

## Exact joins and first gate

The first fixture is `aquo9_central_La__dimer_01`, the existing `La_dimer` field
state, with parent atom indices `[0,1,2,3]`. Its input XYZ is pinned7eff5303,
reporting key96ca811e and parameter filee1cc2bcc; full hashes and paths are in the
artifact table. Its admitted baseline energy is −76.03570895kcal/mol. The
original source key79b82234 differs only by the recorded reporting keywords
`polar-print` and `save-uinduce` plus comments/blank lines. The future export adds
only the proposed no-argument `RESPONSE-MATRIX-EXPORT` flag.

The inventory checks exact numerical Cartesian coordinates, atom order, types,
bonds, parameter bytes and all active key tokens. It preserves full input-file
hashes and original producer paths even when bytes match at different paths.
Only comments/whitespace and the two named reporting flags are ignored for
deduplication. The separate field-input joins verify their actual appended
`EXTERNAL-FIELD` vectors. All54 signatures are distinct. La/Dy inputs share the
same selected17676 parameter lineage; there are exactly two effective key
token lists, one per metal. None of the34 literal dimers is replaced by the
published isolated-water geometry or a rotated representative.

After actual scout parity, the remaining inventory is53 instrumented calls.
Two induction snapshots produced by one `analyze EM` invocation count as one
native call. Planck's fresh uninstrumented comparator, if required by the scout,
is an additional call under that separate protocol. This inventory releases
neither stage and does not bind a new executable.

## Cached evidence and limits

All54 baseline energies are available. The ten original development cells have
forces;48 systems also have admitted per-site `.uind` data from the chosen
moment/trimer authorities. The six original isolated-water radial systems have
no admitted `.uind` in these chosen authorities. The full A, Ed, Ep and separately
identified converged u/v arrays are new export data for every row; an old
`.uind` file alone does not establish an independent v array or the operator.
Keep native ipole/site ordering and original units explicit when joining them.

The46 existing metal-absent source rows (34 waters,4 shells,8 water pairs) stay
cached, with exact inputs and result pointers retained. All seven-term trimer
cycles close on the exported inventory plus these cached rows. Ion terms remain
explicit, including the N-fold subtraction in P; the proposed oxygen-only
correction vanishes on ions, so their cached energies/moments need no candidate
reevaluation. No cancellation of parent-dependent quantum fragment terms is
assumed by this classical-model inventory.

A fixed-coordinate A/Ed/Ep export supports exact finite *operator* response.
It does not supply coordinate derivatives of the operator or driving fields,
and cannot by itself qualify changed forces against the original25 gates.
The baseline energy-conjugate total moment and both native induction solutions
must be retained; do not transfer an induced-array convention to a changed
functional without checking it. Newly consumed moments and susceptibilities
remain descriptive until a prospective loss is declared. Original acceptance
thresholds and every reserved/held geometry remain unchanged.

Construction used one local thread,1.65s wall/1.27CPU-s,496 file-pin checks and
zero native or quantum calls. This was saved metadata extraction, not a new raw
scientific replay, fit, force test or operator solve.
