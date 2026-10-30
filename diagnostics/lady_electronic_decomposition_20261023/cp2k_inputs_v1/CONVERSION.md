# Archived quartet: CP2K basis and ECP conversion

Prepared 23 October 2026. This package performs source/library conversion only.
It has made **zero molecular calls and zero job submissions**. The four archived
geometries and original ORCA inputs remain unchanged.

## Files and use

- `assets/BASIS_LADY`: orbital basis for H, C, O, La and Dy; common basis name
  `LADY-ARCHIVED`. SHA256 `2ab7190ac0bf52f389fed017672c49f747d75d6d0e92f61dbf4e2f5981b0eb66`.
- `assets/POTENTIAL_LADY`: literal La46/Dy55 ECPs and the installed CP2K release's
  H/C/O all-electron records. SHA256
  `ba87181563f5ca2b8e6e4ff30f43d274752b9c663325905ee6c1dc41c7beccd5`.
- `assets/CONVERT.py`: deterministic standard-library-only converter and exact
  token/metadata checks. Existing different outputs are refused.
- `assets/PROVENANCE.json`: source pins, all converted records, four geometry
  joins, library evidence and explicit limitations. SHA256
  `ff923ac1a23676e0cd799a989e17271bf30e10aefb74f7a503a50899b82553b4`.

Use the two files as `BASIS_SET_FILE_NAME` and `POTENTIAL_FILE_NAME`. Every KIND
uses `BASIS_SET LADY-ARCHIVED`. La uses `POTENTIAL ECP-LADY-La46`; Dy uses
`POTENTIAL ECP-LADY-Dy55`. Both metal KIND sections require
`ECP_SEMI_LOCAL TRUE`. H/C/O use `POTENTIAL ALL`. No GTH replacement is involved.
Root owns the remaining deck, numerical settings and execution declaration.

## Exact metal conversion and reader convention

The literal `inputs/assets/{La,Dy}_basis.inc` files supply all metal orbital
exponents, contraction coefficients and ECP terms. No coefficient arithmetic,
rounding, decontraction, rescaling or basis substitution occurs. Each original
NewGTO shell becomes one CP2K set with one contraction; the leading set label
is 1, following the BSE CP2K formatter. Default orbital normalization 2 is
required: CP2K normalizes primitives and contracted functions. The spherical
shell space remains 6s4p3d2f, or 47 orbitals, for each metal.

The installed `src/atom_types.F::read_ecp_potential` reads the removed core
count from `nelec`, then the local UL block, followed by consecutive nonlocal
S/P/D/F blocks. Each numeric row is `n exponent coefficient`. The original
zero-coefficient g channel is retained literally as UL; it is not an extra
nonlocal g projector. Since that local term is zero, no subtraction changes
the S/P/D/F coefficients. The restricted converter rejects a nonzero local
channel rather than generalizing this operation silently.

All original radial tokens are **n=2 and stay 2**. CP2K's LIBGRPP evaluates
`coefficient * r^(n-2) * exp(-exponent*r^2)`, with exponents in inverse squared
Bohr. A cross-format check is the published ORCA SD(10,MDF) vanadium example:
its terms, including n=2, match the installed `data/ECP_POTENTIALS` V record.
This resolves the input convention without interpreting the ORCA manual's
schematic radial exponent as an instruction to add or subtract two.

La and Dy retain effective nuclear charge 11; the M3+ fragments consequently
have 8 explicit electrons. The isolated, reviewed Dy55 reader patch is already
part of the pinned installation. This conversion makes no further source change.

## H/C/O: explicit source-equivalence evidence and limit

All four saved native outputs identify def2-TZVP and print its contraction
patterns, but do **not** print H/C/O primitive coefficients. Accordingly, this
package does not claim byte-level recovery of those internal ORCA tables.

The [official ORCA 6.1 basis documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/basisset.html)
identifies the keyword as Karlsruhe def2-TZVP and describes the BSE/EMSL basis
library provenance. The converter uses the
[BSE def2-TZVP H/C/O data](https://www.basissetexchange.org/api/basis/def2-tzvp/format/json/?elements=1,6,8&version=1),
with its Weigend–Ahlrichs 2005 reference, DOI `10.1039/B508541A`.
The saved BSE version 0 (original BSE) and version 1 (Turbomole 7.3) contain
identical H/C/O `electron_shells` objects. The separately downloaded official
BSE ORCA and CP2K exports contain exactly the same exponent/coefficient tokens.
All download URLs, hashes and bytes are saved in `assets/BSE_DOWNLOADS.json`.

The resulting H pattern `311/1` and C/O pattern `62111/411/11/1` match every
relevant archived output. Orbital counts are H 6, C/O 31 and La/Dy 47, yielding
391 for both aquo8 geometries and 490 for both hydrated-acetate geometries,
exactly matching the old native summaries. This is documented named-library
equivalence plus metadata agreement; it is not a cross-engine energy test.

## Scope retained

The quartet remains effective singlet throughout (88 or 110 total explicit
electrons); this does not qualify the physical Dy spin or spin–orbit treatment.
Full geometry/input/output pins are recorded per state. ORCA def2/J and the
metal AuxJ records belong to its RIJCOSX numerical approximation and are not
silently used as orbital basis functions. CP2K's numerical integration,
exchange treatment, fragment SCF and D4 bookkeeping need their separately
declared settings and actual checks. Formatting/source checks here establish
none of those numerical results.

To reproduce only the conversion:

```bash
/home/jwestrob/.pyenv/versions/3.12.0/bin/python diagnostics/lady_electronic_decomposition_20261023/cp2k_inputs_v1/assets/CONVERT.py
```
