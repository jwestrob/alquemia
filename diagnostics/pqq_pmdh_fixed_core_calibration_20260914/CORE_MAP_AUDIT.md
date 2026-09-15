# PQQ-MDH fixed-core mapping audit

Date: 2026-09-14

Machine-readable selectors and source hashes are frozen in `core_map.tsv`.
This audit prepared no QM inputs and launched no calculations.

## Result

All 25 frozen calibration-panel structures have an unambiguous common mapping for:

- the metal and complete PQQ cofactor;
- the conserved metal-ligating Glu and Asn;
- the conserved catalytic Asp;
- the XoxF/Ln-specific extra Asp, when present; and
- the cationic Arg/Lys that hydrogen-bonds to the catalytic Asp.

The fixed protein core is recorded explicitly per structure in the
`fixed_protein_core_selectors` column. The metal and PQQ must be added to that
protein list. Execution should consume these frozen selectors and fail closed
on a source-hash or residue-identity mismatch; it should not rediscover the
core from a distance cutoff.

## Uniform source layout

All 25 source CIFs have the same entity layout:

- chain A: one contiguous protein, numbered from 1 with no insertion codes;
- chain B, residue `LIG_B 1`, atom `LA`: the only metal;
- chain C, residue `LIG_C 1`: the only PQQ-like nonpolymer; and
- no waters, alternate conformers, sequence-number gaps, or other
  nonpolymers.

The 14 `Ca-verified` labels describe the biological class. Their source models,
like the 11 `La-verified` models, physically contain La and PQQ; every filename
is consequently `*-pqq-la_model.cif`.

Every PQQ has the same 24-heavy-atom generic naming graph. All 25 pass the
current `protenix_generic_pqq_v1` connectivity validator and can be represented
as the production `pqq_ox_3minus_v1` microstate (`C14H3N2O8`, 27 atoms after
restoring its three invariant hydrogens). After historical normalization,
`B:LIG_B1` becomes `B:LA1` and `C:LIG_C1` becomes `C:PQQ1`; the chain-A residue
selectors are unchanged.

## Mapping logic and evidence

For audit purposes, the conserved Glu and Asn were the unique closest residues
of their respective types to the metal. The separation from the next candidate
is large: the next Glu is at least 10.95 A from the metal and the next Asn is at
least 6.59 A away.

The catalytic Asp is the Asp immediately following Trp in the unique
`TP...WD.[DAST]` motif. The residue two positions after that Asp is:

- Asp in 11/11 La-verified controls (`W-D-x-D`); or
- Ala, Ser, or Thr in 14/14 Ca-verified controls.

Only an Asp at this `+2` homologous position belongs to the fixed core. The
nonacidic Ca-control homologs are recorded for validation but are not included.

The catalytic-Asp cationic partner is also unambiguous. It is Arg in 23
structures and Lys in Q4W6G0 and Q8GR64. Every contact uses catalytic Asp OD2
and Arg NH2 or Lys NZ. Asp-O--N distances are 2.584--3.333 A, while the nearest
alternative cationic residue is at least 7.36 A from the catalytic Asp. The
partners occur 23, 27, or 28 residues downstream and their N atoms remain
4.329--5.102 A from the metal, confirming that they are second-shell rather
than direct metal ligands.

Observed metal distances across the panel are:

| Role | Range (A) |
|---|---:|
| conserved Glu side-chain O | 2.156--2.621 |
| conserved Asn OD1 | 2.370--2.848 |
| catalytic Asp side-chain O | 2.433--3.522 |
| Ln-specific extra Asp O | 2.410--2.648 |
| catalytic Asp O to partner N | 2.584--3.333 |

The old normalized and protonated structures were checked in both independent
workspace generations (`00_Lav_`/`01_Cav_` and `colinpqq_la_`/`colinpqq_ca_`).
All 50 workspaces preserve every residue selector and all relevant heavy-atom
coordinates. The active Asp/Glu residues remain deprotonated; Arg and Lys have
their canonical cationic protonation; and the legacy PQQ remains heavy-atom
only. No mapping or protonation anomaly was found.

## Legacy cutoff bias discovered

The legacy 3.2 A carver did not use this chemically fixed core. It omitted the
catalytic Asp from:

- 1/11 La controls: C5AXV8; and
- 8/14 Ca controls: ATQ70401.1, O24759, P12293, P15279, P16027, P38539,
  Q60AR6, and Q9L935.

Both independent legacy carve generations reproduce exactly those omissions.
Accordingly, the historical 25/25 separation was obtained with unequal,
geometry-dependent chemical models: 10/11 La controls but only 6/14 Ca
controls retained the catalytic Asp. That does not invalidate the historical
observation, but it means the old decision threshold cannot be transferred to
the fixed-core method.

A 3.6 A shell happens to recover the Glu/Asn/catalytic-Asp/extra-Asp donor set
for these exact 25 structures without adding an off-core protein donor. It
still omits every cationic partner, and it provides no guarantee on future
structures. The explicit map is therefore required.

## Charge consequence of the proposed core

With fully deprotonated oxidized PQQ(3-), deprotonated Glu/Asp, neutral Asn,
and cationic Arg/Lys, the expected total charges are:

| Biological class | Protein fragments | Scaffold | La arm | Ca arm |
|---|---:|---:|---:|---:|
| La control (Glu + Asn + catalytic Asp + extra Asp + partner) | -2 | -5 | -2 | -3 |
| Ca control (Glu + Asn + catalytic Asp + partner) | -1 | -4 | -1 | -2 |

These charges should be derived from and recorded in each carve manifest, not
accepted as unchecked input constants.
