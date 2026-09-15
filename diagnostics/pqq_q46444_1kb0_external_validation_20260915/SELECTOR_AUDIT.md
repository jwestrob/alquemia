# Q46444 / 1KB0 selector audit

Date: 2026-09-15  
Protocol: `pqq_vertical_swap_r2scan3c_native_cpcm_fixed_core_v3`

This audit inspected only the raw 1KB0 coordinate file. It did not inspect or
calculate a v3 energy. The machine-readable authority is `target_spec.tsv`.

## Frozen site

The selected site is model 1, chain A, `A:CA801/CA` with `A:PQQ1800`.
The Ca atom is fully occupied and the PQQ atom graph matches the canonical PDB
CCD naming used by fixed_core_v3. The typed 3.1 A donor ledger is:

| Donor | Distance, A |
|---|---:|
| `A:ASN263/OD1` | 2.378 |
| `A:PQQ1800/O7A` | 2.420 |
| `A:GLU185/OE1` | 2.443 |
| `A:PQQ1800/N6` | 2.500 |
| `A:PQQ1800/O5` | 2.527 |
| `A:GLU185/OE2` | 2.547 |
| `A:ASP308/OD1` | 2.662 |

The fixed protein roles are `A:GLU185`, `A:ASN263`, `A:ASP308`, and
`A:LYS335`. `A:THR310` is the catalytic-Asp+2 homolog; because it is not
acidic it is excluded from the QM core. Its OG1 is 6.002 A from Ca. The closest
catalytic-Asp/cationic-partner contact is `A:ASP308/OD2--A:LYS335/NZ` at
2.861 A. The resulting expected total charges are La -1 and Ca -2.

## Source inventory and exclusions

The raw source contains 1,016 waters, all on chain A. None is within 3.6 A of
Ca; the nearest is `A:HOH1845/O` at 3.850 A. No noncore O/N/S atom lies inside
the 3.1 A coordination gate. The nearest noncore donor-like atom is the
trifluoroborate buffer atom `A:TFB1810/OXT` at 3.338 A. It is excluded with the
rest of TFB1810 and is not replaced.

Raw alternate conformers occur at `A:ALA23`, `A:ALA24`, `A:ASP657`,
`A:GLU198`, `A:HOH2056`, `A:HOH2067`, `A:HOH2079`, `A:HOH2316`,
`A:HOH2351`, `A:LYS328`, and `A:MET598`. None is a frozen core residue.
Preparation uses the immutable residue-consistent highest-occupancy,
blank-then-A policy and fails closed if this inventory or any core selector,
coordinate, atom name, charge, or source hash changes.

The raw-source SHA-256 is
`2b26ccaecd3ed6c89af5e6d8404e39ee09787d79d1904c7c99c48658685aa916`.

