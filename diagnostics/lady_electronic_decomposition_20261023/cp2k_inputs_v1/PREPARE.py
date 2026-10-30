#!/usr/bin/env python3
"""Translate exact archived coordinates into four native decks; no calculation."""
import hashlib
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_once(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_text() != text:
            raise RuntimeError(f"Refusing changed artifact: {path}")
    else:
        path.write_text(text)


COMMON = """&GLOBAL
  PROJECT LADY_ALMO
  RUN_TYPE ENERGY
  PRINT_LEVEL MEDIUM
&END GLOBAL
&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME /assets/BASIS_LADY
    POTENTIAL_FILE_NAME /assets/POTENTIAL_LADY
    CHARGE {charge}
    MULTIPLICITY 1
    UKS F
    &MGRID
      CUTOFF [Ry] 1200
      REL_CUTOFF [Ry] 100
      NGRIDS 5
    &END MGRID
    &POISSON
      PERIODIC NONE
      POISSON_SOLVER MT
    &END POISSON
    &QS
      METHOD GAPW
      ALMO_SCF T
      EPS_DEFAULT 1.0E-12
      EPS_PGF_ORB 1.0E-10
      EPSFIT 1.0E-6
      EPSSVD 1.0E-10
      EPSISO 1.0E-14
      EPSRHO0 1.0E-10
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-10
      MAX_SCF 300
      IGNORE_CONVERGENCE_FAILURE F
      &OT
        MINIMIZER CG
        PRECONDITIONER FULL_ALL
      &END OT
    &END SCF
    &ALMO_SCF
      ALMO_SCF_GUESS MOLECULAR
      ALMO_EXTRAPOLATION_ORDER 0
      XALMO_EXTRAPOLATION_ORDER 0
      ALMO_ALGORITHM DIAG
      XALMO_ALGORITHM PCG
      DELOCALIZE_METHOD FULL_SCF
      MO_OVERLAP_INV_ALG DENSE_CHOLESKY
      EPS_FILTER 1.0E-12
      CONSTRUCT_NLMOS F
      &ALMO_OPTIMIZER_DIIS
        EPS_ERROR 1.0E-9
        MAX_ITER 300
        N_DIIS 8
        EPS_ERROR_EARLY -1
        MAX_ITER_EARLY -1
      &END ALMO_OPTIMIZER_DIIS
      &XALMO_OPTIMIZER_PCG
        EPS_ERROR 1.0E-9
        MAX_ITER 300
        MAX_ITER_OUTER_LOOP 0
        CONJUGATOR HAGER_ZHANG
        PRECONDITIONER FULL
        LIN_SEARCH_EPS_ERROR 1.0E-3
        LIN_SEARCH_STEP_SIZE_GUESS 0.1
        EPS_ERROR_EARLY -1
        MAX_ITER_EARLY -1
      &END XALMO_OPTIMIZER_PCG
      &ANALYSIS F
      &END ANALYSIS
      &FRAGMENT 1
        CHARGE 3
        MULTIPLICITY 1
      &END FRAGMENT
      &FRAGMENT 2
        CHARGE {shell_charge}
        MULTIPLICITY 1
      &END FRAGMENT
    &END ALMO_SCF
    &XC
      &XC_FUNCTIONAL NO_SHORTCUT
        &PBE
          SCALE_X 0.75
          SCALE_C 1.0
        &END PBE
      &END XC_FUNCTIONAL
      &HF
        FRACTION 0.25
        &INTERACTION_POTENTIAL
          POTENTIAL_TYPE COULOMB
        &END INTERACTION_POTENTIAL
        &SCREENING
          EPS_SCHWARZ 1.0E-12
          SCREEN_ON_INITIAL_P F
        &END SCREENING
        &MEMORY
          MAX_MEMORY 512
        &END MEMORY
      &END HF
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC [angstrom] 30.0 30.0 30.0
      PERIODIC NONE
    &END CELL
    &TOPOLOGY
      CONNECTIVITY USER
      &CENTER_COORDINATES F
      &END CENTER_COORDINATES
    &END TOPOLOGY
    &COORD
{coordinates}
    &END COORD
{kinds}
  &END SUBSYS
&END FORCE_EVAL
"""


def main():
    ledger = json.loads((BASE/'ENERGY_LEDGER.json').read_text())
    records = {}
    for gid, system in ledger['systems'].items():
        xyz = Path(system['provenance']['xyz']['path'])
        assert digest(xyz) == system['provenance']['xyz']['sha256']
        lines = xyz.read_text().splitlines()
        atom_lines = lines[2:]
        assert len(atom_lines) == int(lines[0]) == system['geometry_state']['atoms']
        assert all(len(row.split()) == 4 for row in atom_lines)
        metal = system['geometry_state']['metal']
        assert atom_lines[0].split()[0] == metal
        coordinates = '\n'.join('      '+row+(' METAL 1' if i == 0 else ' SHELL 2')
                                for i, row in enumerate(atom_lines))
        kinds = []
        for element in dict.fromkeys(row.split()[0] for row in atom_lines):
            potential = ('ECP-LADY-La46' if metal == 'La' else 'ECP-LADY-Dy55') if element == metal else 'ALL'
            kinds.append(f'''    &KIND {element}
      ELEMENT {element}
      BASIS_SET LADY-ARCHIVED
      POTENTIAL {potential}
      RADIAL_GRID 200
      LEBEDEV_GRID 590
'''+('      ECP_SEMI_LOCAL TRUE\n' if element == metal else '')+'    &END KIND')
        charge = system['geometry_state']['charge']
        deck = COMMON.format(charge=charge, shell_charge=charge-3,
                             coordinates=coordinates, kinds='\n'.join(kinds))
        path = BASE/'decks'/f'{gid}.inp'
        write_once(path, deck)
        # Recover every original coordinate token from the completed native deck.
        written = deck.split('    &COORD\n')[1].split('    &END COORD')[0].splitlines()
        assert [row.split()[:4] for row in written] == [row.split() for row in atom_lines]
        records[gid] = {'deck': str(path), 'deck_sha256': digest(path),
                        'xyz': str(xyz), 'xyz_sha256': digest(xyz),
                        'atom_count': len(atom_lines), 'fragment_count': 2,
                        'geometry_state': system['geometry_state'],
                        'coordinate_tokens_unchanged': True}
    for motif in ['aquo8_central', 'hydrated_acetate_mono_central']:
        a = Path(records[motif+'_La']['xyz']).read_text().splitlines()
        b = Path(records[motif+'_Dy']['xyz']).read_text().splitlines()
        assert a[3:] == b[3:]
        assert a[2].split()[1:] == b[2].split()[1:]
    result = {'status': 'four_decks_prepared_no_molecular_calls',
              'generator_sha256': digest(Path(__file__)),
              'ledger_sha256': digest(BASE/'ENERGY_LEDGER.json'),
              'assets': {n: digest(BASE/'assets'/n) for n in ['BASIS_LADY','POTENTIAL_LADY','PROVENANCE.json']},
              'systems': records, 'remaining_release': 'blocked_until_pilot_admission',
              'molecular_calls': 0}
    write_once(BASE/'DECK_PROVENANCE.json', json.dumps(result, indent=2)+'\n')
    print('Prepared and checked four decks; no molecular calls.')


if __name__ == '__main__':
    main()
