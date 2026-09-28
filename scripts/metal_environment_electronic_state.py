"""Explicit native electronic states for the separate finite-response research path.

No production parser/default changes. Only one CaII, LaIII or DyIII and closed-shell
ligands are admitted. DyIII uses a physical sextet with explicit 4f electrons;
this is a declared scalar-relativistic hypothesis, not a validated ground state.
"""
from pathlib import Path
import re
import numpy as np
from affordable_common import (InvalidArtifact, BOHR_TO_A, HA_TO_KCAL, energy,
                               record, verify, xyz)
from affordable_response import read_engrad
from mace_omol_vacuum import METHOD, embedded_input, scientific_input

NUMBERS = {'H':1, 'C':6, 'N':7, 'O':8, 'F':9, 'P':15, 'S':16, 'Cl':17,
           'Ca':20, 'La':57, 'Dy':66}
PBE0_PROFILE = 'isolated_pbe0_d4_def2tzvpp_v1'
PBE0_METHOD = 'PBE0 D4 def2-TZVPP def2/J RIJCOSX NoAutostart DefGrid3 TightSCF EnGrad'

STATES = {'Ca': (2, 1, 0), 'La': (3, 1, 46), 'Dy': (3, 6, 28)}


def describe(path, metal, charge, multiplicity):
    if metal not in STATES or type(charge) is not int or type(multiplicity) is not int:
        raise InvalidArtifact('unsupported metal or noninteger electronic state')
    atoms = xyz(path)
    symbols = [a[0] for a in atoms]
    if any(s not in NUMBERS for s in symbols):
        raise InvalidArtifact('unsupported electronic-region element')
    if symbols.count(metal) != 1 or sum(s in STATES for s in symbols) != 1:
        raise InvalidArtifact('single target metal required; spectator metals need a separate explicit model')
    oxidation, expected_mult, core = STATES[metal]
    if multiplicity != expected_mult:
        raise InvalidArtifact('state differs from declared closed-ligand metal hypothesis')
    all_electrons = sum(NUMBERS[s] for s in symbols) - charge
    explicit = all_electrons - core
    unpaired = multiplicity - 1
    if explicit < unpaired or (explicit - unpaired) % 2:
        raise InvalidArtifact('electron count and multiplicity parity differ')
    return dict(metal=metal, oxidation_state_hypothesis=oxidation, charge=charge,
                physical_multiplicity=multiplicity, orca_spin_convention='2S+1',
                alpha_minus_beta=unpaired, all_electron_count=all_electrons,
                ecp_core_electrons=core, explicit_electrons=explicit,
                alpha_electrons=(explicit+unpaired)//2, beta_electrons=(explicit-unpaired)//2,
                ecp_name='def2-ECP' if core else None, f_in_core=False if metal=='Dy' else None,
                spin_orbit_included=False, state_qualification='hypothesis_not_ground_state_validation',
                xyz=record(path))


def input_text(charge, multiplicity, embedded=True, guess=None, method_profile=None):
    if type(charge) is not int or type(multiplicity) is not int or multiplicity not in (1, 6):
        raise InvalidArtifact('unsupported input state')
    if guess not in (None, 'HCore', 'PModel'):
        raise InvalidArtifact('unsupported native initial guess')
    if method_profile not in (None, PBE0_PROFILE) or (method_profile and embedded):
        raise InvalidArtifact('unsupported electronic method/environment combination')
    original = (f'! {PBE0_METHOD}\n* xyzfile {charge} 1 core.xyz\n' if method_profile else
                embedded_input(charge) if embedded else scientific_input(charge))
    if guess is not None:
        first, rest = original.split('\n', 1)
        original = first + f'\n%scf\n Guess {guess}\nend\n' + rest
    return original.replace(f'* xyzfile {charge} 1 core.xyz',
                            f'* xyzfile {charge} {multiplicity} core.xyz')


_FLOAT = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'


def electronic_evidence(text, symbols, multiplicity):
    """Read printed diagnostics; absence is not a zero or a stability pass.

    Multiple SCF records are retained. The last atomic population block must
    map exactly onto the supplied coordinates before atomwise values are used.
    This parser also accepts archived non-target elements for format tests;
    `describe` remains the only elemental capability gate.
    """
    def values(pattern):
        return [float(v.replace('D', 'E').replace('d', 'e'))
                for v in re.findall(pattern, text, re.M)]
    s2 = values(r'^\s*Expectation value of <S\*\*2>\s*:\s*(' + _FLOAT + r')')
    ideal = (multiplicity - 1) * (multiplicity + 1) / 4
    populations = {}
    for scheme in ('MULLIKEN', 'LOEWDIN'):
        blocks = list(re.finditer(r'^' + scheme + r' ATOMIC CHARGES AND SPIN POPULATIONS\s*$', text, re.M))
        if not blocks:
            populations[scheme.lower()] = None
            continue
        tail = text[blocks[-1].end():]
        rows = []
        for line in tail.splitlines():
            if not rows and (not line.strip() or set(line.strip()) == {'-'}):
                continue
            match = re.fullmatch(r'\s*(\d+)\s+([A-Z][a-z]?)\s*:\s*(' + _FLOAT + r')\s+(' + _FLOAT + r')\s*', line)
            if not match:
                break
            idx, element, charge, spin = match.groups()
            rows.append(dict(atom_index=int(idx), element=element,
                             charge=float(charge.replace('D','E')),
                             spin_population=float(spin.replace('D','E'))))
        if [r['atom_index'] for r in rows] != list(range(len(symbols))) or [r['element'] for r in rows] != symbols:
            raise InvalidArtifact('printed atomic spin population mapping differs')
        populations[scheme.lower()] = dict(atoms=rows,
            spin_sum=sum(r['spin_population'] for r in rows),
            charge_sum=sum(r['charge'] for r in rows),
            block_count=len(blocks))
    # Deliberately anchored: the credits contain the words "stability analysis".
    # Retain genuine headings/verdicts, but never infer stability from convergence.
    stability_lines = [line for line in text.splitlines() if re.search(
        r'^\s*(?:SCF STABILITY ANALYSIS|STABILITY ANALYSIS|(?:The\s+)?wavefunction (?:is|was) (?:stable|unstable))\b', line, re.I)]
    stability = 'not_run_or_not_printed' if not stability_lines else 'printed_unqualified'
    alpha = values(r'^\s*N\(Alpha\)\s*:\s*(' + _FLOAT + r')')
    beta = values(r'^\s*N\(Beta\)\s*:\s*(' + _FLOAT + r')')
    return dict(s2_values=s2, s2_last=s2[-1] if s2 else None,
                spin_pure_s2=ideal, s2_deviation=(s2[-1]-ideal) if s2 else None,
                s2_interpretation='UKS diagnostic, not proof of a spin-pure eigenstate',
                population_analyses=populations,
                integrated_alpha_values=alpha, integrated_beta_values=beta,
                stability_status=stability, stability_evidence_lines=stability_lines,
                local_f_occupation_status='not_extracted',
                spin_localization_qualified=False)


def parse(task, output, engrad, *, embedded=True):
    """Parse actual native output; retain missing spin/stability evidence explicitly."""
    state = describe(verify(task['xyz']), task['metal'], task['charge'], task['multiplicity'])
    if verify(task['input']).read_text() != input_text(task['charge'], task['multiplicity'], embedded, task.get('scf_guess'), task.get('method_profile')):
        raise InvalidArtifact('native research input differs')
    text = Path(output).read_text(); value = energy(output)
    if task.get('scf_guess'):
        wanted = 'HCORE' if task['scf_guess']=='HCore' else 'MODEL POTENTIAL'
        if f'INITIAL GUESS: {wanted}' not in text:
            raise InvalidArtifact('executed initial guess differs')
    if not re.search(r'Program Version\s+6\.1\.1\b',text):
        raise InvalidArtifact('ORCA version differs')
    if re.search(r'^\s*(?:CPCM SOLVATION MODEL|SMD SOLVATION(?: MODEL)?|COSMO SOLVATION(?: MODEL)?)\s*$',text,re.M|re.I):
        raise InvalidArtifact('finite dry reference contains undeclared solvent')
    gradient_text = re.sub(r'(?m)^  ===> : Will NOT make the numerical gradients translationally invariant,\n         in case numerical gradients are calculated!\n','',text)
    if re.search(r'numerical (?:gradient|differentiation)',gradient_text,re.I):
        raise InvalidArtifact('numerical gradients unsupported')
    if embedded:
        pc=verify(task['pointcharges']); count=int(pc.read_text().splitlines()[0])
        counts=re.findall(r'Reading point charge file\s+\.{2,}\s+ok\s+\((\d+) point charges\)',text)
        if not counts or any(int(n)!=count for n in counts) or 'environment.pc' not in text:
            raise InvalidArtifact('executed point-charge inventory differs')
    actual_ecp=re.findall(r'Type\s+(\w+)\s+ECP\s+(\S+)\s+\(replacing\s+(\d+)\s+core electrons',text)
    expected_ecp=[(task['metal'],'Def2-ECP',str(state['ecp_core_electrons']))] if state['ecp_core_electrons'] else []
    if actual_ecp!=expected_ecp:
        raise InvalidArtifact('executed native ECP differs')
    for pattern,wanted in [(r'Total Charge\s+Charge\s+\.{2,}\s+(-?\d+)',task['charge']),
                           (r'Multiplicity\s+Mult\s+\.{2,}\s+(\d+)',task['multiplicity']),
                           (r'Number of Electrons\s+NEL\s+\.{2,}\s+(\d+)',state['explicit_electrons'])]:
        if re.findall(pattern,text)!=[str(wanted)]:
            raise InvalidArtifact('executed electronic state differs')
    patterns={'analytic_scf':r'ORCA SCF GRADIENT CALCULATION','native_dispersion':r'DISPERSION GRADIENT',
              'native_gcp':r'gCP correction\s+\.{2,}\s+done','total_cartesian':r'CARTESIAN GRADIENT',
              'native_D4':r'DFTD4','native_gcp_energy':r'gCP correction\s+[-+0-9.]'}
    if task.get('method_profile') == PBE0_PROFILE:
        del patterns['native_gcp']; del patterns['native_gcp_energy']
        if re.search(r'gCP correction\s+[-+0-9.]',text):
            raise InvalidArtifact('undeclared composite gCP in PBE0 reference')
    if state['ecp_core_electrons']:
        patterns['native_ecp']=r'ECP gradient\s+\(SHARK\)\s+\.{2,}\s+done'
    components={k:bool(re.search(p,text,re.I)) for k,p in patterns.items()}
    if not all(components.values()):raise InvalidArtifact('analytic native component absent')
    raw=read_engrad(engrad); atoms=xyz(verify(task['xyz'])); numbers=[NUMBERS[a[0]] for a in atoms]
    if raw['atom_count']!=len(atoms) or abs(raw['energy_Ha']-value)>1e-8 or raw['atomic_numbers'].tolist()!=numbers:
        raise InvalidArtifact('gradient energy or atom inventory differs')
    if not np.allclose(raw['coordinates_bohr']*BOHR_TO_A,[a[1:] for a in atoms],atol=1e-6,rtol=0):
        raise InvalidArtifact('gradient geometry differs')
    hftypes=re.findall(r'Hartree-Fock type\s+HFTyp\s+\.{2,}\s+(\w+)', text)
    expected_hftype='UHF' if task['multiplicity'] > 1 else 'RHF'
    if hftypes != [expected_hftype]:
        raise InvalidArtifact('executed restricted/unrestricted wavefunction differs')
    evidence=electronic_evidence(text, [a[0] for a in atoms], task['multiplicity'])
    if not np.isfinite(raw['gradient_Ha_per_bohr']).all():
        raise InvalidArtifact('nonfinite native gradient')
    return dict(energy_hartree=value, electronic_state=state, native_components=components,
                gradient_kcal_mol_per_A=(raw['gradient_Ha_per_bohr']*HA_TO_KCAL/BOHR_TO_A).tolist(),
                quantity='gradient_not_force', engrad=record(engrad), output=record(output),
                spin_evidence=evidence, executed_hftype=expected_hftype,
                electronic_state_qualified=False, full_hybrid_qualified=False,
                physical_boundary_force_status='not_projected', classification=None)
