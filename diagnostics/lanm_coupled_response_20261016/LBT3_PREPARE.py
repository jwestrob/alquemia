#!/usr/bin/env python3
"""Two consumed LBT3 model-1 whole peptides; preparation only, Slurm required.

OpenMM adds/minimizes two missing terminal H per source. No MACE/QM scoring.
All deposited heavy atoms remain fixed; no solvent or fragment caps are added.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import random
import sys
import time

PROTOCOL = 'lbt3_model1_whole17_pH6_ff19sb_terminal_completion_LaDy_v1'
SEQUENCE = 'FIDTNNDGWIEGDELLA'
RESNAMES = 'PHE ILE ASP THR ASN ASN ASP GLY TRP ILE GLU GLY ASP GLU LEU LEU ALA'.split()
SEED = 20260918


def pdb_atoms(path):
    """Read only model 1, retaining literal deposited identities/coordinates."""
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('MODEL') and int(line[10:14]) != 1:
            break
        if line.startswith('ENDMDL'):
            break
        if not line.startswith(('ATOM  ', 'HETATM')):
            continue
        if line[16].strip() or line[26].strip():
            raise ValueError('alternate/insertion not supported in this fixed source')
        rows.append({'id':f'{line[21]}/{int(line[22:26])}//{line[12:16].strip()}',
                     'resname':line[17:20].strip(), 'element':line[76:78].strip().title(),
                     'xyz_A':[float(line[a:b]) for a,b in ((30,38),(38,46),(46,54))]})
    return rows


def prepare_source(name, pins, manifest_row, ffpath, legacy, out, work):
    import numpy as np
    import openmm as mm
    from openmm import app, unit
    from affordable_common import InvalidArtifact, record, verify, write_new
    from mace_global_prepare import FF, atom_id, terminal_position
    from mace_hybrid import write_xyz

    raw, model1 = verify(pins['raw']), verify(pins['model1'])
    if (pins['raw']['sha256'] != manifest_row['source_sha256'] or
            pins['model1']['sha256'] != manifest_row['first_model']['sha256'] or
            manifest_row['model'] != 1 or record(ffpath) != record(FF)):
        raise InvalidArtifact('source/forcefield differs from pinned established recipe')
    source = pdb_atoms(raw)
    if source != pdb_atoms(model1) or len(source) != 253:
        raise InvalidArtifact('verbatim model-1 inventory differs')
    original_metal = 'La' if name == '7CCO' else 'Lu'
    target = source[-1]
    if (len({a['id'] for a in source}) != 253 or
            target['id'] != f'A/101//{original_metal.upper()}' or
            target['element'] != original_metal or
            any(a['element'] not in ('H','C','N','O') for a in source[:-1])):
        raise InvalidArtifact('unexpected source metal, waters, spectators or duplicate identity')
    pdb = app.PDBFile(str(model1)); model = app.Modeller(pdb.topology, pdb.positions)
    model.delete([a for a in model.topology.atoms() if a.residue.name == original_metal.upper()])
    residues = list(model.topology.residues())
    if ([(r.chain.id,r.id,r.name) for r in residues] !=
            [('A',str(i+1),r) for i,r in enumerate(RESNAMES)]):
        raise InvalidArtifact('whole 17-residue construct differs')
    before = {atom_id(a):(a.element.symbol, np.array(model.positions[a.index].value_in_unit(unit.angstrom)))
              for a in model.topology.atoms()}
    source_map = {}
    for aid,(element,position) in before.items():
        matches = [s for s in source[:-1] if s['element']==element and
                   np.allclose(s['xyz_A'], position, rtol=0, atol=1e-8)]
        if len(matches) != 1 or aid.split('/')[:3] != matches[0]['id'].split('/')[:3]:
            raise InvalidArtifact('source-to-native atom mapping is ambiguous')
        source_map[aid] = matches[0]
    if len(source_map) != 252 or len({r['id'] for r in source_map.values()}) != 252:
        raise InvalidArtifact('incomplete deposited-atom mapping')
    first = {a.name for a in residues[0].atoms()}
    last = {a.name:a for a in residues[-1].atoms()}
    if 'H' not in first or {'H2','H3'} & first or 'OXT' in last:
        raise InvalidArtifact('terminal inventory differs from fixed completion policy')
    carbon = last['C']
    neighbors = [b if a==carbon else a for a,b in model.topology.bonds() if carbon in (a,b)]
    if {a.name for a in neighbors} != {'CA','O'} or any(a.residue != residues[-1] for a in neighbors):
        raise InvalidArtifact('terminal completion would change existing covalent connectivity')
    position = np.array(model.positions.value_in_unit(unit.angstrom))
    oxt = terminal_position(*(position[last[k].index] for k in ('C','CA','O')))
    oxygen = model.topology.addAtom('OXT',app.element.oxygen,residues[-1])
    model.topology.addBond(carbon,oxygen)
    model.positions = np.concatenate((position,oxt[None,:]))*.1*unit.nanometer
    ff = app.ForceField(str(ffpath))
    variants = [r.name if r.name in ('ASP','GLU') else None for r in residues]
    random.seed(SEED); np.random.seed(SEED)
    # Native API fixes every existing atom and minimizes only newly added H.
    work['native_missing_H_calls_attempted'] += 1
    chosen = model.addHydrogens(ff,pH=6.0,variants=variants,
                               platform=mm.Platform.getPlatformByName('Reference'))
    work['native_missing_H_calls_completed'] += 1
    after = {atom_id(a):(a.element.symbol,np.array(model.positions[a.index].value_in_unit(unit.angstrom)))
             for a in model.topology.atoms()}
    additions = set(after)-set(before)
    if additions != {'A/1//H2','A/1//H3','A/17//OXT'} or set(before)-set(after):
        raise InvalidArtifact('completion altered more than the three declared terminal atoms')
    if any(after[k][0] != e or not np.allclose(after[k][1],p,rtol=0,atol=1e-12)
           for k,(e,p) in before.items()):
        raise InvalidArtifact('native missing-H preparation moved a deposited atom')
    if not np.allclose(after['A/17//OXT'][1],oxt,rtol=0,atol=1e-12):
        raise InvalidArtifact('native missing-H preparation moved completed OXT')
    out.mkdir(parents=True,exist_ok=False)
    completed = out/'completed_peptide.pdb'
    with completed.open('x') as handle:
        app.PDBFile.writeFile(model.topology,model.positions,handle,keepIds=True)
    # Native export operations match the established full-scaffold exporter.
    # Do not call protein(): that helper always projects all source H bonds.
    native = app.PDBFile(str(completed))
    system = ff.createSystem(native.topology,nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
    atoms = [{'id':atom_id(a),'element':a.element.symbol,
              'xyz_A':list(native.positions[a.index].value_in_unit(unit.angstrom)),
              'resname':a.residue.name,'kind':'protein_source' if atom_id(a) in source_map else 'terminal_completion'}
             for a in native.topology.atoms()]
    if len(atoms) != 255 or len({a['id'] for a in atoms}) != 255:
        raise InvalidArtifact('completed peptide inventory differs from 255')
    for a in atoms:
        aid = a['id']
        if aid in source_map:
            a['deposited_source'] = source_map[aid]
            if not np.allclose(a['xyz_A'],source_map[aid]['xyz_A'],rtol=0,atol=1e-12):
                raise InvalidArtifact('deposited H or heavy atom changed after PDB serialization')
        else:
            a['terminal_completion'] = True
    # Native full-peptide system: no metal force-field parameter is invented.
    ids = [atom_id(a) for a in native.topology.atoms()]
    if ids != [a['id'] for a in atoms]:
        raise InvalidArtifact('native Amber and physical atom order differ')
    nb = next(f for f in system.getForces() if isinstance(f,mm.NonbondedForce))
    native_atoms = [dict(a,charge_e=nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge))
                    for i,a in enumerate(atoms)]
    if abs(sum(a['charge_e'] for a in native_atoms)+5) > 1e-6:
        raise InvalidArtifact('native atom charges do not sum to -5')
    bondforce = next(f for f in system.getForces() if isinstance(f,mm.HarmonicBondForce))
    bonds = []
    for i in range(bondforce.getNumBonds()):
        a,b,length,_ = bondforce.getBondParameters(i)
        bonds.append({'atom_a_id':ids[a],'atom_b_id':ids[b],
                      'native_particle_indices':[a,b],
                      'distance_A':float(np.linalg.norm(np.array(atoms[a]['xyz_A'])-atoms[b]['xyz_A'])),
                      'forcefield_equilibrium_A':length.value_in_unit(unit.angstrom),
                      'involves_H':atoms[a]['element']=='H' or atoms[b]['element']=='H'})
    peptide_bonds = [b for b in bonds if b['atom_a_id'].split('/')[1] != b['atom_b_id'].split('/')[1]]
    if len(peptide_bonds) != 16:
        raise InvalidArtifact('incomplete full-peptide covalent graph')
    if any(not 1.0 <= b['distance_A'] <= 1.8 for b in peptide_bonds):
        raise InvalidArtifact('existing full-chain peptide-connection gate failed')
    hydrogen_bonds = [b for b in bonds if b['involves_H']]
    audit = {'all_deposited_H_preserved':True,'all_deposited_heavy_preserved':True,
             'H_bond_lengths':hydrogen_bonds,'H_bond_length_extrema_A':
             [min(b['distance_A'] for b in hydrogen_bonds),max(b['distance_A'] for b in hydrogen_bonds)],
             'native_equilibrium_deviations_are_diagnostics_not_repair_triggers':True,
             'source_H_normalization':False}
    native_xml = out/'peptide_system.xml'
    native_xml.write_text(mm.XmlSerializer.serialize(system))
    write_new(out/'PEPTIDE_ATOMS.json',native_atoms)
    write_new(out/'BONDS.json',bonds)
    write_new(out/'HYDROGEN_AUDIT.json',audit)
    physical = atoms+[{'id':'metal','element':'M','xyz_A':target['xyz_A'],
                      'kind':'selected_metal','deposited_source':target}]
    distances = np.linalg.norm(np.array([a['xyz_A'] for a in physical])[:,None,:]-
                               np.array([a['xyz_A'] for a in physical])[None,:,:],axis=-1)
    np.fill_diagonal(distances,np.inf)
    minimum = float(np.min(distances))
    if minimum < .5:
        raise InvalidArtifact('existing whole-chain absolute overlap gate below 0.5 A')
    write_new(out/'PHYSICAL.json',physical)
    # The exact legacy cap builder consumes these prepared physical atoms,
    # mapped back to deposited names. No old XYZ is used as a new score input.
    from dataclasses import replace
    old_atoms = legacy.parse_first_model(model1)[2]
    by_source = {a['deposited_source']['id']:a for a in physical if 'deposited_source' in a}
    prepared_legacy = [replace(a,xyz=tuple(by_source[f'{a.chain}/{a.resnum}//{a.name}']['xyz_A']))
                       for a in old_atoms]
    endpoint = next(e for e in legacy.ENDPOINTS if e['pdb_id']==name)
    core,fragments,contacts = legacy.build_cluster(prepared_legacy,endpoint)
    if len(core) != 48 or sum(f['formal_charge'] for f in fragments) != -4:
        raise InvalidArtifact('legacy matched-core inventory or formal charge differs')
    core_map = [{'core_index':0,'kind':'source','physical_index':255,'physical_id':'metal'}]
    by_id = {a['id']:i for i,a in enumerate(physical)}
    offset = 1
    for f in fragments:
        for atom_name in f['source_atom_names']:
            source_id = f"A/{f['resnum']}//{atom_name}"
            a = by_source[source_id]; c = core[offset]
            if c.element != a['element'] or list(c.xyz) != a['xyz_A']:
                raise InvalidArtifact('shared real atoms differ between whole peptide and matched core')
            core_map.append({'core_index':offset,'kind':'source','physical_index':by_id[a['id']],
                             'physical_id':a['id'],'deposited_source_id':source_id})
            offset += 1
        for cap in f['generated_caps']:
            core_map.append({'core_index':offset,'kind':'synthetic_cap','definition':cap,
                             'xyz_A':list(core[offset].xyz),
                             'whole_physical_atom':False})
            offset += 1
    if offset != 48:
        raise InvalidArtifact('matched-core source/cap map incomplete')
    write_new(out/'MATCHED_CORE_MAP.json',{'mapping':core_map,'fragments':fragments,'contacts':contacts,
                                         'legacy_builder':record(legacy.__file__),
                                         'shared_real_atoms_exactly_equal':True})
    core_physical = [{'id':c.label,'element':'M' if i==0 else c.element,'xyz_A':list(c.xyz)}
                     for i,c in enumerate(core)]
    write_new(out/'MATCHED_CORE_PHYSICAL.json',core_physical)
    endpoints = {}; tasks = []; atomic_numbers = {'H':1,'C':6,'N':7,'O':8,'La':57,'Dy':66}
    for representation,particle_rows,q,metal_index in (('whole',physical,-2,255),('matched_core',core_physical,-1,0)):
        endpoints[representation] = {}
        for metal,mult in (('La',1),('Dy',6)):
            coords = [(metal if a['element']=='M' else a['element'],*a['xyz_A']) for a in particle_rows]
            electrons = sum(atomic_numbers[a[0]] for a in coords)-q
            if (electrons-(mult-1))%2:
                raise InvalidArtifact('physical electron/spin parity differs')
            path = out/f'{representation}_{metal}.xyz'; write_xyz(path,coords)
            task = {'task_id':f'{name}_{representation}_{metal}','source_id':name,
                    'representation':representation,'metal':metal,'xyz':record(path),'charge':q,
                    'multiplicity':mult,'spin_multiplicity':mult,'physical_spin_2S':mult-1,
                    'atom_count':len(coords),'all_electron_count':electrons,'oxidation_state':3,
                    'metal_index':metal_index,'physical_source':record(out/'PHYSICAL.json'),
                    'source_cap_map':record(out/'MATCHED_CORE_MAP.json') if representation=='matched_core' else None,
                    'electronic_state':'LaIII closed shell' if metal=='La' else
                        'DyIII high-spin sextet hypothesis; no spin-orbit qualification'}
            endpoints[representation][metal] = task; tasks.append(task)
    result = {'protocol_id':PROTOCOL,'source_id':name,'status':'prepared','sequence':SEQUENCE,
              'source_pins':pins,'source_conditioned_metal':original_metal,'model':1,
              'physical_atoms':record(out/'PHYSICAL.json'),'atoms':256,'protein_atoms':255,
              'protein_charge':-5,'charge':-2,'endpoints':endpoints,'tasks':tasks,
              'waters':[],'spectators':[],'caps':[],'completed_peptide':record(completed),
              'source_to_native':source_map,'native_peptide_system':record(native_xml),
              'native_peptide_atoms':record(out/'PEPTIDE_ATOMS.json'),
              'covalent_bonds':record(out/'BONDS.json'),'hydrogen_audit':record(out/'HYDROGEN_AUDIT.json'),
              'terminal_additions':[{'id':k,'element':after[k][0],'xyz_before_PDB_serialization_A':after[k][1].tolist()}
                                    for k in sorted(additions)],
              'terminal_chemistry':'free N-terminal NH3+ / C-terminal COO-; Asp3/7/13 and Glu11/14 deprotonated',
              'hydrogen_preparation':{'pH':6.0,'variants':chosen,'seed':SEED,'platform':'Reference',
                  'new_H_only_native_minimization':True,'native_max_iterations':50,
                  'actual_iterations':'unknown_API_returns_no_iteration_count','source_H_normalization':False},
              'deposited_heavy_coordinates_preserved':True,'deposited_H_coordinates_preserved':True,
              'minimum_all_atom_distance_A':minimum,
              'peptide_bond_count':16,'target_has_no_covalent_graph_edges':True,
              'forcefield':record(ffpath),'scoring_calls':0,'quantum_calls':0,
              'parent_mechanics_scope':'complete native peptide only; target interaction belongs to eventual whole-system scorer',
              'interpretation':'consumed opposing-preference peptide control; fixed dry source states, not solution affinity or LanM success'}
    from collections import Counter
    result['source_element_counts'] = dict(Counter(a['element'] for a in source))
    result['prepared_peptide_element_counts'] = dict(Counter(a['element'] for a in atoms))
    write_new(out/'PREPARATION.json',result)
    return {'source_id':name,'status':'prepared','preparation':record(out/'PREPARATION.json'),'tasks':tasks}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('repo-root','config','output'):
        parser.add_argument('--'+key,required=True,type=Path)
    args=parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('all molecular preparation, including missing-H placement, requires Slurm')
    sys.path.insert(0,str(args.repo_root.resolve()/'scripts'))
    from affordable_common import InvalidArtifact, read_json, record, verify, write_new
    cfg=read_json(args.config); verify(cfg['plan'])
    manifest=read_json(verify(cfg['benchmark_manifest'])); ffpath=verify(cfg['forcefield'])
    if manifest['sequence'] != SEQUENCE or set(cfg['sources']) != {'7CCO','7CCN'}:
        raise InvalidArtifact('exact two-source LBT3 scope required')
    rows={r['pdb_id']:r for r in manifest['endpoints']}
    sys.path.insert(0,str(verify(cfg['benchmark_manifest']).parent))
    import build_benchmark as legacy
    out=args.output.resolve(); out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic(); results=[]
    work={'native_missing_H_calls_attempted':0,'native_missing_H_calls_completed':0,
          'actual_missing_H_iterations':'unknown_API_returns_no_iteration_count'}
    for name in ('7CCO','7CCN'):
        try:
            results.append(prepare_source(name,cfg['sources'][name],rows[name],ffpath,legacy,out/name,work))
        except Exception as exc:
            results.append({'source_id':name,'status':'unavailable','reason':str(exc),'exception':type(exc).__name__})
        print(json.dumps(results[-1]),flush=True)
    result={'protocol_id':PROTOCOL,'config':record(args.config),'plan':cfg['plan'],
            'implementation':record(__file__),'rows':results,
            'status':'complete' if all(r['status']=='prepared' for r in results) else 'incomplete',
            'elapsed_s':time.monotonic()-start,'slurm_job_id':os.environ['SLURM_JOB_ID'],
            'software':{p:importlib.metadata.version(p) for p in ('openmm','numpy')},
            'dependencies':{p:record(args.repo_root/'scripts'/p) for p in
                ('mace_global_prepare.py','mace_hybrid.py','affordable_common.py')},
            'legacy_builder':record(legacy.__file__),
            'preparation_work':work,
            'quantum_calls':0,'MACE_calls':0,'prepared_endpoint_cells':4*sum(r['status']=='prepared' for r in results),
            'native_missing_H_minimizations_per_successful_source':1,'integration_steps':0}
    write_new(out/'RESULT.json',result)
    if result['status'] != 'complete':
        raise SystemExit(1)


if __name__=='__main__':
    main()
