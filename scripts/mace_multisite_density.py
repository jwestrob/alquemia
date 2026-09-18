"""Explicit physical-identity bridge and generating field for real multisite inputs."""
from __future__ import annotations
import argparse
import copy
import math
from pathlib import Path
import shutil
from affordable_common import InvalidArtifact,read_json,record,verify,write_new

PROTOCOL='source_graph_multisite_normalized_responsive_density_inputs_v2'
ALIAS='source_selected_metal_identity_alias_v1'


def aliased(parent):
    p=copy.deepcopy(parent);selected=[a for a in p['physical_atoms'] if a['kind']=='selected_metal']
    if len(selected)!=1 or selected[0]['id']=='metal' or any(a['id']=='metal' for a in p['physical_atoms']):
        raise InvalidArtifact('unique real source selected-metal identity required')
    original=selected[0]['id'];selected[0]['id']='metal'
    p['selected_metal_identity_alias']=dict(protocol=ALIAS,physical_alias='metal',source_id=original,
        coordinate_changes=0,element_changes=0)
    return p


def alias(preparation,output,water_ion_forcefield=None):
    parent=read_json(preparation)
    for pin in [parent['source'],parent['source_preparation'],parent['raw_source'],*parent['implementation'].values()]:verify(pin)
    p=aliased(parent);root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,root/'alias_implementation.py')
    p.update(identity_parent=record(preparation),identity_implementation=record(root/'alias_implementation.py'))
    if water_ion_forcefield:p['generating_water_ion_forcefield']=record(water_ion_forcefield)
    write_new(root/'preparation.json',p);audit_alias(p);return record(root/'preparation.json')


def audit_alias(p):
    parent=read_json(verify(p['identity_parent']));verify(p['identity_implementation'])
    expected=aliased(parent)
    if 'generating_water_ion_forcefield' in p:
        verify(p['generating_water_ion_forcefield']);expected['generating_water_ion_forcefield']=p['generating_water_ion_forcefield']
    if {k:v for k,v in p.items() if k not in ('identity_parent','identity_implementation')}!=expected:
        raise InvalidArtifact('physical identity bridge changed a source state')
    return True


def background(prep):
    """Fixed ff19SB plus the pinned neutral-water/background-Ca templates."""
    from openmm import app,unit,NonbondedForce
    from mace_amoeba_capability import topology
    audit_alias(prep)
    top,positions,ids,bonds,water_bonds,metals=topology(prep,(prep['case_id'],),allow_background_calcium=True)
    ff=app.ForceField(str(verify(prep['forcefield'])),str(verify(prep['generating_water_ion_forcefield'])))
    system=ff.createSystem(top,nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
    nb=next(f for f in system.getForces() if isinstance(f,NonbondedForce))
    q={pid:float(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)) for i,pid in enumerate(ids)}
    physical={a['id']:a for a in prep['physical_atoms']}
    from mace_global_prepare import atom_id
    water_model=app.Modeller(top,positions*unit.angstrom)
    water_model.delete([a for a in water_model.topology.atoms() if physical[atom_id(a)]['kind']!='retained_site_water'])
    water_system=app.ForceField(str(verify(prep['water_forcefield']))).createSystem(water_model.topology,
        nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False)
    water_nb=next(f for f in water_system.getForces() if isinstance(f,NonbondedForce))
    for atom in water_model.topology.atoms():
        actual=float(water_nb.getParticleParameters(atom.index)[0].value_in_unit(unit.elementary_charge))
        if q[atom_id(atom)]!=actual:raise InvalidArtifact('generating water charge differs from original template')
    if set(q)!=(set(physical)-{'metal'}):raise InvalidArtifact('incomplete generating-field inventory')
    groups={}
    for pid in ids:groups.setdefault(pid.rsplit('/',1)[0],[]).append(pid)
    for a in prep['background_metals']:
        if q[a['id']]!=2 or a['formal_charge_e']!=2:raise InvalidArtifact('background Ca template charge differs')
    for w in prep['explicit_waters']:
        names=groups[w['oxygen_id'].rsplit('/',1)[0]]
        if len(names)!=3 or sorted(physical[i]['element'] for i in names)!=['H','H','O'] or abs(math.fsum(q[i] for i in names))>1e-12:
            raise InvalidArtifact('water template not complete/neutral')
    expected=prep['protein_charge_e']+prep['background_metal_charge_e']
    if abs(math.fsum(q.values())-expected)>1e-6:raise InvalidArtifact('generating-field full charge fails')
    return dict(charges_e=q,bonds=[list(b) for b in bonds+water_bonds],forcefield=prep['forcefield'],
        water_ion_forcefield=prep['generating_water_ion_forcefield'],water_reference=prep['water_forcefield'],source=prep['source'],background_metals=prep['background_metals'],
        formal_charge_e=expected,scope='fixed_generating_field_only; not_AMOEBA_environment')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--preparation',required=True);p.add_argument('--output',required=True)
    p.add_argument('--water-ion-forcefield');a=p.parse_args();print(alias(a.preparation,a.output,a.water_ion_forcefield))
