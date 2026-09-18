"""Parameter-only AMOEBA support for declared, source-backed multisite cases."""
from __future__ import annotations
import argparse
from pathlib import Path
import shutil
import time
import numpy as np
from affordable_common import InvalidArtifact,read_json,record,verify,write_new
from mace_amoeba_capability import prepare_case,topology
from mace_tinker_framework_solver import native_call,parse_parameters,check_parameters
PROTOCOL='source_multisite_background_Ca_AMOEBA2018_capability_v1'
NATIVE_PROTOCOL='source_multisite_background_Ca_AMOEBA2018_native_GK_capability_v2'


def prepare(config,output):
    c=read_json(config)
    if c['protocol'] not in (PROTOCOL,NATIVE_PROTOCOL) or not c['cases']:raise InvalidArtifact('declared multisite capability configuration required')
    for pin in [c['plan'],*c['forcefields'],*c['cases'].values()]:verify(pin)
    root=Path(output).resolve();root.mkdir(parents=True,exist_ok=False);impl=root/'implementation';impl.mkdir()
    for n in ('mace_multisite_amoeba.py','mace_amoeba_capability.py','affordable_common.py','mace_global_prepare.py','mace_hybrid.py','mace_tinker_framework_solver.py'):
        shutil.copyfile(Path(__file__).with_name(n),impl/n)
    # The public helper imports additional frozen dependencies from its source
    # installation. Preserve the full upstream snapshot so the CLI can replay.
    parent=read_json(verify(c['implementation_parent']))
    for n,p in parent['implementation'].items():
        if not (impl/n).exists():shutil.copyfile(verify(p),impl/n)
    m=dict(protocol=c['protocol'],config=record(config),plan=c['plan'],cases=c['cases'],forcefields=c['forcefields'],
        implementation={p.name:record(p) for p in impl.glob('*.py')},new_energy_calls=0,new_force_calls=0,new_DFT_calls=0,new_MACE_calls=0)
    native=c['protocol']==NATIVE_PROTOCOL
    if native:
        m['native_source']=c['native_source'];verify(m['native_source'])
    write_new(root/'manifest.json',m);start=time.monotonic();cpu=time.process_time();rows=[]
    for name,pin in c['cases'].items():
        prep=read_json(verify(pin))
        if prep['case_id']!=name:raise InvalidArtifact('source/config identity differs')
        r=prepare_case(pin,[verify(p) for p in (c['forcefields'][:1] if native else c['forcefields'])],root/name,tuple(c['cases']),allow_background_calcium=True,native_gk_only=native)
        if r['status']=='framework_parameterized_no_energy':
            mapping=read_json(verify(r['mapping']));byid={p['id']:p for p in mapping['parameters']}
            backgrounds=[]
            for atom in prep['background_metals']:
                p=byid[atom['id']];multipole=p['multipole_md_units'];gk=p['gk_md_units']
                backgrounds.append(dict(id=atom['id'],charge_e=multipole[0],polarizability_nm3=multipole[-1],GK_parameters_MD=gk))
                if abs(multipole[0]-2)>1e-12 or multipole[-1]<=0 or (gk is not None and gk[1]<=0):raise InvalidArtifact('unsupported calcium parameters')
            r['actual_background_parameters']=backgrounds
            actual={frozenset(b) for b in mapping['protein_bonds']}
            if any(frozenset(b) not in actual for b in prep['source_covalent_connections']):raise InvalidArtifact('actual acetyl/source bond lost')
            r['source_covalent_connections']=prep['source_covalent_connections']
            if native:r['native_GK']=native_parameters(c,prep,mapping,root/name,tuple(c['cases']))
        rows.append(dict(case_id=name,result=record(root/name/'result.json'),status=r['status'],audit=r))
    result=dict(protocol=c['protocol'],manifest=record(root/'manifest.json'),complete=all(r['status']=='framework_parameterized_no_energy' for r in rows),
        rows=rows,wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        new_energy_calls=0,new_force_calls=0,new_DFT_calls=0,new_MACE_calls=0,numerical_score=None,full_model_qualified=False,baseline_changed=False)
    write_new(root/'result.json',result);return result


def native_parameters(config,prep,mapping,root,names):
    from openmm import app
    parent=read_json(verify(config['native_source']));software=read_json(verify(parent['software']))
    top,pos,ids,bonds,water_bonds,metals=topology(prep,names,allow_background_calcium=True)
    if ids!=mapping['system_atom_ids'] or pos.tolist()!=mapping['positions_A']:raise InvalidArtifact('native export physical state changed')
    ff=app.ForceField(str(verify(config['forcefields'][0])));bonded=ff._buildBondedToAtomList(top);types=[]
    for residue in top.residues():
        template,matches=ff._getResidueTemplateMatches(residue,bonded)
        if matches is None:raise InvalidArtifact('native export template unavailable')
        types.extend(int(template.atoms[i].type) for i in matches)
    ix={pid:i+1 for i,pid in enumerate(ids)};neighbors={pid:[] for pid in ids}
    for a,b in bonds+water_bonds:neighbors[a].append(ix[b]);neighbors[b].append(ix[a])
    d=root/'native';d.mkdir();rows=[str(len(ids))+' source multisite parameter-only framework']
    for i,(a,pid,xyz,typ) in enumerate(zip(top.atoms(),ids,pos,types),1):
        rows.append(f'{i} {a.element.symbol} '+' '.join(format(v,'.17g') for v in xyz)+f' {typ} '+' '.join(map(str,sorted(neighbors[pid]))))
    (d/'framework.xyz').write_text('\n'.join(rows)+'\n');(d/'framework.freeze').write_text('0\n')
    opts=['parameters '+str(verify(parent['native_parameters'])),'multipoleterm ONLY','polarizeterm','solvateterm',
          'solvate GK','dielectric 1.0','gk-radius SOLUTE','gkc 2.455','polarization MUTUAL','polar-iter 100','polar-eps 1e-5']
    (d/'framework.key').write_text('\n'.join(opts)+'\n')
    t=dict(xyz=record(d/'framework.xyz'),key=record(d/'framework.key'),mask=record(d/'framework.freeze'),
        frozen_indices=[],poleps=1e-5,rotation=np.eye(3).tolist(),translation_A=[0.,0.,0.])
    receipt=native_call(software['executable'],t,'preflight',d/'native.log',1);write_new(d/'receipt.json',receipt)
    text=verify(receipt['log']).read_text()
    if receipt['returncode'] or 'ALQUEMIA_PREFLIGHT_COMPLETE' not in text:raise InvalidArtifact('native parameter read failed')
    params=parse_parameters(text);check_parameters(params,t,mapping)
    errors=[]
    for a,p in zip(params['atoms'],mapping['parameters']):
        q,alpha=p['multipole_md_units'][0],p['multipole_md_units'][-1]*1000
        errors.append(max(abs(q-a['charge_e']),abs(alpha-a['polarizability_A3'])))
    if max(errors)>1e-10:raise InvalidArtifact('native permanent charge or polarization differs from same XML')
    write_new(d/'parameters.json',params)
    return dict(task=t,parameters=record(d/'parameters.json'),receipt=record(d/'receipt.json'),
        source_parameters=parent['native_parameters'],software=parent['software'],maximum_charge_polarizability_difference=max(errors),
        background_calcium=[dict(physical_id=pid,parameters=params['atoms'][ix[pid]-1]) for pid in ids if pid in {a['id'] for a in prep['background_metals']}],
        new_energy_calls=0,new_response_solves=0)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();r=prepare(a.config,a.output);print({'complete':r['complete'],'cases':[(x['case_id'],x['status']) for x in r['rows']]})
