"""Balanced, fixed-scaffold water-to-carboxylate proton-location paths."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
from affordable_common import HA_TO_KCAL, InvalidArtifact, cache_key, read_json, record, verify, write_new, xyz
from hydration_square import endpoint
from hydration_water_motion import METHOD, checked_gradient
from mace_hybrid import check_atoms, write_xyz

PROTOCOL='native_r2scan3c_fixed_scaffold_internal_proton_transfer_v1'
SETTINGS={'contact_cutoff_A':2.60,'acid_OH_A':.98,'acid_angle_deg':[90.,130.],
          'fractions':[.25,.50,.75,1.],'other_H_min_A':.70,'other_heavy_min_A':.75,
          'coordinate':'donor_centered_shortest_spherical_arc_linear_radius',
          'H_selection':'nearest_acceptor_then_original_index'}


def angle(u,v):
    return float(np.degrees(np.arccos(np.clip(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v),-1,1))))


def path_position(atoms,contact,fraction):
    """One real proton moves continuously; return coordinate and analytic dx/dt."""
    if not 0<=fraction<=1:raise InvalidArtifact('path fraction outside frozen segment')
    q=np.array([a[1:] for a in atoms]);od=q[contact['donor_oxygen_index']]
    oa=q[contact['acceptor_index']];h=q[contact['hydrogen_index']]
    target=oa+SETTINGS['acid_OH_A']*(od-oa)/np.linalg.norm(od-oa)
    v0=h-od;v1=target-od;r0=np.linalg.norm(v0);r1=np.linalg.norm(v1)
    u0=v0/r0;u1=v1/r1;theta=np.arccos(np.clip(np.dot(u0,u1),-1,1));t=fraction
    if theta>np.pi-1e-8:raise InvalidArtifact('antipodal path has no unique shortest arc')
    if theta<1e-9:u=u0;du=np.zeros(3)
    else:
        u=(np.sin((1-t)*theta)*u0+np.sin(t*theta)*u1)/np.sin(theta)
        du=theta*(-np.cos((1-t)*theta)*u0+np.cos(t*theta)*u1)/np.sin(theta)
    radius=(1-t)*r0+t*r1
    return od+radius*u,(r1-r0)*u+radius*du


def path_rows(atoms,contact,fraction):
    coord,_=path_position(atoms,contact,fraction);rows=list(atoms)
    rows[contact['hydrogen_index']]=('H',*coord)
    return rows


def geometry_gate(atoms,contact,fraction):
    rows=path_rows(atoms,contact,fraction);q=np.array([a[1:] for a in rows]);h=contact['hydrogen_index']
    excluded={h,contact['donor_oxygen_index'],contact['acceptor_index']}
    for i,a in enumerate(rows):
        if i in excluded:continue
        cutoff=SETTINGS['other_H_min_A'] if a[0]=='H' else SETTINGS['other_heavy_min_A']
        if np.linalg.norm(q[h]-q[i])<cutoff:raise InvalidArtifact(f'proton path overlaps atom {i} at fraction {fraction}')
    return rows


def source_centers(collection):
    c=read_json(collection);m=read_json(verify(c['manifest']))
    if c['status']!='complete':raise InvalidArtifact('complete real occupancy collection required')
    preparations={p['case']:p for p in (read_json(verify(pin)) for pin in m['preparations'])}
    rows={r['task_id']:r for r in c['rows']};centers={}
    for t in m['reused']:
        p=preparations[t['case']]
        if t['pattern']!='1'*sum(w['role']=='variable' for w in p['water_groups']):
            raise InvalidArtifact('full-water source required')
        r=rows[t['task_id']];actual=endpoint(r['result']['output'],r['result']['receipt'],t['xyz'],t['input'])
        if actual!=r['result'] or verify(t['input']).read_text()!=METHOD+f"\n* xyzfile {t['charge']} 1 core.xyz\n":
            raise InvalidArtifact('source endpoint chemistry/energy changed')
        atoms=xyz(verify(t['xyz']));checked_gradient(r['gradient'],actual,atoms)
        centers[t['task_id']]={'case':t['case'],'metal':t['metal'],'charge':t['charge'],
            'xyz':t['xyz'],'input':t['input'],'result':actual,'gradient':r['gradient'],
            'preparation':next(pin for pin in m['preparations'] if read_json(verify(pin))['case']==t['case'])}
    return m,centers


def contacts(atoms,preparation):
    """Select by source O--O geometry only, before any proton-state energy."""
    q=np.array([a[1:] for a in atoms]);mapping=preparation['atom_graph']['source_to_qm'];found=[]
    for water in preparation['water_groups']:
        if water['role']!='variable':continue
        for atom in mapping:
            if atom['kind']!='source':continue
            src=atom['source']
            if src['canonical_resname'] not in ('ASP','GLU') or src['atom'] not in ('OD1','OD2','OE1','OE2'):continue
            d=float(np.linalg.norm(q[water['oxygen_index']]-q[atom['qm_index']]))
            if d>SETTINGS['contact_cutoff_A']:continue
            same=lambda a:all(a['source'][k]==src[k] for k in ('chain','resnum','insertion_code'))
            carbon_name='CG' if src['canonical_resname']=='ASP' else 'CD'
            carbon=[a['qm_index'] for a in mapping if a['kind']=='source' and same(a) and a['source']['atom']==carbon_name]
            if len(carbon)!=1:raise InvalidArtifact('carboxyl parent carbon unavailable or ambiguous')
            h=min(water['hydrogen_indices'],key=lambda i:(np.linalg.norm(q[i]-q[atom['qm_index']]),i))
            label=f"w{water['source']['chain']}{water['source']['resnum']}_to_{src['chain']}{src['resnum']}{src['atom']}"
            found.append({'contact_id':label,'donor_oxygen_index':water['oxygen_index'],'hydrogen_index':h,
                'acceptor_index':atom['qm_index'],'carboxyl_carbon_index':carbon[0],
                'water_source':water['source'],'acceptor_source':src,'OO_distance_A':d,
                'acid_COH_angle_deg':angle(q[carbon[0]]-q[atom['qm_index']],q[water['oxygen_index']]-q[atom['qm_index']])})
    return sorted(found,key=lambda a:a['contact_id'])


def prepare(collection,agreement,output):
    parent,centers=source_centers(collection);out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    impl=out/'implementation'
    # Reuse an existing complete runner closure; add only this state's analyzer.
    closure=Path(verify(parent['full_collection'])).parent/'implementation'
    shutil.copytree(closure,impl)
    for name in ('chemical_states_proton.py','hydration_water_motion.py','hydration_square.py','affordable_common.py','affordable_response.py','mace_hybrid.py'):
        shutil.copyfile(Path(__file__).with_name(name),impl/name)
    tasks=[];branches=[];definitions={}
    for identity,c in centers.items():
        atoms=xyz(verify(c['xyz']));p=read_json(verify(c['preparation']));selected=contacts(atoms,p)
        for contact in selected:
            key=identity+'__'+contact['contact_id'];definitions[key]={'center_id':identity,'contact':contact}
            for fraction in SETTINGS['fractions']:
                tid=key+'__q'+f'{round(fraction*100):03d}'
                branch={'task_id':tid,'path_id':key,'fraction':fraction,'status':'prepared'}
                try:
                    if not SETTINGS['acid_angle_deg'][0]<=contact['acid_COH_angle_deg']<=SETTINGS['acid_angle_deg'][1]:
                        raise InvalidArtifact('unsupported acid O-H geometry')
                    rows=geometry_gate(atoms,contact,fraction)
                except InvalidArtifact as exc:
                    branch.update(status='unsupported',reason=str(exc));branches.append(branch);continue
                td=out/'tasks'/tid;td.mkdir(parents=True);xp=td/'core.xyz';write_xyz(xp,rows)
                ip=td/'endpoint.inp';ip.write_text(METHOD+f"\n* xyzfile {c['charge']} 1 core.xyz\n")
                task={**branch,'case':c['case'],'metal':c['metal'],'charge':c['charge'],'multiplicity':1,
                    'input':record(ip),'xyz':record(xp),'output_path':str(td/'endpoint.out'),'engrad_path':str(td/'endpoint.engrad'),
                    'state':check_atoms(rows,c['charge'])}
                task['cache_key']=cache_key(task|{'protocol_id':PROTOCOL,'settings':SETTINGS,'source':record(collection)})
                tasks.append(task);branches.append(branch)
    prep={'protocol_id':PROTOCOL,'source_collection':record(collection),'agreement':record(agreement),
        'settings':SETTINGS,'centers':centers,'paths':definitions,'branches':branches,
        'reaction':{'reactants':{'carboxylate':-1,'water':0},'products':{'carboxylic_acid':0,'hydroxide':-1},
                    'delta_H_atoms':0,'delta_total_charge':0,'proton_reservoir_coefficient':0,'water_reservoir_coefficient':0},
        'occupancy_probabilities':None,'pKa':None,'baseline_changed':False}
    write_new(out/'preparation.json',prep)
    manifest={'protocol_id':PROTOCOL,'tasks':tasks,'preparation':record(out/'preparation.json'),
        'agreement':record(agreement),'orca':parent['orca'],'execution_policy':parent['execution_policy'],
        'implementation':{p.name:record(p) for p in impl.glob('*.py')},'execution_resources':{'mpi_ranks':16,'concurrent_tasks':4},
        'expected_new_endpoints':len(branches),'unsupported_endpoints':len(branches)-len(tasks),'baseline_changed':False}
    write_new(out/'manifest.json',manifest);return validate(out/'manifest.json')


def validate(manifest):
    m=read_json(manifest);p=read_json(verify(m['preparation']))
    if p['settings']!=SETTINGS or m['protocol_id']!=PROTOCOL:raise InvalidArtifact('state/path policy changed')
    verify(m['agreement'])
    for pin in m['implementation'].values():verify(pin)
    for t in m['tasks']:
        path=p['paths'][t['path_id']];c=p['centers'][path['center_id']];source=xyz(verify(c['xyz']))
        actual=xyz(verify(t['xyz']));wanted=geometry_gate(source,path['contact'],t['fraction']);h=path['contact']['hydrogen_index']
        if (len(actual)!=len(source) or any(actual[i]!=source[i] for i in range(len(source)) if i!=h)
                or not np.allclose(actual[h][1:],wanted[h][1:],atol=1e-12,rtol=0)
                or t['charge']!=c['charge'] or check_atoms(actual,t['charge'])!=t['state']):raise InvalidArtifact('state atom/charge/coordinate mismatch')
        if verify(t['input']).read_text()!=METHOD+f"\n* xyzfile {c['charge']} 1 core.xyz\n":raise InvalidArtifact('state Hamiltonian changed')
    verify(m['orca'])
    for pin in m['execution_policy'].values():verify(pin)
    return {'status':'physical_manifest_validation_pass','tasks':len(m['tasks']),
            'expected_endpoints':len(p['branches']),'unsupported':m['unsupported_endpoints'],
            'paths':len(p['paths']),'manifest':record(manifest)}


def collect(manifest,output):
    validate(manifest);m=read_json(manifest);p=read_json(verify(m['preparation']));rows=[]
    for branch in p['branches']:
        row=dict(branch);path=p['paths'][branch['path_id']];c=p['centers'][path['center_id']]
        row.update(case=c['case'],metal=c['metal'],energy_change_kcal_mol=None,path_gradient_kcal_mol=None)
        if branch['status']=='unsupported':rows.append(row);continue
        t=next(t for t in m['tasks'] if t['task_id']==branch['task_id'])
        try:
            r=endpoint(record(t['output_path']),record(t['output_path']+'.execution.json'),t['xyz'],t['input'])
            gp=read_json(verify(r['receipt']))['artifacts']['engrad'];g=checked_gradient(gp,r,xyz(verify(t['xyz'])))
            _,jac=path_position(xyz(verify(c['xyz'])),path['contact'],branch['fraction']);h=path['contact']['hydrogen_index']
            row.update(status='complete',result=r,gradient=gp,energy_change_kcal_mol=(r['energy_hartree']-c['result']['energy_hartree'])*HA_TO_KCAL,
                path_gradient_kcal_mol=float(g[h]@jac),proton_gradient_kcal_mol_A=g[h].tolist(),
                proton_gradient_norm_kcal_mol_A=float(np.linalg.norm(g[h])))
        except (OSError,ValueError,KeyError) as exc:row.update(status='unavailable',reason=str(exc))
        rows.append(row)
    r={'manifest':record(manifest),'analysis_implementation':record(__file__),'rows':rows,'expected_endpoints':len(p['branches']),
       'status':'complete' if all(r['status']=='complete' for r in rows) else 'incomplete',
       'occupancy_probabilities':None,'pKa':None,'baseline_changed':False}
    write_new(output,r);return {'status':r['status'],'complete':sum(a['status']=='complete' for a in rows),'expected':len(rows)}


def report(collection,output):
    c=read_json(collection);m=read_json(verify(c['manifest']));p=read_json(verify(m['preparation']));out=Path(output);out.mkdir(parents=True,exist_ok=False)
    pairs=[]
    for case,contact in sorted({(a['case'],p['paths'][a['path_id']]['contact']['contact_id']) for a in c['rows']}):
        for fraction in SETTINGS['fractions']:
            rows={a['metal']:a for a in c['rows'] if a['case']==case and p['paths'][a['path_id']]['contact']['contact_id']==contact and a['fraction']==fraction}
            ok=set(rows)=={'Ca','La'} and all(a['status']=='complete' for a in rows.values())
            pairs.append({'case':case,'contact':contact,'fraction':fraction,'status':'complete' if ok else 'unavailable',
                'delta_transfer_Ca_minus_La_kcal_mol':rows['Ca']['energy_change_kcal_mol']-rows['La']['energy_change_kcal_mol'] if ok else None})
    neutral=[]
    for key,path in p['paths'].items():
        center=p['centers'][path['center_id']];atoms=xyz(verify(center['xyz']))
        gradient=checked_gradient(center['gradient'],center['result'],atoms)
        _,jac=path_position(atoms,path['contact'],0.)
        neutral.append({'path_id':key,'case':center['case'],'metal':center['metal'],
            'contact':path['contact']['contact_id'],'fraction':0.,'energy_change_kcal_mol':0.,
            'zero_origin':'exact same endpoint subtraction, not an unavailable correction',
            'path_gradient_kcal_mol':float(gradient[path['contact']['hydrogen_index']]@jac),
            'result':center['result'],'gradient':center['gradient']})
    r={**c,'report_implementation':record(__file__),'balanced_differences':pairs,'reused_neutral_path_origins':neutral}
    write_new(out/'result.json',r)
    lines=['# Fixed-scaffold internal proton-location paths','','Electronic differences only; no pH population or relaxed-basin claim.','',
           '| Case | Contact | Fraction | Metal | ΔE kcal/mol | dE/dt kcal/mol | Status |','|---|---|---:|---|---:|---:|---|']
    for a in c['rows']:lines.append(f"| {a['case']} | {p['paths'][a['path_id']]['contact']['contact_id']} | {a['fraction']} | {a['metal']} | {a['energy_change_kcal_mol']} | {a['path_gradient_kcal_mol']} | {a['status']} |")
    (out/'TABLES.md').write_text('\n'.join(lines)+'\n');return {'status':c['status'],'balanced_pairs':len(pairs)}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='op',required=True)
    for op,fields in [('prepare',('collection','agreement','output')),('validate',('manifest',)),('collect',('manifest','output')),('report',('collection','output'))]:
        q=sub.add_parser(op)
        for name in fields:q.add_argument('--'+name,required=True)
    args=vars(parser.parse_args());operation=args.pop('op');print(json.dumps(globals()[operation](**args),indent=2))
