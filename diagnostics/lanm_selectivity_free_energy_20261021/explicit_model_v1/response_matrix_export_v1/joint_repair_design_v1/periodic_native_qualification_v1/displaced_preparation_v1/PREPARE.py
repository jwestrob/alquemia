"""Exact finite geometry preparation only. No model, PME or solver imports/calls."""
import argparse, ast, hashlib, json, math, os, shutil, socket, time
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
PARENT=HERE.parent
ROOT=next(p for p in HERE.parents if (p/'diagnostics').is_dir() and (p/'workspaces').is_dir())
def require(ok,msg):
    if not ok: raise ValueError(msg)
def pin(path):
    p=Path(path).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def verify(p):
    require(pin(p['path'])==p,'changed input '+p['path']);return Path(p['path'])
def load(p):return json.loads(verify(p).read_text())
def save(p,x):
    with p.open('x') as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
def nearest(x,cell):
    # Literal orthogonal ANINT convention: exact half ties away from zero.
    q=x/cell
    return x-cell*np.copysign(np.floor(np.abs(q)+.5),q)
def parse_xyz(path):
    lines=Path(path).read_text().splitlines(keepends=True)
    n=int(lines[0].split()[0]);cell=np.array(list(map(float,lines[1].split())))
    rows=[s.split() for s in lines[2:]]
    require(len(rows)==n and cell.shape==(6,),'XYZ dimensions')
    require([int(x[0]) for x in rows]==list(range(1,n+1)),'contiguous atom map')
    xyz=np.array([[float(x[k]) for k in (2,3,4)] for x in rows])
    ids=[(x[0],x[1],*x[5:]) for x in rows]
    require(np.all(np.isfinite(xyz)) and np.all(np.isfinite(cell)),'finite geometry')
    return lines,rows,xyz,cell,ids
def write_xyz(path,lines,rows,origin,xyz,cell,oldcell):
    output=[lines[0],lines[1] if np.array_equal(cell,oldcell) else ' '.join(format(v,'.17g') for v in cell)+'\n']
    for i,row in enumerate(rows):
        if np.array_equal(xyz[i],origin[i]):output.append(lines[i+2])
        else:output.append(f'{row[0]:>6s} {row[1]:<4s} '+ ' '.join(format(v,'.17g') for v in xyz[i])+' '+' '.join(row[5:])+'\n')
    with path.open('x') as f:f.write(''.join(output))
def rotate(v,n,theta):
    return v*math.cos(theta)+np.cross(n,v)*math.sin(theta)+np.outer(v@n,n)*(1-math.cos(theta))
def tokens(s):return [x.split('#',1)[0].split() for x in s.splitlines() if x.split('#',1)[0].split()]

def prepare(out):
    tick,cpu=time.monotonic(),time.process_time();out.mkdir(parents=True,exist_ok=False)
    result=dict(status='failed_geometry_preparation',implementation=pin(__file__),checks=[],error=None,new_native_calls=0,
                model_evaluations=0,PME_initializations=0,matrix_assemblies=0,factorizations=0,energy_gradient_calls=0)
    def check(name,ok):result['checks'].append(dict(name=name,passed=bool(ok)));require(ok,name)
    try:
        check('one numerical thread',all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')))
        inp=json.loads((HERE/'INPUTS.json').read_text());result['inputs']=pin(HERE/'INPUTS.json')
        for p in inp['authorities'].values():verify(p)
        deps=load(inp['authorities']['dependencies']);root=load(inp['authorities']['root_contract_review'])
        check('root geometry-only scope',root['status']=='root_adopted_staged_numerical_contract_with_geometry_preparation_only' and root['contract']==deps['contract']==inp['authorities']['contract'] and root['dependencies']==inp['authorities']['dependencies'])
        dep=deps['dependencies']
        for k in ('central_geometry','central_key','BIO_parameters','probe_definition','composition_manifest','central_input_manifest','central_root_static_admission','candidate','candidate_terminal_review'):verify(dep[k])
        state=load(dep['central_input_manifest']);admission=load(dep['central_root_static_admission'])
        check('old La engineering state only',state['candidate']==dep['candidate'] and state['composition_manifest']==dep['composition_manifest'] and state['geometry']==dep['central_geometry'] and state['metal']=='La' and admission['input_manifest']==dep['central_input_manifest'] and admission['status']=='admitted_nominal_saved_La_static_numerical_test')
        lines,rows,xyz,cell,ids=parse_xyz(verify(dep['central_geometry']));n=len(rows)
        types=np.array([int(r[5]) for r in rows]);water=np.flatnonzero(types==402)
        m=int(np.flatnonzero(types==900)[0]);check('identity counts',n==11958 and m==0 and np.sum(types==900)==1 and len(water)==3980 and np.sum(types==403)==7960 and np.sum(types==405)==7 and np.sum(types==414)==10 and set(types)=={900,402,403,405,414})
        check('central box and all-active map',np.array_equal(cell,np.array([49.324]*3+[90.]*3)) and state['shape']==[3*n,3*n] and state['all_alpha_positive'] and state['metal_coordinate_indices']==[0,1,2])
        waters={}
        for oi in water:
            hs=sorted(int(v)-1 for v in rows[oi][6:]);require(len(hs)==2 and all(types[h]==403 and list(map(int,rows[h][6:]))==[oi+1] for h in hs),'water topology')
            waters[int(oi)]=hs
        check('all water topology explicit',len(set(h for hs in waters.values() for h in hs))==7960)
        check('ANINT half ties and signed images',np.array_equal(nearest(np.array([22.,-22.,66.,-66.]),np.array([44.]*4)),np.array([-22.,22.,-22.,22.])))
        displacement=nearest(xyz[water]-xyz[m],cell[:3]);dist2=np.einsum('ij,ij->i',displacement,displacement)
        oi=int(min(zip(dist2,water),key=lambda z:(float(z[0]),int(z[1])))[1]);rot_ids=[oi]+waters[oi]
        local=np.vstack([np.zeros(3)]+[nearest(xyz[h]-xyz[oi],cell[:3]) for h in waters[oi]])
        centroid=np.mean(local,axis=0);u=nearest(xyz[oi]-xyz[m],cell[:3]);axis=np.cross(u,local[1]-local[2]);norm=float(np.linalg.norm(axis))
        check('rotation axis defined',math.isfinite(norm) and norm>0);axis/=norm;lever=local-centroid
        check('fixed rotation central identity',np.array_equal(xyz[rot_ids]+(rotate(lever,axis,0.)-lever),xyz[rot_ids]))
        # Reuse only the pinned pure probe function AST; never import the backend module.
        tree=ast.parse(verify(dep['probe_definition']).read_text());func=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='ten_probes')
        ns={'np':np};exec(compile(ast.Module(body=[func],type_ignores=[]),str(dep['probe_definition']['path']),'exec'),ns)
        probes=ns['ten_probes'](3*n,(0,1,2));whole=dict(probes)['rademacher_0'].reshape(n,3)/math.sqrt(3.)
        check('frozen whole-box coefficients',len(probes)==10 and set(np.unique(dict(probes)['rademacher_0']))=={-1.,1.} and np.allclose(np.linalg.norm(whole,axis=1),1.,rtol=0,atol=4*np.finfo(float).eps))
        prefix=verify(dep['central_key']).read_text();derive=ast.parse(verify(inp['authorities']['key_settings']).read_text())
        settings=next(ast.literal_eval(x.value) for x in derive.body if isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='SETTINGS' for t in x.targets))
        appended='\n# Prepared finite periodic E/F geometry only; execution remains unbound.\n'+''.join(line+'\n' for line in settings.splitlines() if line.strip() and not line.startswith('#') and not line.upper().startswith('PME-NORM-EXPORT'))
        key=prefix+appended;kt=tokens(key);base=tokens(prefix);added=tokens(appended)
        check('exact key prefix and scientific settings',kt==base+added and len([t for t in base if t[0].lower()=='polar-local-response'])==2 and not any(t[0].lower() in {'ewald-boundary','pme-norm-export'} for t in kt))
        expected=[['ewald'],['ewald-cutoff','9.0'],['vdw-cutoff','9.0'],['vdw-taper','8.1'],['ewald-alpha','0.4'],['pewald-alpha','0.4'],['pme-grid','64','64','64'],['pme-order','5'],['ppme-order','5'],['fft-package','FFTW'],['polarization','MUTUAL'],['polar-iter','100']]
        check('only frozen periodic settings appended',added==expected and ['polar-eps','1.0e-8'] in base and ['openmp-threads','1'] in base)
        center=out/'central';center.mkdir();shutil.copyfile(verify(dep['central_geometry']),center/'system.xyz');shutil.copyfile(verify(dep['BIO_parameters']),center/'amoebabio09.prm');(center/'system.prepared.key').write_text(key)
        check('literal central geometry/parameters retained',pin(center/'system.xyz')['sha256']==dep['central_geometry']['sha256'] and pin(center/'amoebabio09.prm')['sha256']==dep['BIO_parameters']['sha256'])
        directions=[]
        for a in range(3):
            v=np.zeros_like(xyz);v[m,a]=1.;directions.append(dict(id='metal_'+str(a),kind='translation',axis=a,tangent=v))
        v=np.zeros_like(xyz);v[rot_ids]=np.cross(axis,lever);directions.append(dict(id='water_rotation',kind='rotation',atom_indices_1based=[i+1 for i in rot_ids],axis=axis.tolist(),local_centroid_A=centroid.tolist(),local_vectors_A=local.tolist(),tangent=v))
        directions.append(dict(id='whole_box',kind='linear',probe='rademacher_0',normalization='divide by sqrt(3)',tangent=whole))
        for a in range(3):
            v=np.zeros_like(xyz);v[:,a]=xyz[:,a]-xyz[m,a];S=np.zeros((3,3));S[a,a]=1.;directions.append(dict(id='strain_'+str(a),kind='strain',axis=a,strain_matrix=S.tolist(),anchor_A=xyz[m].tolist(),L0_A=1.,tangent=v))
        prepared=[];dd=out/'directions';dd.mkdir();eps=np.finfo(float).eps
        for direction in directions:
            v=direction.pop('tangent');fp=dd/(direction['id']+'.npy');np.save(fp,v,allow_pickle=False);direction['coefficients_dr_dq']=pin(fp)
            check(direction['id']+' fixed finite tangent',v.shape==(n,3) and np.all(np.isfinite(v)) and np.array_equal(np.load(fp,allow_pickle=False),v))
            for step in (1e-4,5e-5):
                for sign in (-1,1):
                    q=sign*step;new=xyz+q*v;newcell=cell.copy()
                    if direction['kind']=='rotation':new=xyz.copy();new[rot_ids]=xyz[rot_ids]+(rotate(lever,axis,q)-lever)
                    if direction['kind']=='strain':newcell[direction['axis']]*=1+q
                    name=f"{direction['id']}__{step:.0e}__{'minus' if sign<0 else 'plus'}";folder=out/'decks'/name;folder.mkdir(parents=True)
                    write_xyz(folder/'system.xyz',lines,rows,xyz,new,newcell,cell);(folder/'system.prepared.key').write_text(key);shutil.copyfile(center/'amoebabio09.prm',folder/'amoebabio09.prm')
                    _,rr,got,gcell,gids=parse_xyz(folder/'system.xyz')
                    check(name+' exact coordinate/cell roundtrip and identity',np.array_equal(got,new) and np.array_equal(gcell,newcell) and gids==ids)
                    check(name+' domain and cutoff margins',np.all(gcell[3:]==90.) and min(gcell[:3])>=44 and 21<min(gcell[:3])/2 and 9<min(gcell[:3])/2)
                    if direction['kind']!='rotation':check(name+' signed fixed linear coefficients',np.max(np.abs((got-xyz)-q*v))<=128*eps*max(1.,float(np.max(np.abs(xyz))),float(np.max(np.abs(q*v)))))
                    else:
                        dq=got[rot_ids]-xyz[rot_ids];rot_local=local+dq
                        pairs=[(0,1),(0,2),(1,2)];err=max(abs(float(np.linalg.norm(rot_local[a]-rot_local[b])-np.linalg.norm(local[a]-local[b]))) for a,b in pairs)
                        check(name+' intact rotation and fixed centroid',err<=128*eps*max(1.,float(np.max(np.abs(xyz)))) and np.max(np.abs(np.mean(dq,axis=0)))<=128*eps*max(1.,float(np.max(np.abs(xyz)))))
                        radial=lever-np.outer(lever@axis,axis);ri=int(np.argmax(np.linalg.norm(radial,axis=1)));a=radial[ri];b=(lever+dq)[ri];b=b-axis*float(b@axis)
                        observed=math.atan2(float(axis@np.cross(a,b)),float(a@b));angular_bound=128*eps*max(1.,float(np.max(np.abs(xyz))))/float(np.linalg.norm(a))
                        check(name+' signed angle equals declared q/L0',abs(observed-q)<=angular_bound)
                    check(name+' unchanged key and parameter bytes',pin(folder/'system.prepared.key')['sha256']==pin(center/'system.prepared.key')['sha256'] and pin(folder/'amoebabio09.prm')['sha256']==dep['BIO_parameters']['sha256'])
                    prepared.append(dict(id=name,direction=direction['id'],q_A=q,half_step_A=step,sign=sign,cell_A=gcell.tolist(),input={k:pin(folder/f) for k,f in [('xyz','system.xyz'),('key','system.prepared.key'),('parameters','amoebabio09.prm')]},minimum_cell_side_A=float(min(gcell[:3])),compact_halfbox_clearance_A=float(min(gcell[:3])/2-21),native_initialization=None,PME_export=None,operator_assembly=None,stability_certificate=None,native_energy_gradient=None,native_output=None,root_execution_release=None))
        save(dd/'DIRECTIONS.json',dict(status='fixed_central_direction_coefficients',units='dr/dq dimensionless; q in Angstrom; rotation theta=q/(1 Angstrom)',central_geometry=dep['central_geometry'],active_atom_order='native XYZ index order, all sites positive alpha',probe_source=dep['probe_definition'],probe_function_AST_sha256=hashlib.sha256(ast.dump(func,include_attributes=False).encode()).hexdigest(),selected_water_O_atom=oi+1,nearest_water_distance_A=math.sqrt(float(min(dist2))),axis_definition='normalize[minimum_image(O-M) cross (local_H1-local_H2)]',directions=directions))
        check('exact four endpoints each eight directions',len(prepared)==32 and len({r['id'] for r in prepared})==32 and all(sorted(r['q_A'] for r in prepared if r['direction']==d['id'])==[-1e-4,-5e-5,5e-5,1e-4] for d in directions))
        manifest=dict(status='prepared_geometry_only_not_executable',implementation=pin(__file__),inputs=pin(HERE/'INPUTS.json'),contract=deps['contract'],root_preparation_release=inp['authorities']['root_contract_review'],candidate=dep['candidate'],composition=dep['composition_manifest'],central_original_geometry=dep['central_geometry'],central_prepared_input={k:pin(center/f) for k,f in [('xyz','system.xyz'),('key','system.prepared.key'),('parameters','amoebabio09.prm')]},directions=pin(dd/'DIRECTIONS.json'),rows=prepared,count=32,directions_count=8,all_endpoint_native_and_certificate_bindings_null=True,new_native_calls=0,new_physical_reference_calls=0,PME_initializations=0,assemblies=0,factorizations=0,energy_gradient_calls=0,stage2_execution_released=False)
        save(out/'MANIFEST.json',manifest);result.update(status='passed_geometry_only_preparation',manifest=pin(out/'MANIFEST.json'),directions=manifest['directions'],coordinate_cell_decks=32,fixed_directions=8,selected_water_O_atom=oi+1,minimum_cell_side_A=min(r['minimum_cell_side_A'] for r in prepared),minimum_compact_halfbox_clearance_A=min(r['compact_halfbox_clearance_A'] for r in prepared),future_bindings_null=True)
    except Exception as exc:result['error']=dict(type=type(exc).__name__,message=str(exc))
    result['cost']=dict(wall_seconds=time.monotonic()-tick,CPU_seconds=time.process_time()-cpu,numerical_threads=1,hostname=socket.gethostname(),PID=os.getpid(),job_id=None)
    save(out/'RESULT.json',result);print(json.dumps({k:result[k] for k in ('status','error','cost')}));print(json.dumps(pin(out/'RESULT.json')))
    return 0 if result['error'] is None else 1
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);args=p.parse_args();raise SystemExit(prepare(args.out.resolve()))
