#!/usr/bin/env python3
"""One-thread saved force projections; no molecular engine imports or calls."""
import argparse, hashlib, json, math, os, platform, time
from pathlib import Path

THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
ROOT_RESULT='e096e732d70a7e352466ab8a1dc46ff19de8427416a0f21e6dcce3a960352c7c'
EV_KCAL=23.060547830619
def pin(p):
    p=Path(p).resolve(); return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path'])==p, p; return Path(p['path'])
def read(p): return json.loads(Path(p).read_text())
def save(p,r): Path(p).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def flat(a): return [x for v in a for x in v]
def sub(a,b): return [x-y for x,y in zip(a,b)]
def norm(a): return math.sqrt(dot(a,a))
def cross(a,b): return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def close(a,b): assert abs(a-b)<=1e-8+1e-8*max(abs(a),abs(b)),(a,b)
def metrics(f):
    v=flat(f);s=dot(v,v);return {'component_rmse_kcal_mol_A':math.sqrt(s/len(v)),'squared_norm':s}
def modes(x,f):
    n=len(x);center=[sum(v[k] for v in x)/n for k in range(3)]
    r=[sub(v,center) for v in x];mean=[sum(v[k] for v in f)/n for k in range(3)]
    t=mean*n;v=flat(f);q=sub(v,t);basis=[]
    for e in ([1,0,0],[0,1,0],[0,0,1]):
        b=flat([cross(e,ri) for ri in r])
        for u in basis:b=sub(b,[dot(b,u)*z for z in u])
        length=norm(b)
        if length>1e-12:basis.append([z/length for z in b])
    assert len(basis)==3
    rot=[sum(dot(q,u)*u[k] for u in basis) for k in range(3*n)];internal=sub(q,rot)
    parts={'translation':t,'rotation':rot,'internal':internal}
    for i in range(3*n):close(v[i],sum(z[i] for z in parts.values()))
    for a,b in [('translation','rotation'),('translation','internal'),('rotation','internal')]:close(dot(parts[a],parts[b]),0)
    close(dot(v,v),sum(dot(z,z) for z in parts.values()))
    return {k:{'squared_norm':dot(z,z),'component_rmse_kcal_mol_A':math.sqrt(dot(z,z)/(3*n))} for k,z in parts.items()}
def groups(g):
    n=g['atoms'];metal=[g['metal_index']];ligand=g.get('ligand_indices',[])
    rem=set(range(n))-set(metal)-set(ligand);waters=[];b=g['source_bonds_zero_based']
    while rem:
        seed=min(rem);todo=[seed];comp=set()
        while todo:
            i=todo.pop()
            if i not in comp:comp.add(i);todo.extend(j for j in b[i] if j not in comp)
        assert comp<=rem
        oxygen=[i for i in comp if g['symbols'][i]=='O'];hydrogen=sorted(i for i in comp if g['symbols'][i]=='H')
        assert len(oxygen)==1 and len(hydrogen)==2 and len(comp)==3
        waters.append(oxygen+hydrogen);rem-=comp
    assert len(waters)==g['water_count']
    cat={'metal':metal,'ligand':ligand,'water':[i for w in waters for i in w]}
    assert sorted(i for v in cat.values() for i in v)==list(range(n))
    return cat,waters
def water_bonded(x):
    # Exact SI water classes34/35 and native polynomial conventions.
    f={k:[[0.0]*3 for _ in range(3)] for k in ('Bond Stretching','Angle Bending','Urey-Bradley')}
    energy={k:0.0 for k in f};u=[sub(x[i],x[0]) for i in (1,2)];r=list(map(norm,u));unit=[[v/z for v in a] for a,z in zip(u,r)]
    for i,(rr,uu) in enumerate(zip(r,unit),1):
        d=rr-.9572;energy['Bond Stretching']+=556.85*d*d*(1-2.55*d+3.793125*d*d)
        deriv=556.85*(2*d-3*2.55*d*d+4*3.793125*d**3)
        for k in range(3):f['Bond Stretching'][i][k]-=deriv*uu[k];f['Bond Stretching'][0][k]+=deriv*uu[k]
    cosine=max(-1,min(1,dot(*unit)));theta=math.acos(cosine);radian=180/math.pi;dt=theta*radian-108.5
    coeff=[1,-.014,.000056,-.0000007,.000000022];K=48.70/(radian*radian)
    energy['Angle Bending']=K*sum(c*dt**(j+2) for j,c in enumerate(coeff));deriv=K*radian*sum((j+2)*c*dt**(j+1) for j,c in enumerate(coeff))
    for i in (0,1):
        grad=[(cosine*unit[i][k]-unit[1-i][k])/(r[i]*math.sin(theta)) for k in range(3)]
        for k in range(3):f['Angle Bending'][i+1][k]-=deriv*grad[k];f['Angle Bending'][0][k]+=deriv*grad[k]
    hh=sub(x[1],x[2]);d=norm(hh);delta=d-1.5337;energy['Urey-Bradley']=-7.60*delta*delta
    for k in range(3):v=-2*(-7.60)*delta*hh[k]/d;f['Urey-Bradley'][1][k]+=v;f['Urey-Bradley'][2][k]-=v
    total=[[sum(a[i][k] for a in f.values()) for k in range(3)] for i in range(3)]
    return {'OH_A':r,'HOH_deg':theta*radian,'HH_A':d},energy,total,{k:metrics(a) for k,a in f.items()}
def describe(g,f):
    cat,waters=groups(g);total=metrics(f);out={'all':total,'categories':{},'water_modes':{},'waters':[]}
    for k,ids in cat.items():
        if ids:out['categories'][k]={'atoms':len(ids),**metrics([f[i] for i in ids])};out['categories'][k]['fraction_of_total_squared_norm']=out['categories'][k]['squared_norm']/total['squared_norm'] if total['squared_norm'] else None
    close(sum(v['squared_norm'] for v in out['categories'].values()),total['squared_norm'])
    for ids in waters:
        m=modes([g['positions_A'][i] for i in ids],[f[i] for i in ids]);out['waters'].append({'indices':ids,'modes':m})
    for name in ('translation','rotation','internal'):
        s=sum(w['modes'][name]['squared_norm'] for w in out['waters']);out['water_modes'][name]={'squared_norm':s,'component_rmse_kcal_mol_A':math.sqrt(s/(9*len(waters))),'fraction_of_total_squared_norm':s/total['squared_norm'] if total['squared_norm'] else None}
    close(sum(v['squared_norm'] for v in out['water_modes'].values()),out['categories']['water']['squared_norm'])
    if cat['ligand']:out['ligand_modes']=modes([g['positions_A'][i] for i in cat['ligand']],[f[i] for i in cat['ligand']])
    return out
def main():
    a=argparse.ArgumentParser();a.add_argument('--root',required=True);a.add_argument('--output',required=True);args=a.parse_args();start=time.monotonic()
    assert all(os.environ.get(k)=='1' for k in THREADS)
    b=Path(args.root).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
    result={'status':'unavailable','implementation':pin(__file__),'specification':pin(Path(__file__).with_name('WATER_FORCE_ATTRIBUTION_SPEC.md')),'new_molecular_calls':0,'threads':1,'rows':[],'paired':[]}
    try:
        cp=b/'hydrated_comparison_v1/partial_v1/RESULT.json';assert pin(cp)['sha256']==ROOT_RESULT;c=read(cp);result['parent']=pin(cp)
        for p in c['inputs']:verify(p)
        result['inputs']=c['inputs'];assert not c['missing_references'] and not c['failed_references']
        ar=read(b/'explicit_model_v1/hydrated_mixed_run_v1/RESULT.json');assert ar['complete_rows']==20
        rows={r['id']:r for r in ar['rows']};ids=[f'hydrated_{d}_central_{m}' for d in ('acetate_mono','acetate_bi','acetamide') for m in ('La','Dy')]+[f'aquo8_central_{m}' for m in ('La','Dy')]
        gs={gid:rows[gid]['input']['source_geometry'] for gid in ids}
        # Original aquo metadata carries the same graph under its original field name.
        for gid,g in gs.items():
            if 'source_bonds_zero_based' not in g:g['source_bonds_zero_based']=rows[gid]['input']['bonds_zero_based']
        data={'AMOEBA_transfer':{},'MACE':{},'ecp_pbe0':{}}
        for gid in ids:
            r=rows[gid];assert r['status']=='complete'
            for key in ('xyz','key','parameters'):verify(r['input'][key])
            for key in ('analyze_output','gradient_output'):verify(r['native'][key])
            data['AMOEBA_transfer'][gid]=r['native']['force_kcal_mol_A']
        for rel in ('mace_chemistry_v1/submission_v1/execution/RESULT.json','hydrated_chemistry_run_v1/mace/execution/RESULT.json'):
            for r in read(b/rel)['rows']:
                if r['id'] in ids:assert r['status']=='complete';data['MACE'][r['id']]=[[v*EV_KCAL for v in f] for f in r['forces_eV_A']]
        for p in c['inputs']:
            if '/tasks/ecp_pbe0__' in p['path']:
                r=read(p['path']);assert r['status']=='complete' and all(r['checks'].values());gid=r['cell_id'].removeprefix('ecp_pbe0__')
                for key in ('output','engrad','input','geometry'):verify(r[key])
                data['ecp_pbe0'][gid]=r['forces_kcal_mol_A']
        assert all(set(v)==set(ids) for v in data.values())
        # Validate source constants, no key overrides of water bonded terms.
        inp=rows[ids[0]]['input'];prm=verify(inp['parameters']).read_text();key=verify(inp['key']).read_text();required=['bond         34   35          556.85     0.9572','angle        35   34   35      48.70     108.50','ureybrad     35   34   35      -7.60     1.5337']
        assert all(s in prm for s in required);assert not any(s.split()[:4]==['strbnd','35','34','35'] for s in prm.splitlines())
        for line in key.splitlines():
            s=line.split()
            if s and s[0] in ('bond','angle','ureybrad','strbnd'):assert '34' not in s[1:4] and '35' not in s[1:4]
        result['water_bonded_parameter_lines']=required;src=b/'explicit_model_v1/tinker_source/source';result['formula_sources']=[pin(src/(name+'.f')) for name in ('initprm','ebond1','eangle1','eurey1')]
        result['water_geometry_and_bonded']=[];errors={m:{} for m in ('AMOEBA_transfer','MACE')}
        for gid in ids:
            g=gs[gid];cat,waters=groups(g);bonded=[[0.0]*3 for _ in range(g['atoms'])];energies={k:0.0 for k in ('Bond Stretching','Angle Bending','Urey-Bradley')};waterrecords=[]
            for w in waters:
                geo,e,f,terms=water_bonded([g['positions_A'][i] for i in w]);waterrecords.append({'indices':w,'geometry':geo,'bonded_force_terms':terms})
                for k,v in e.items():energies[k]+=v
                for i,ff in zip(w,f):bonded[i]=ff
            if gid.startswith('aquo8'):
                for k,v in energies.items():assert abs(v-rows[gid]['native']['components'][k]['energy_kcal_mol'])<=1e-7,(gid,k,v)
            result['water_geometry_and_bonded'].append({'id':gid,'waters':waterrecords,'analytic_water_bonded_components_kcal_mol':energies,'water_bonded_force':describe(g,bonded)})
            for method in data:result['rows'].append({'id':gid,'kind':'force','model':method,**describe(g,data[method][gid])})
            for method in errors:
                f=[sub(a,z) for a,z in zip(data[method][gid],data['ecp_pbe0'][gid])];errors[method][gid]=f;desc=describe(g,f)
                if gid.startswith('hydrated'):
                    donor=gid.removeprefix('hydrated_').rsplit('_central_',1)[0];old=next(r for r in c['force_errors'] if r['model']==method and r['donor']==donor);close(desc['all']['component_rmse_kcal_mol_A'],old[g['metal']+'_component_rmse_kcal_mol_A'])
                result['rows'].append({'id':gid,'kind':'error','model':method,**desc})
        for family in ('hydrated_acetate_mono','hydrated_acetate_bi','hydrated_acetamide','aquo8'):
            la,dy=[family+'_central_'+m for m in ('La','Dy')];g=gs[la];assert g['positions_A']==gs[dy]['positions_A']
            for method in errors:
                f=[sub(a,z) for a,z in zip(errors[method][dy],errors[method][la])];desc=describe(g,f)
                if family.startswith('hydrated'):
                    old=next(r for r in c['force_errors'] if r['model']==method and r['donor']==family.removeprefix('hydrated_'));close(desc['all']['component_rmse_kcal_mol_A'],old['paired_component_rmse_kcal_mol_A'])
                result['paired'].append({'family':family,'model':method,'kind':'Dy_minus_La_error',**desc})
        result.update(status='complete',reference_accuracy_qualified=False,fit_or_production_release=False)
    except Exception as e:result['error']=repr(e)
    result['elapsed_seconds']=time.monotonic()-start;result['python']=platform.python_version();save(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json'),'error':result.get('error')}))
    return 0 if result['status']=='complete' else 1
if __name__=='__main__':raise SystemExit(main())
