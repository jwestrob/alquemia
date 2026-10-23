#!/usr/bin/env python3
"""Root-released, serial native representability pilot; no automatic extension."""
import argparse, hashlib, importlib.util, json, math, os, shutil, time
from pathlib import Path

THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path'])==p,'Changed pin: '+p['path'];return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(path,obj):Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def key_id(p):return hashlib.sha256(json.dumps(p,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def jacobi_eigenvalues(a):
    a=[r[:] for r in a];n=len(a)
    for _ in range(100*n*n):
        i,j=max(((i,j) for i in range(n) for j in range(i+1,n)),key=lambda z:abs(a[z[0]][z[1]]))
        if abs(a[i][j])<=1e-14*max(1.,max(abs(a[k][k]) for k in range(n))):break
        angle=.5*math.atan2(2*a[i][j],a[j][j]-a[i][i]);c=math.cos(angle);s=math.sin(angle)
        x,y,z=a[i][i],a[j][j],a[i][j]
        a[i][i]=c*c*x-2*c*s*z+s*s*y;a[j][j]=s*s*x+2*c*s*z+c*c*y;a[i][j]=a[j][i]=0.
        for k in range(n):
            if k in (i,j):continue
            x,y=a[k][i],a[k][j];a[k][i]=a[i][k]=c*x-s*y;a[k][j]=a[j][k]=s*x+c*y
    return sorted((max(0.,a[i][i]) for i in range(n)),reverse=True)

def parameter_key(base,family,values):
    for p,v in zip(family['parameters'],values):assert p['lower']<=v<=p['upper']
    v={p['name']:x for p,x in zip(family['parameters'],values)}
    lines=[]
    for name,cl,ty in (('water',34,36),('carboxylate',403,403),('amide',55,105)):
        eps=v['water_epsilon'] if name=='water' else family['fixed_pair_epsilon'][name]
        lines += [f"vdwpr {cl} 501 {v[name+'_radius']:.17g} {eps:.17g}",f"polpair {ty} 501 {v[name+'_damping']:.17g}"]
    return base.rstrip()+'\n\n# Frozen seven-parameter Dy-only representability pilot\n'+'\n'.join(lines)+'\n'

def run(config_path,out):
    cfg=json.loads(Path(config_path).read_text());out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    result={'status':'unavailable','implementation':pin(__file__),'config':pin(config_path),'native_program_invocations':0,
            'native_evaluations':{},'stages':{},'validation_outcomes_opened':False,'physical_accuracy_qualified':False,
            'production_released':False,'candidate_frozen':False,'job_id':os.environ.get('SLURM_JOB_ID'),'error':None}
    try:
        assert cfg['implementation']==pin(__file__) and os.environ.get('SLURM_JOB_ID','').isdigit()
        assert int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1 and int(os.environ.get('SLURM_NTASKS','1'))==1
        assert all(os.environ.get(k)=='1' for k in THREADS)
        data=read(cfg['training']);assert data['status']=='training_complete' and data['complete']==28 and data['saved_controls_passed']
        family=read(cfg['family']);assert data['family']==cfg['family'] and data['objective']==cfg['objective'] and data['specification']==cfg['specification']
        verify(cfg['specification']);assert family['specification']==cfg['specification']
        assert family['native_program_invocation_cap']==1276 and family['water_vector_cap']==25 and family['donor_vector_cap']==49
        training_cfg=read(data['config']);verify(data['implementation'])
        assert training_cfg['native_helper']==cfg['native_helper']
        # Hash all admitted raw sources again. No reference engine or preparation rerun.
        for r in data['rows'].values():
            verify(r['baseline_result']);verify(r['reference_result'])
            for k in ('xyz','key','parameters'):verify(r['input'][k])
            for k in ('analyze_output','gradient_output'):verify(r['baseline']['native'][k])
            for k in ('output','engrad','input','panel'):verify(r['reference'][k])
        release=read(cfg['root_release'])
        for k in ('implementation','training','family','objective','specification'):
            assert release[k]==cfg[k]
        assert release['status']=='approved_finite_development_representability_pilot'
        assert release['maximum_native_program_invocations']==1276 and release['validation_scoring_authorized'] is False
        assert release['uncorrected_PBE0_energy_and_forces_only'] is True and release['physical_accuracy_claimed'] is False
        assert release['radial_screen_all_passed']==data['radial_screen_all_passed']
        assert release['La_controls_passed']=={s:data['baseline_objectives'][s]['La_control_passed'] for s in ('water','donor')}
        assert isinstance(release['reference_eligibility_rationale'],str) and release['reference_eligibility_rationale'].strip()
        if not data['radial_screen_all_passed']:
            assert release['coarse_screen_disposition']=='explicit_truncation_limited_development_only'
        obj=module(cfg['objective'],'pilot_objective');native=module(cfg['native_helper'],'pilot_native')
        for p in cfg['executables'].values():verify(p)
        base=verify(family['base_key']).read_text();verify(family['base_parameters'])
        parameters=family['parameters'];assert len(parameters)==7
        initial=[p['initial'] for p in parameters];model_base=obj.baseline(data)
        cache={}
        # Cache authoritative original candidate values; no repeated chemistry.
        for gid,row in data['rows'].items():
            if gid.endswith('_Dy'):cache[key_id({'parameters':initial,'cell':gid})]=row['baseline']
        def evaluate_native(values,gid):
            ck=key_id({'parameters':values,'cell':gid})
            if ck in cache:return cache[ck]
            row=data['rows'][gid];assert gid.endswith('_Dy')
            directory=out/'native'/ck;directory.mkdir(parents=True,exist_ok=False)
            rec={'cell':gid,'parameters':values[:],'status':'unavailable','receipts':[],'source_input':row['input']}
            result['native_evaluations'][ck]=rec
            try:
                shutil.copyfile(verify(row['input']['xyz']),directory/'system.xyz')
                shutil.copyfile(verify(row['input']['parameters']),directory/'amoeba09.prm')
                (directory/'system.key').write_text(parameter_key(base,family,values))
                rec['input']={k:pin(directory/name) for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm'))}
                assert rec['input']['xyz']['sha256']==row['input']['xyz']['sha256'] and rec['input']['parameters']['sha256']==family['base_parameters']['sha256']
                for name,argv in [('analyze',[cfg['executables']['analyze']['path'],'system.xyz','E']),
                                  ('gradient',[cfg['executables']['testgrad']['path'],'system.xyz','Y','N'])]:
                    assert result['native_program_invocations']<1276
                    result['native_program_invocations']+=1
                    native.invoke(argv,directory,directory/(name+'.out'),rec['receipts'])
                metadata=dict(row['input'],derivative_test=False)
                raw=native.parse(directory/'analyze.out',directory/'gradient.out',metadata)
                rec.update(status='complete',native=raw)
                value={'energy':raw['total_energy_kcal_mol'],'forces':raw['force_kcal_mol_A'],'native':raw}
                cache[ck]=value;return value
            except Exception as exc:
                rec['error']=repr(exc);raise
            finally:save(out/'RESULT.json',result)
        selected=initial[:];water_model={}
        for stage in ('water','donor'):
            indices=[i for i,p in enumerate(parameters) if p['stage']==stage];cap=family[stage+'_vector_cap'];search_cap=cap-2*len(indices)
            qs=obj.WATER if stage=='water' else obj.DONOR
            records={};stage_result={'status':'running','vector_cap':cap,'search_cap':search_cap,'evaluations':records,'polls':[]}
            result['stages'][stage]=stage_result
            def from_unit(z):
                values=selected[:]
                for i,v in zip(indices,z):values[i]=parameters[i]['lower']+v*(parameters[i]['upper']-parameters[i]['lower'])
                return values
            z0=[(selected[i]-parameters[i]['lower'])/(parameters[i]['upper']-parameters[i]['lower']) for i in indices]
            def evaluate(z,role):
                token=key_id(z)
                if token in records:return records[token]
                assert len(records)<cap
                values=from_unit(z)
                # Exact original baseline is an authoritative cache entry, not a
                # decimal-roundtrip-generated new parameter set.
                if z==z0: values=selected[:]
                rec={'normalized_parameters':z[:],'parameters':values,'role':role,'status':'unavailable','cells':{q+'_Dy':'unavailable' for q in qs}}
                records[token]=rec;model=dict(model_base);model.update(water_model)
                for q in qs:
                    gid=q+'_Dy';model[gid]=evaluate_native(values,gid);rec['cells'][gid]='complete'
                score=obj.evaluate(data,model,stage);rec.update(status='complete',score=score)
                save(out/'RESULT.json',result);return rec
            best=evaluate(z0,'initial');step=.10;stop='budget_limited'
            while len(records)<search_cap:
                z=best['normalized_parameters'];polled=[];complete=True
                for j in range(len(indices)):
                    for sign in (-1,1):
                        q=z[:];q[j]=max(0.,min(1.,q[j]+sign*step))
                        if q==z:continue
                        if key_id(q) not in records and len(records)>=search_cap:complete=False;break
                        polled.append(evaluate(q,'search'))
                    if not complete:break
                chosen=best
                for rec in polled:
                    if rec['score']['objective']<chosen['score']['objective']:chosen=rec
                improved=chosen is not best;best=chosen
                stage_result['polls'].append({'step':step,'complete':complete,'improved':improved,'best_objective':best['score']['objective']})
                if not complete:break
                if not improved:step/=2
                if step<.00625:stop='local_step_rule_met';break
                if not polled:stop='no_distinct_poll_points';break
            search_best=best;cols=[]
            for j in range(len(indices)):
                lo=best['normalized_parameters'][:];hi=lo[:];lo[j]=max(0.,lo[j]-.02);hi[j]=min(1.,hi[j]+.02)
                a=evaluate(lo,'sensitivity')['score']['weighted_residuals'];b=evaluate(hi,'sensitivity')['score']['weighted_residuals']
                cols.append([(v-u)/(hi[j]-lo[j]) for u,v in zip(a,b)])
            gram=[[sum(a*b for a,b in zip(x,y)) for y in cols] for x in cols]
            singular=[math.sqrt(v) for v in jacobi_eigenvalues(gram)];top=singular[0]
            rank=sum(v>top*1e-6 for v in singular) if top>0 else 0
            condition=top/singular[-1] if singular[-1]>0 else None
            bound=any(min(z,1-z)<=.01 for z in best['normalized_parameters'])
            # Diagnostic perturbations never replace the selected search point.
            stage_result.update(status=stop,best=search_best,distinct_vectors=len(records),
                near_bound=bound,sensitivity={'singular_values':singular,'rank':rank,'columns':len(indices),
                    'column_norms':[math.sqrt(sum(v*v for v in c)) for c in cols],'condition':condition,
                    'weakly_identified':rank<len(indices) or condition is None or condition>1e6},
                numerical_optimization_complete=stop=='local_step_rule_met')
            selected=best['parameters'][:]
            if stage=='water':
                for q in obj.WATER:water_model[q+'_Dy']=evaluate_native(selected,q+'_Dy')
            save(out/'RESULT.json',result)
        result.update(status='finite_pilot_complete',candidate_parameters=selected,
                      parameter_names=[p['name'] for p in parameters],candidate_frozen=True,
                      optimizer_budget_limited=any(s['status']=='budget_limited' for s in result['stages'].values()),
                      controls_passed=all(s['best']['score']['La_control_passed'] and s['best']['score']['candidate_development_passed'] for s in result['stages'].values()),
                      identifiability_and_bounds_passed=all(not s['near_bound'] and not s['sensitivity']['weakly_identified'] for s in result['stages'].values()),
                      automatic_validation_release=False)
        (out/'CANDIDATE.key').write_text(parameter_key(base,family,selected));result['candidate_key']=pin(out/'CANDIDATE.key')
    except Exception as exc:result['error']=repr(exc);result['status']='failed_or_unavailable'
    finally:
        result['elapsed_seconds']=time.monotonic()-start;save(out/'RESULT.json',result)
    print(json.dumps({'status':result['status'],'native_program_invocations':result['native_program_invocations'],'result':pin(out/'RESULT.json')}))
    return 0 if result['status']=='finite_pilot_complete' else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();raise SystemExit(run(a.config,a.output))
