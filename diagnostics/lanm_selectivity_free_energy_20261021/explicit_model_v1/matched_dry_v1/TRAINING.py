#!/usr/bin/env python3
"""Admit fixed saved training E/F, compare coarse derivatives. Zero native calls."""
import argparse, datetime, hashlib, importlib.util, json, math, os, re, time
from pathlib import Path

THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')
def pin(path):
    p=Path(path).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path'])==p, 'Changed pin: '+p['path'];return Path(p['path'])
def read(p): return json.loads(verify(p).read_text())
def save(path,obj): Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def module(p,name):
    s=importlib.util.spec_from_file_location(name,verify(p));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def finite(a): return all(math.isfinite(x) for v in a for x in v)

def prepare(config_path,out):
    started=time.monotonic();cfg=json.loads(Path(config_path).read_text());out.mkdir(parents=True,exist_ok=False)
    result={'status':'unavailable','implementation':pin(__file__),'config':pin(config_path),'new_molecular_calls':0,
            'local_threads':1,'rows':{k:{'status':'unavailable'} for k in cfg['reference_results']},
            'fit_released':False,'validation_outcomes_opened':False,'physical_accuracy_qualified':False}
    try:
        assert cfg['implementation']==pin(__file__) and len(result['rows'])==28
        assert all(os.environ.get(k)=='1' for k in THREADS)
        verify(cfg['specification']);family=read(cfg['family']);verify(cfg['reserved_geometry_manifest'])
        obj=module(cfg['objective'],'repair_objective');native=module(cfg['native_helper'],'repair_native');ref=module(cfg['reference_helper'],'repair_reference')
        results={}
        for rp in cfg['baseline_results']:
            r=read(rp); verify(r['config']);verify(r['implementation'])
            assert r.get('status') in ('panel_complete','complete','hydrated_complete','matched_dry_complete'),r.get('status')
            for row in r['rows']:
                if row['id'] in result['rows']:
                    if row['id'] in results:
                        assert results[row['id']][0]['native']==row['native'];continue
                    results[row['id']]=(row,rp)
        assert set(results)==set(result['rows'])
        for gid,rp in cfg['reference_results'].items():
            b,bp=results[gid];assert b['status']=='complete';i=b['input'];g=i['source_geometry']
            assert i['id']==gid==g['geometry_id'] and g['metal'] in ('La','Dy') and g['metal_index']==0
            assert g['physical_multiplicity']==(1 if g['metal']=='La' else 6)
            for k in ('xyz','key','parameters'): verify(i[k])
            assert i['key']['sha256']==family['base_key']['sha256'] and i['parameters']['sha256']==family['base_parameters']['sha256']
            replay=native.parse(verify(b['native']['analyze_output']),verify(b['native']['gradient_output']),i)
            assert replay==b['native'],'Native AMOEBA saved replay differs: '+gid
            t=[r.split() for r in verify(i['xyz']).read_text().splitlines()[1:] if r.strip()]
            assert len(t)==g['atoms']==i['atoms']==len(i['types'])
            assert [int(v[5]) for v in t]==i['types']
            assert [[int(a)-1 for a in v[6:]] for v in t]==i['bonds_zero_based']
            # Preserve original printed native precision; no geometry regeneration.
            assert max(abs(float(v[k+2])-g['positions_A'][j][k]) for j,v in enumerate(t) for k in range(3))<=5.1e-11
            rr=read(rp);assert rr['status']=='complete' and rr['returncode']==0 and all(rr['checks'].values())
            assert rr['cell_id']=='ecp_pbe0__'+gid and rr['adapter']==cfg['reference_helper']
            p=ref.load_panel(verify(rr['panel']));matches=[c for c in p['cells'] if c['cell_id']==rr['cell_id']];assert len(matches)==1;c=matches[0]
            cg=c['geometry'];assert cg['symbols']==g['symbols'] and cg['positions_A']==g['positions_A']
            assert c['recipe']=='ecp_pbe0' and c['charge']==i['charge']==g['charge']
            assert c['multiplicity']==1 and c['physical_multiplicity']==g['physical_multiplicity']
            assert c['ecp_core']==(46 if g['metal']=='La' else 55)
            for k in ('input','output','engrad','geometry'): verify(rr[k])
            deck=Path(rr['input']['path']).read_text();ranks=int(re.search(r'nprocs\s+(\d+)',deck)[1]);mc=int(re.search(r'%maxcore\s+(\d+)',deck)[1])
            assert deck==ref.deck(p,c,ranks,mc),'Reference input recipe differs'
            again=ref.collect_one(p,c,Path(rr['output']['path']).parent,0)
            assert all(rr[k]==v for k,v in again.items()),'Reference raw replay differs'
            assert finite(rr['forces_kcal_mol_A']) and len(rr['forces_kcal_mol_A'])==g['atoms']
            lines=Path(rr['output']['path']).read_text(errors='replace').splitlines()
            stopping=[{'line':n+1,'text':s.strip()} for n,s in enumerate(lines) if re.search(r'Last (Energy|MAX-Density|RMS-Density|DIIS|Orbital)|SCF CONVERGED|TOTAL SCF ENERGY|Total Energy Change|1-Electron Energy Change',s)]
            row={'status':'complete','geometry':g,'input':i,'baseline_result':bp,'reference_result':rp,
                 'baseline':{'energy':replay['total_energy_kcal_mol'],'forces':replay['force_kcal_mol_A'],'native':replay},
                 'reference':{'energy':rr['energy_hartree']*ref.HA,'energy_hartree':rr['energy_hartree'],'forces':rr['forces_kcal_mol_A'],
                    'native_checks':rr['checks'],'output':rr['output'],'engrad':rr['engrad'],'input':rr['input'],'panel':rr['panel'],
                    'state':{k:c[k] for k in ('charge','multiplicity','physical_multiplicity','ecp_core','explicit_electrons')},
                    'printed_scf_criteria':rr['printed_scf_criteria'],'stopping_evidence':stopping,
                    'physical_state_qualified':rr['physical_state_qualified'],'method_accuracy_qualified':rr['method_accuracy_qualified']}}
            row['bodies']=obj.bodies(row);result['rows'][gid]=row
        for q in obj.WATER+obj.DONOR:
            a,b=[result['rows'][q+'_'+m] for m in ('La','Dy')]
            assert a['geometry']['positions_A']==b['geometry']['positions_A'] and a['bodies']==b['bodies']
            assert a['input']['bonds_zero_based']==b['input']['bonds_zero_based']
            assert a['geometry']['symbols'][1:]==b['geometry']['symbols'][1:]
        result['specification']=cfg['specification'];result['family']=cfg['family'];result['objective']=cfg['objective']
        result['radial_derivative_screen']=obj.radial_screen(result)
        result['radial_screen_all_passed']=all(r['screen_passed'] for r in result['radial_derivative_screen'])
        model=obj.baseline(result)
        result['baseline_objectives']={s:obj.evaluate(result,model,s) for s in ('water','donor')}
        # Meaningful saved-only controls: a self-target has zero residual, and an
        # internal equal/opposite central pair force has zero rigid response.
        refs=obj.references(result)
        assert all(obj.evaluate(result,refs,s)['objective']==0 for s in ('water','donor'))
        test=result['rows']['hydrated_acetamide_central_La'];f=[[0.,0.,0.] for _ in test['geometry']['positions_A']]
        body=next(b for b in test['bodies'] if b['category']=='ligand');a,b=body['indices'][:2]
        delta=[test['geometry']['positions_A'][a][k]-test['geometry']['positions_A'][b][k] for k in range(3)]
        f[a]=delta;f[b]=[-v for v in delta];v,_=obj.rigid_residual(test,f,[[0.,0.,0.] for _ in f]);assert max(abs(x) for x in v)<1e-12
        result.update(status='training_complete',complete=28,saved_controls_passed=True,
                      reference_precision_qualified=False,fit_requires_root_reference_eligibility=True)
    except Exception as exc:
        result['error']=repr(exc)
    finally:
        result['elapsed_seconds']=time.monotonic()-started;result['utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
        save(out/'RESULT.json',result)
    print(json.dumps({'status':result['status'],'result':pin(out/'RESULT.json'),'seconds':result['elapsed_seconds']}))
    return 0 if result['status']=='training_complete' else 1

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();raise SystemExit(prepare(a.config,a.output))
