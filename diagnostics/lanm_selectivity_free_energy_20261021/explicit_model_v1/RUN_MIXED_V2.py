#!/usr/bin/env python3
"""Frozen hydrated donor geometries, unchanged native Tinker Hamiltonian."""
import argparse, hashlib, importlib.util, json, os, shutil, socket, time
from pathlib import Path

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def verify(p):
    assert pin(p['path'])==p;return Path(p['path'])

def read(p):return json.loads(verify(p).read_text())
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')

def binding(c):
    assert c['implementation']==pin(__file__);verify(c['specification'])
    failed=read(c['failed_preparation']);assert failed['status']=='unavailable' and failed['new_native_program_invocations']==0
    sp=importlib.util.spec_from_file_location('original_tinker_panel',verify(c['helper']));h=importlib.util.module_from_spec(sp);sp.loader.exec_module(h)
    assert all(os.environ.get(k)=='1' for k in h.THREAD_KEYS)
    parent=read(c['parent_result']);original=read(c['original_manifest']);panel=read(c['panel'])
    assert parent['status']=='panel_complete' and parent['complete_rows']==32 and parent['manifest']==c['original_manifest']
    assert parent['implementation']==c['helper'] and original['transfer_qualified'] is False
    assert panel['status']=='prepared' and len(panel['geometries'])==18 and len(panel['reused_controls'])==2
    assert panel['contact_failures']==[], 'Prepared panel must have no contact failures'
    for k in ('specification','implementation','source_panel','helper'):verify(panel[k])
    for p in c['executables'].values():verify(p)
    pc=read(parent['config']);assert c['executables']==pc['executables']
    for row in original['rows']:
        for k in ('xyz','key','parameters'):verify(row[k])
    for g in panel['geometries']+panel['reused_controls']:verify(g['xyz'])
    return h,parent,original,panel

def prepare(c,out,h,parent,original,panel):
    template=original['rows'][0];rows=[]
    for g in panel['geometries']:
        assert g['metal_index']==0 and g['metal'] in ('La','Dy')
        n=g['atoms'];assert n==len(g['positions_A'])==len(g['symbols'])==len(g['source_tinker_types'])==len(g['source_bonds_zero_based'])
        types=list(g['source_tinker_types']);assert types[0]==500;types[0]=500 if g['metal']=='La' else 501
        bonds=g['source_bonds_zero_based'];assert bonds[0]==[]
        assert all(i in bonds[j] for i,b in enumerate(bonds) for j in b)
        assert g['charge']==(3 if g['donor_family']=='acetamide' else 2)
        d=out/g['geometry_id'];d.mkdir()
        text=f'{n}  {g["geometry_id"]}\n'+''.join(f'{i+1:5d} {s:3s} '+''.join(f'{v:20.12f}' for v in xyz)+f' {t:5d}'+''.join(f' {j+1:5d}' for j in b)+'\n' for i,(s,xyz,t,b) in enumerate(zip(g['symbols'],g['positions_A'],types,bonds)))
        (d/'system.xyz').write_text(text)
        for k,name in [('key','system.key'),('parameters','amoeba09.prm')]:shutil.copyfile(verify(template[k]),d/name)
        rows.append({'id':g['geometry_id'],'atoms':n,'metal':g['metal'],'charge':g['charge'],'derivative_test':False,
                     'source_geometry':g,'types':types,'bonds_zero_based':bonds,'xyz':pin(d/'system.xyz'),
                     'key':pin(d/'system.key'),'parameters':pin(d/'amoeba09.prm'),'reuse':False})
    for g in panel['reused_controls']:
        r=next(r for r in original['rows'] if r['id']==g['geometry_id']);old=r['source_geometry']
        assert all(old[k]==g[k] for k in ('symbols','positions_A','atoms','charge','metal','water_count'))
        rows.append({**r,'reuse':True})
    assert len(rows)==20 and len({r['id'] for r in rows})==20
    m={'protocol':'explicit_model_v1_hydrated_mixed','status':'prepared','config':pin(c['_path']),
       'implementation':pin(__file__),'parent_result':c['parent_result'],'panel':c['panel'],
       'rows':rows,'new_geometry_slots':18,'reused_aquo_slots':2,'energy_calls':0,'force_calls':0,
       'ALPB':False,'transfer_qualified':False}
    save(out/'MANIFEST.json',m);return {'status':'prepared','manifest':pin(out/'MANIFEST.json'),'new_native_program_invocations':0}

def run(c,out,h,parent,original,panel,result):
    assert os.environ.get('SLURM_JOB_ID'),'Root-owned Slurm run required'
    m=read(c['manifest']);assert m['status']=='prepared' and m['parent_result']==c['parent_result'] and m['panel']==c['panel'] and m['implementation']==pin(__file__)
    prepared_config=read(m['config']);assert all(prepared_config[k]==c[k] for k in prepared_config if k not in ('manifest','_path'))
    rows=m['rows'];assert len(rows)==20 and sum(not r['reuse'] for r in rows)==18
    result['rows']=[{'id':r['id'],'input':r,'status':'unavailable','native':None,'receipts':[]} for r in rows]
    stopped=False
    for r,t in zip(rows,result['rows']):
        if stopped:continue
        try:
            for k in ('xyz','key','parameters'):verify(r[k])
            if r['reuse']:
                old=next(p for p in parent['rows'] if p['id']==r['id']);assert old['status']=='complete'
                native=h.parse(verify(old['native']['analyze_output']),verify(old['native']['gradient_output']),r)
                assert native==old['native'];t.update(native=native,reused_from=c['parent_result'])
            else:
                d=out/r['id'];d.mkdir()
                for k,name in [('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')]:shutil.copyfile(verify(r[k]),d/name)
                for name,args in [('analyze',[c['executables']['analyze']['path'],'system.xyz','E']),('gradient',[c['executables']['testgrad']['path'],'system.xyz','Y','N'])]:
                    result['new_native_program_invocations']+=1;h.invoke(args,d,d/(name+'.out'),t['receipts'])
                t['native']=h.parse(d/'analyze.out',d/'gradient.out',r)
                g=r['source_geometry'];grad=t['native']['gradient_kcal_mol_A']
                t['ligand_radial_derivative_kcal_mol_A']=sum(grad[i][k]*g['radial_direction'][k] for i in g['ligand_indices'] for k in range(3))
            t['status']='complete'
        except Exception as e:t['error']=repr(e);stopped=True
        save(out/'RESULT.json',result)
    result['complete_rows']=sum(t['status']=='complete' for t in result['rows'])
    result['status']='complete' if result['complete_rows']==20 else 'incomplete_or_failed'
    by={t['id']:t for t in result['rows']};contrasts=[]
    for family in ('acetate_mono','acetate_bi','acetamide'):
      for displacement in ('minus','central','plus'):
        ids=[f'hydrated_{family}_{displacement}_{m}' for m in ('La','Dy')];p=[by[i] for i in ids];a=[by[f'aquo8_central_{m}'] for m in ('La','Dy')]
        row={'family':family,'distortion':displacement,'ids':ids,'status':'unavailable'}
        if all(t['status']=='complete' for t in p+a):
            row.update(status='complete',Dy_minus_La_ligand_replacement_kcal_mol=(p[1]['native']['total_energy_kcal_mol']-a[1]['native']['total_energy_kcal_mol'])-(p[0]['native']['total_energy_kcal_mol']-a[0]['native']['total_energy_kcal_mol']),
                Dy_minus_La_radial_derivative_kcal_mol_A=p[1]['ligand_radial_derivative_kcal_mol_A']-p[0]['ligand_radial_derivative_kcal_mol_A'])
        contrasts.append(row)
    result.update(contrasts=contrasts,absolute_replacement_energies=None,free_water_evaluated=False,
       interpretation='Fixed balanced metal-differential; free ligand and water cancel. Not free energy or validated selectivity.',transfer_qualified=False,ALPB=False)

def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=('prepare','run'),required=True);p.add_argument('--config',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    c=json.loads(Path(a.config).read_text());c['_path']=a.config
    r={'status':'unavailable','stage':a.stage,'config':pin(a.config),'implementation':pin(__file__),'host':socket.gethostname(),
       'slurm_job_id':os.environ.get('SLURM_JOB_ID'),'new_native_program_invocations':0,'rows':[]}
    try:
        h,parent,original,panel=binding(c)
        if a.stage=='prepare':r.update(prepare(c,out,h,parent,original,panel))
        else:run(c,out,h,parent,original,panel,r)
    except Exception as e:r['error']=repr(e)
    r['elapsed_seconds']=time.monotonic()-start;save(out/'RESULT.json',r);print(json.dumps({'status':r['status'],'result':pin(out/'RESULT.json')}))
    return 0 if r['status'] in ('prepared','complete') else 1

if __name__=='__main__':raise SystemExit(main())
