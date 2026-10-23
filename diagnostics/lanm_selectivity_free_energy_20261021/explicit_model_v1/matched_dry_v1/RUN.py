#!/usr/bin/env python3
"""Six exact water-deletion cells, unchanged AMOEBA; no engines in prepare."""
import argparse, hashlib, importlib.util, json, os, shutil, socket, time
from pathlib import Path

PROTOCOL='explicit_model_v1_matched_dry_v1'
IDS=[f'matched_dry_{d}_central_{m}' for d in ('acetate_mono','acetate_bi','acetamide') for m in ('La','Dy')]
FIXED={'panel':'0352a54293309d02e84be9a99fd34835f5f0ba2641eef52be078b5f7b295d004',
       'hydrated_result':'2d55465328abcd60dcf0ca5b14761d3f4c96314310c7ec0381a2696d8c3f8674',
       'helper':'4189602145853bc39d8dfaf30915b2a6b1d4632984e0c50fc0bf321644fa306c'}
def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):assert pin(p['path'])==p,p;return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def save(p,r):Path(p).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
def xyz(p):
    lines=verify(p).read_text().splitlines();a=[s.split() for s in lines[2:] if s.strip()]
    assert len(a)==int(lines[0]) and all(len(s)==4 for s in a)
    return a
def binding(c):
    assert c['protocol']==PROTOCOL and c['implementation']==pin(__file__);verify(c['specification'])
    for k,h in FIXED.items():assert c[k]['sha256']==h;verify(c[k])
    spec=importlib.util.spec_from_file_location('original_tinker_panel',verify(c['helper']));h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    assert all(os.environ.get(k)=='1' for k in h.THREAD_KEYS)
    panel=read(c['panel']);parent=read(c['hydrated_result']);pc=read(parent['config'])
    assert panel['status']=='prepared_no_molecular_calls' and panel['new_molecular_calls']==0
    assert [g['geometry_id'] for g in panel['geometries']]==IDS
    assert len(panel['checks'])==6 and all(t['surviving_line_equality'] for t in panel['checks'])
    assert parent['status']=='complete' and parent['complete_rows']==20 and len(parent['rows'])==20
    verify(parent['implementation']);verify(panel['implementation']);verify(panel['specification'])
    assert pc['executables']==c['executables'] and pc['helper']==c['helper']
    for p in c['executables'].values():verify(p)
    parents={r['id']:r for r in parent['rows']};rows=[]
    for g in panel['geometries']:
        old=parents[g['source_geometry_id']];assert old['status']=='complete';inp=old['input'];wet=inp['source_geometry']
        for k in ('xyz','key','parameters'):verify(inp[k])
        for k in ('analyze_output','gradient_output'):verify(old['native'][k])
        assert g['source_geometry_id']==g['geometry_id'].replace('matched_dry_','hydrated_')
        kept=g['retained_source_indices'];deleted=g['deleted_water_source_indices'];n=g['atoms']
        assert sorted(kept+deleted)==list(range(wet['atoms'])) and len(set(kept+deleted))==wet['atoms']
        assert kept==[wet['metal_index']]+wet['ligand_indices'] and g['metal_index']==0
        assert len(kept)==n and g['water_count']==0 and g['removed_water_count']==wet['water_count']
        assert sum(wet['symbols'][i]=='O' for i in deleted)==g['removed_water_count']
        assert sum(wet['symbols'][i]=='H' for i in deleted)==2*g['removed_water_count'] and len(deleted)==3*g['removed_water_count']
        a=xyz(g['xyz']);orig=xyz(g['source_xyz']);assert a==[orig[i] for i in kept]
        assert [v[0] for v in a]==g['symbols']==[wet['symbols'][i] for i in kept]
        assert [[float(v) for v in r[1:]] for r in a]==g['positions_A']==[wet['positions_A'][i] for i in kept]
        assert all(g[k]==wet[k] for k in ('metal','charge','physical_multiplicity','radial_direction'))
        assert g['charge']==(3 if g['donor_family']=='acetamide' else 2)
        assert g['physical_multiplicity']==(1 if g['metal']=='La' else 6)
        types=[inp['types'][i] for i in kept];expected=list(g['source_tinker_types']);expected[0]=500 if g['metal']=='La' else 501;assert types==expected
        remap={j:i for i,j in enumerate(kept)};bonds=[[remap[j] for j in inp['bonds_zero_based'][i]] for i in kept]
        assert bonds==g['source_bonds_zero_based'] and bonds[0]==[] and all(i in bonds[j] for i,b in enumerate(bonds) for j in b)
        rows.append({'id':g['geometry_id'],'atoms':n,'metal':g['metal'],'charge':g['charge'],'derivative_test':False,
                     'source_geometry':g,'types':types,'bonds_zero_based':bonds,'hydrated_parent_id':old['id'],
                     'parent_key':inp['key'],'parent_parameters':inp['parameters']})
    for i in range(0,6,2):
        a,b=rows[i:i+2];assert a['source_geometry']['positions_A']==b['source_geometry']['positions_A'];assert a['types'][1:]==b['types'][1:] and a['bonds_zero_based']==b['bonds_zero_based']
    return h,rows
def validate_files(row):
    for k in ('xyz','key','parameters'):verify(row[k])
    assert row['key']['sha256']==row['parent_key']['sha256'] and row['parameters']['sha256']==row['parent_parameters']['sha256']
    lines=verify(row['xyz']).read_text().splitlines();assert int(lines[0].split()[0])==row['atoms']
    a=[s.split() for s in lines[1:] if s.strip()];g=row['source_geometry'];assert len(a)==row['atoms']
    for i,s in enumerate(a):
        assert int(s[0])==i+1 and s[1]==g['symbols'][i]
        assert [float(v) for v in s[2:5]]==g['positions_A'][i] and int(s[5])==row['types'][i]
        assert [int(v)-1 for v in s[6:]]==row['bonds_zero_based'][i]
def prepare(c,out,rows):
    result=[]
    for row in rows:
        d=out/row['id'];d.mkdir();g=row['source_geometry'];a=xyz(g['xyz'])
        content=f'{row["atoms"]}  {row["id"]}\n'
        for i,(s,t,b) in enumerate(zip(a,row['types'],row['bonds_zero_based'])):
            content+=f'{i+1:5d} {s[0]:3s} '+' '.join(s[1:])+f' {t:5d}'+''.join(f' {j+1:5d}' for j in b)+'\n'
        (d/'system.xyz').write_text(content)
        shutil.copyfile(verify(row['parent_key']),d/'system.key');shutil.copyfile(verify(row['parent_parameters']),d/'amoeba09.prm')
        r={**row,'xyz':pin(d/'system.xyz'),'key':pin(d/'system.key'),'parameters':pin(d/'amoeba09.prm')};validate_files(r);result.append(r)
    m={'protocol':PROTOCOL,'status':'prepared','implementation':pin(__file__),'config':pin(c['_path']),
       'specification':c['specification'],'panel':c['panel'],'hydrated_result':c['hydrated_result'],'rows':result,
       'new_native_program_invocations':0,'energy_calls':0,'force_calls':0,'ALPB':False,'transfer_qualified':False}
    save(out/'MANIFEST.json',m);return {'status':'prepared','manifest':pin(out/'MANIFEST.json'),'prepared_rows':6}
def run(c,out,h,expected,r):
    assert os.environ.get('SLURM_JOB_ID'),'Root-owned Slurm run required'
    prep=read(c['preparation']);assert prep['status']=='prepared' and prep['new_native_program_invocations']==0 and prep['manifest']==c['manifest']
    m=read(c['manifest']);pc=read(m['config']);assert prep['config']==m['config']
    assert {k:v for k,v in pc.items() if k not in ('manifest','preparation')}=={k:v for k,v in c.items() if k not in ('manifest','preparation','_path')}
    assert m['implementation']==pin(__file__) and m['panel']==c['panel'] and m['hydrated_result']==c['hydrated_result'] and m['status']=='prepared'
    assert len(m['rows'])==6
    for actual,e in zip(m['rows'],expected):
        assert {k:v for k,v in actual.items() if k not in ('xyz','key','parameters')}==e;validate_files(actual)
    r['manifest']=c['manifest'];r['rows']=[{'id':row['id'],'input':row,'status':'unavailable','native':None,'receipts':[]} for row in m['rows']]
    for t in r['rows']:
        row=t['input']
        try:
            d=out/row['id'];d.mkdir()
            for k,name in (('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')):
                shutil.copyfile(verify(row[k]),d/name);assert pin(d/name)['sha256']==row[k]['sha256']
            for name,argv in [('analyze',[c['executables']['analyze']['path'],'system.xyz','E']),('gradient',[c['executables']['testgrad']['path'],'system.xyz','Y','N'])]:
                r['new_native_program_invocations']+=1;assert r['new_native_program_invocations']<=12;h.invoke(argv,d,d/(name+'.out'),t['receipts'])
            t['native']=h.parse(d/'analyze.out',d/'gradient.out',row);g=row['source_geometry']
            t['ligand_radial_derivative_kcal_mol_A']=sum(t['native']['gradient_kcal_mol_A'][i][k]*g['radial_direction'][k] for i in g['ligand_indices'] for k in range(3))
            t['status']='complete'
        except Exception as e:t['error']=repr(e);break
        finally:save(out/'RESULT.json',r)
    r['complete_rows']=sum(t['status']=='complete' for t in r['rows']);r['status']='complete' if r['complete_rows']==6 else 'incomplete_or_failed'
def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=('prepare','run'),required=True);p.add_argument('--config',required=True);p.add_argument('--output',required=True);a=p.parse_args();start=time.monotonic()
    out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=False);c=json.loads(Path(a.config).read_text());c['_path']=a.config
    r={'protocol':PROTOCOL,'status':'unavailable','stage':a.stage,'config':pin(a.config),'implementation':pin(__file__),'host':socket.gethostname(),
       'slurm_job_id':os.environ.get('SLURM_JOB_ID'),'new_native_program_invocations':0,'rows':[{'id':i,'status':'unavailable','input':None,'native':None,'receipts':[]} for i in IDS],
       'ALPB':False,'transfer_qualified':False}
    try:
        h,rows=binding(c)
        if a.stage=='prepare':r.update(prepare(c,out,rows))
        else:run(c,out,h,rows,r)
    except Exception as e:r['error']=repr(e)
    r['elapsed_seconds']=time.monotonic()-start;save(out/'RESULT.json',r);print(json.dumps({'status':r['status'],'result':pin(out/'RESULT.json'),'error':r.get('error')}))
    return 0 if r['status'] in ('prepared','complete') else 1
if __name__=='__main__':raise SystemExit(main())
