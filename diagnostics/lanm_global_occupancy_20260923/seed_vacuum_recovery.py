"""Two same-metal ALPB-seeded vacuum retries; original four cells are immutable."""
import argparse,json,re,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
W=ROOT/'workspaces/lanm_global_occupancy_20260923'
sys.path.insert(0,str(W/'scoring_v2/implementation'))
from affordable_common import read_json,record,verify,write_new,cache_key
from affordable_workflow import dry_run,execute
from compact_solvation import completed,diagnostics
from structure_informed_starts import scf_details

def prepare(output,ranks,memory,maxcore):
    original=W/'native_feasibility_retry_v2/manifest.json';m=read_json(original)
    if ranks<1 or maxcore<6000 or 2*ranks*maxcore>memory*1.048576*.76:raise ValueError('invalid resource sizing')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False);tasks=[]
    for metal in ('La','Dy'):
        source=next(t for t in m['tasks'] if t['metal']==metal and t['medium']=='alpb')
        target=next(t for t in m['tasks'] if t['metal']==metal and t['medium']=='vacuum')
        pin=completed(original,source['task_id'])
        if pin is None:raise ValueError('missing converged seed')
        audit=diagnostics(pin,source)
        if audit['charge_sanity_status']!='pass':raise ValueError('seed state failed')
        for k in ('xyz','charge','multiplicity','physical_multiplicity'):
            if source[k] != target[k] and not(k=='xyz' and source[k]['sha256']==target[k]['sha256']):
                raise ValueError('seed/target physical state differs')
        td=Path(target['output_path']).parent;sd=Path(source['output_path']).parent
        if record(td/'endpoint.runtime.xtb.json')['sha256']!=audit['parameter_export']['sha256']:
            raise ValueError('native parameters differ across media')
        d=out/'tasks'/metal;d.mkdir(parents=True)
        shutil.copyfile(verify(target['xyz']),d/'core.xyz')
        old=verify(target['input']).read_text()
        body=re.sub(r'%maxcore \d+',f'%maxcore {maxcore}',old.replace(' NoAutostart',''))
        if 'ALPB' in body or 'TolE 1e-10' not in body:raise ValueError('vacuum recipe changed')
        (d/'endpoint.inp').write_text(body)
        seeds={}
        for suffix in ('gbw','xtbw'):
            original_seed=record(sd/f'endpoint.runtime.{suffix}')
            shutil.copyfile(verify(original_seed),d/f'seed.immutable.{suffix}')
            seeds[suffix]={'source':original_seed,'immutable':record(d/f'seed.immutable.{suffix}')}
        t=dict(target,input=record(d/'endpoint.inp'),xyz=record(d/'core.xyz'),output_path=str(d/'endpoint.out'),
               seed_files=seeds,seed_execution=pin,source_parameter_export=audit['parameter_export'])
        t['scientific_key']=cache_key({'original_key':target['scientific_key'],'input':t['input']['sha256'],
                                     'seeds':seeds,'mpi_ranks':ranks})
        tasks.append(t)
    result=dict(m,protocol_id='nikasha_LanM_native_ALPB_seeded_vacuum_recovery_v1',tasks=tasks,
                execution_resources={'mpi_ranks':ranks,'concurrent_tasks':2,'maxcore_mb_per_rank':maxcore,'scheduler_memory_mib':memory},
                primary_manifest=record(original),preparer=record(Path(__file__)),
                interpretation='initialization recovery only; single seed per metal, no ground-state qualification')
    write_new(out/'manifest.json',result);v=dry_run(out/'manifest.json');write_new(out/'PREFLIGHT.json',v);return v

def run(mp):
    m=read_json(mp);before=[]
    for t in m['tasks']:
        d=Path(t['output_path']).parent
        for suffix,seeds in t['seed_files'].items():
            verify(seeds['source']);shutil.copyfile(verify(seeds['immutable']),d/f'endpoint.runtime.{suffix}')
            before.append({'task':t['task_id'],'suffix':suffix,'active':record(d/f'endpoint.runtime.{suffix}')})
    write_new(Path(mp).parent/'SEEDS_BEFORE.json',before)
    try:return execute(mp)
    finally:collect(mp)

def collect(mp):
    m=read_json(mp);rows=[]
    for t in m['tasks']:
        row={'task_id':t['task_id'],'metal':t['metal'],'medium':'vacuum','status':'unavailable','energy_hartree':None}
        try:
            pin=completed(mp,t['task_id'])
            if pin is None:raise ValueError('no successful receipt')
            text=verify(pin['output']).read_text();audit=diagnostics(pin,t);scf=scf_details(text)
            if 'INITIAL GUESS: XTBRESTART' not in text:raise ValueError('native restart not confirmed')
            if audit['charge_sanity_status']!='pass' or audit['parameter_export']['sha256']!=t['source_parameter_export']['sha256']:
                raise ValueError('native charge/parameters changed')
            if scf['energy']['tolerance']!=1e-10 or not scf['native_mixer_observed']:raise ValueError('SCF policy changed')
            row.update(pin,status='complete',audit=audit,restart_confirmed=True)
        except Exception as exc:row['reason']=str(exc)
        rows.append(row)
    result={'manifest':record(mp),'rows':rows,'complete_cells':sum(r['status']=='complete' for r in rows),
            'cell_denominator':2,'accommodation_status':'unavailable_both_proposals_rejected_geometry',
            'primary_result_changed':False,'ground_state_qualified':False}
    write_new(Path(mp).parent/'COLLECTION.json',result);return result

def main():
    p=argparse.ArgumentParser();p.add_argument('op',choices=['prepare','execute','collect'])
    p.add_argument('--output',type=Path);p.add_argument('--manifest',type=Path)
    p.add_argument('--ranks',type=int);p.add_argument('--memory-mib',type=int);p.add_argument('--maxcore-mb',type=int)
    a=p.parse_args()
    r=prepare(a.output,a.ranks,a.memory_mib,a.maxcore_mb) if a.op=='prepare' else run(a.manifest) if a.op=='execute' else collect(a.manifest)
    print(json.dumps(r,default=str))
if __name__=='__main__':main()
