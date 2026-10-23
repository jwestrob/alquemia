#!/usr/bin/env python3
"""Saved-only admission of original DLPNO scout outputs. Never runs an engine."""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import re
import time

HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('original_local_cc',HERE.parent/'REFERENCE.py')
old=importlib.util.module_from_spec(s);s.loader.exec_module(old)
b=old.b
ORIGINAL_SHA='ea95da207a0a205f02922e1d7fe144e2f4efa0cb9636591c98d8c885934b7748'
CHANGED={'local_'+k for k in old.TIERS['tight']}|{'local_MP2_converged'}

def numbers(section,key,separator):
    tokens=re.findall(r'^\s*'+re.escape(key)+r'\s*'+separator+r'\s*(\S+)\s*$',section,re.M)
    return [old.cp.number(x) if re.fullmatch(old.F,x) else None for x in tokens]

def exact(values,value):
    return len(values)==1 and values[0] is not None and math.isfinite(values[0]) and abs(values[0]-value)<=abs(value)*1e-9

def stage_evidence(text,tier):
    native=text.split('****END OF INPUT****')[-1]
    headers=list(re.finditer(r'^DLPNO SETTINGS \(2015 fully linear scaling implementation\)\s*$',native,re.M))
    primary=native[headers[0].end():].split('INITIAL GUESS AND PNO GENERATION',1)[0] if len(headers)==1 else ''
    target=old.TIERS[tier]
    main={k:numbers(primary,k,':') for k in target}
    scales={k:numbers(primary,k,':') for k in ('TScaleDOMP2PreScr','TScaleMKNMP2PreScr','TScalePNOMP2PreScr')}
    parts=native.split('INITIAL GUESS AND PNO GENERATION')[1:]
    stages=[]
    for part in parts:
        kinds=re.findall(r'^Type of local MP2 treatment\s*\.{2,}\s*([^\n]+)',part,re.M)
        kind=kinds[0].strip() if len(kinds)==1 else None
        vals={k:numbers(part,k,'=') for k in target}
        item=dict(kind=kind,thresholds=vals)
        if kind=='full local MP2':
            et=numbers(part,'Energy convergence tolerance',r'\.{2,}')
            rt=numbers(part,'Residual convergence tolerance',r'\.{2,}')
            lines=re.findall(r'^\s*LMP2-Iter=[^\n]*',part,re.M)
            rows=[];unsupported=[]
            for line in lines:
                match=re.fullmatch(r'\s*LMP2-Iter=\s*(\d+):\s*ELMP2=\s*('+old.F+r')\s+DE=\s*('+old.F+r')\s+RMAX=\s*('+old.F+r')\s*(=> CONVERGED)?\s*',line)
                if match:
                    rows.append(dict(iteration=int(match[1]),energy=old.cp.number(match[2]),delta=old.cp.number(match[3]),residual=old.cp.number(match[4]),converged=bool(match[5])))
                else:unsupported.append(line)
            terminal=rows[-1] if rows else None
            passed=(exact(et,1e-7) and exact(rt,5e-7) and bool(rows) and not unsupported
                and [r['iteration'] for r in rows]==list(range(len(rows)))
                and all(math.isfinite(r[k]) for r in rows for k in ('energy','delta','residual'))
                and sum(r['converged'] for r in rows)==1 and terminal['converged']
                and abs(terminal['delta'])<=et[0] and 0<=terminal['residual']<=rt[0])
            item.update(energy_tolerance=et,residual_tolerance=rt,iterations=rows,unsupported_lines=unsupported,terminal=terminal,converged=bool(passed))
        stages.append(item)
    # This is the observed native two-stage TightPNO algorithm. Unknown grammar
    # remains unavailable; the coarse stage is not relabeled as the final stage.
    ordering=[x['kind'] for x in stages]==['semi-local MP2','full local MP2']
    scale_ok=exact(scales['TScaleDOMP2PreScr'],2.) and exact(scales['TScaleMKNMP2PreScr'],10.) and exact(scales['TScalePNOMP2PreScr'],1.)
    factors={'TCutPairs':1.,'TCutPNO':1.,'TCutMKN':10.,'TCutDO':2.}
    maps=[dict(stage=kind,printed=section,thresholds={k:numbers(section,k,r'\.{2,}') for k in target})
        for kind,section in re.findall(r'^Thresholds for map construction and integral transformation for ([^\n]+):\n((?:[ \t]+TCut[^\n]+\n)+)',native,re.M)]
    fine=[x for x in maps if x['stage']=='fine MP2 and CCSD(T) calculation']
    checks={}
    for key,value in target.items():
        checks['local_'+key]=(len(headers)==1 and exact(main[key],value) and ordering and scale_ok
            and exact(stages[0]['thresholds'][key],value*factors[key])
            and exact(stages[1]['thresholds'][key],value)
            and (key=='TCutPNO' or (bool(fine) and all(exact(x['thresholds'][key],value) for x in fine))))
    checks['local_MP2_converged']=ordering and stages[1]['converged']
    return dict(primary_thresholds=main,prescreen_scales=scales,stages=stages,map_threshold_reports=maps,
        other_threshold_lines=[dict(line=i+1,text=x) for i,x in enumerate(native.splitlines()) if re.search(r'TCut|TScale|ScaleTCut',x)],checks=checks)

def parse(text,c,returncode):
    q=old.parse(text,c,returncode)
    original_checks=dict(q['checks']);stages=stage_evidence(text,c['tier'])
    q['checks'].update(stages['checks'])
    q.update(original_parser_checks=original_checks,stage_specific_local_evidence=stages)
    q['status']='complete' if all(q['checks'].values()) else 'failed_admission'
    return q

def load(panel):
    b.check(b.pin(old.__file__)['sha256']==ORIGINAL_SHA,'Frozen producer source differs')
    return old.load(panel)

def admit_original(path,c,p,panel):
    path=Path(path);r=json.loads(path.read_text())
    b.check(r['status'] in ('complete','failed_admission') and r['returncode']==0,'No normal native result to recover')
    b.check(r['cell_id']==c['cell_id'] and r['cell']==c and r['panel']==b.pin(panel) and r['adapter']==p['adapter'],'Original cell/source differs')
    st=old.read(r['start']);rel=old.read(r['release']);old.release_contract(p,panel,rel)
    b.check(rel['stage']=='controls_scout' and c['cell_id'] in old.CONTROL_SCOUTS,'This recovery admits only original four scouts')
    b.check(r['start']==b.pin(path.parent/'START.json') and st['cell']==c and st['panel']==r['panel'] and st['release']==r['release'] and st['input']==r['input'],'START lineage differs')
    b.check(c['cell_id'] in rel['released_cell_ids'] and path.resolve()==Path(rel['execution_root']).resolve()/'tasks'/c['cell_id']/'execution/RESULT.json','Original one-attempt path differs')
    b.check(st['ranks']==4 and st['memory_MiB']>0 and st['maxcore_MiB']==math.floor(.75*st['memory_MiB']/4),'Allocation differs')
    b.check(st['command']==[p['orca']['path'],'endpoint.inp','--bind-to none'] and r['native_program_invocations_attempted']==1,'Original invocation differs')
    b.check(b.verify(r['input']).read_text()==old.deck(p,c,4,st['maxcore_MiB']),'Executed deck differs');b.verify(r['stderr'])
    raw=b.verify(r['output']).read_text();previous=old.parse(raw,c,r['returncode'])
    # Native basis_groups uses integer keys; the saved JSON necessarily uses
    # string keys. Compare the complete parser record in its serialized form.
    b.check(all(r[k]==v for k,v in json.loads(json.dumps(previous)).items()),'Original failed result/raw replay differs')
    b.check(all(v or k in CHANGED for k,v in previous['checks'].items()),'An unrelated native gate failed')
    q=parse(raw,c,r['returncode']);b.check(q['status']=='complete','Corrected native admission failed: '+','.join(k for k,v in q['checks'].items() if not v))
    return dict(q,original_result=b.pin(path),original_status=r['status'],original_adapter=r['adapter'],parser=b.pin(__file__),
        cell_id=c['cell_id'],cell=c,panel=b.pin(panel),output=r['output'],input=r['input'],stderr=r['stderr'],start=r['start'],release=r['release'],
        elapsed_seconds=r['elapsed_seconds'],allocated_cpu_seconds=r['allocated_cpu_seconds'],job_id=r['job_id'],new_native_calls=0)

def recover(panel,result,out):
    tick,cpu=time.monotonic(),time.process_time();p=load(panel);r=json.loads(Path(result).read_text())
    c=next(c for c in p['cells'] if c['cell_id']==r['cell_id']);q=admit_original(result,c,p,panel)
    out.mkdir(parents=True,exist_ok=False)
    q.update(status='saved_native_recovery_complete',implementation=b.pin(__file__),declaration=b.pin(HERE/'DECLARATION.md'),
        wall_seconds=time.monotonic()-tick,parent_cpu_seconds=time.process_time()-cpu,local_threads=1,physical_accuracy_qualified=False,
        all_eighteen_method_controls_qualified=False)
    pin=b.write(out/'RESULT.json',q)
    print(json.dumps(dict(result=pin,cell_id=q['cell_id'],status=q['status'],new_native_calls=0)))
    return q

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--panel',required=True);ap.add_argument('--result',required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();recover(a.panel,a.result,a.out)
if __name__=='__main__':main()
