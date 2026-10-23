#!/usr/bin/env python3
"""Independent saved-only review of the completed 18636 prefix; no native API."""
import hashlib
import importlib.util
import json
import re
import time
from decimal import Decimal
from pathlib import Path

ROOT = Path.cwd()
REL = Path('explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/qualification_v1')
DIAG = ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'/REL
WORK = ROOT/'workspaces/lanm_selectivity_free_energy_20261021'/REL

def pin(p):
    p=Path(p).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def main():
    tick,cpu=time.monotonic(),time.process_time()
    spec=importlib.util.spec_from_file_location('prefix_qualification',DIAG/'RUN.py')
    q=importlib.util.module_from_spec(spec);spec.loader.exec_module(q)
    rp=WORK/'run_v1/execution/RESULT.json';r=json.loads(rp.read_text());cfg=r['config']
    checks=[]
    def ck(ok,name):
        checks.append(name)
        if not ok:raise ValueError(name)
    ck(r['implementation']==pin(DIAG/'RUN.py'),'executed source pin')
    m,h=q.context(cfg,False)
    ck(r['status']=='local_response_numerical_qualification_failed' and r['error']=="RuntimeError('dictionary changed size during iteration')",'original failure preserved')
    ck(r['native_program_invocations']==27 and r['internal_FD_energies_requested']==r['internal_FD_energies_completed']==468,'original ledger 27/468')
    jobs={j['id']:j for j in m['programs']}
    expected={k for k,j in jobs.items() if j['stage'] in ('absent','zero','combined')}
    done={k for k,v in r['rows'].items() if v['status']=='complete'}
    ck(done==expected and len(done)==27,'exact first three completed stages')
    ck(all(v['status']=='not_attempted' for k,v in r['rows'].items() if k not in done),'remaining 30 untouched')
    ck({p.parent.name for p in rp.parent.glob('cells/*/testgrad.out')}==done,'raw output denominator')
    rows={};fd_count=0;identity_count=0;induction_count=0;max_res=0.;max_fd=Decimal(0)
    for cid in sorted(done):
        rec=r['rows'][cid];job=jobs[cid];row=m['rows'][job['cell']]
        d=rp.parent/'cells'/cid
        ck(json.loads((d/'RESULT.json').read_text())==rec,cid+' independent per-cell receipt')
        ck(rec['job']==job and len(rec['receipts'])==1,cid+' exact job/single receipt')
        for k,name in q.FILES:
            ck(pin(d/name)==rec['runtime_input'][k],cid+' runtime '+k)
            ck((d/name).read_bytes()==q.verify(job['input'][k]).read_bytes(),cid+' prepared byte equality '+k)
        receipt=rec['receipts'][0];out=d/'testgrad.out';text=out.read_text()
        argv=[cfg['testgrad']['path'],'system.xyz','Y']+(['Y','0.00001'] if job['fd'] else ['N'])
        ck(receipt['returncode']==0 and receipt['argv']==argv and receipt['cwd']==str(d),cid+' exact native invocation')
        ck(receipt['output']==rec['actual_output']==pin(out),cid+' raw output pins')
        ck(receipt['threads']=={k:'1' for k in q.THREADS} and receipt['timeout_seconds']==300,cid+' single-thread execution')
        ck(not re.search(r'\b(?:NaN|Infinity)\b',text,re.I),cid+' finite raw tokens')
        rb=q.response_readback(out,job,row);ind=h['induction'].induction(out)
        grad=h['gradient'].gradient(out,row['atoms'],job['fd'],h['native'])
        ck(rb==rec['response_readback'] and ind==rec['induction'] and grad==rec['value'],cid+' exact saved parser replay')
        section=text.rsplit('Cartesian Gradient Breakdown over Individual Atoms :',1)[-1]
        tokens=re.findall(r'^\s*(Anlyt|Numer)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s*$',section,re.M)
        parsed={kind:{} for kind in ('Anlyt','Numer')}
        for kind,index,*v in tokens:
            ck(int(index) not in parsed[kind],cid+' raw gradient unique '+kind+index)
            parsed[kind][int(index)]=[Decimal(x.replace('D','E')) for x in v]
        ck(set(parsed['Anlyt'])==set(range(1,row['atoms']+1)),cid+' independent gradient denominator')
        ck([[float(v) for v in parsed['Anlyt'][i]] for i in range(1,row['atoms']+1)]==grad['gradient'],cid+' independent raw analytic gradient')
        error=None
        if job['fd']:
            ck(set(parsed['Numer'])==set(parsed['Anlyt']),cid+' independent numerical denominator')
            error=max(abs(a-b) for i in parsed['Anlyt'] for a,b in zip(parsed['Anlyt'][i],parsed['Numer'][i]))
            ck(error<=Decimal('0.001') and grad['fd_passed'],cid+' independent Cartesian FD threshold')
            fd_count+=6*row['atoms'];max_fd=max(max_fd,error)
        else:ck(not parsed['Numer'],cid+' no undeclared finite differences')
        if 'identity' in rec:
            old=row['baseline_gradient'];identity_count+=1
            ck(grad['energy']==old['energy'] and grad['forces']==old['forces'],cid+' exact original E/F identity')
        induction_count+=ind['solutions'];max_res=max(max_res,ind['maximum_printed_RMS_residual_Debye'])
        rows[cid]={'result':pin(d/'RESULT.json'),'output':pin(out),'native_returncode':0,'FD_max_decimal':str(error) if error is not None else None,'induction':ind,'identity_exact':True if 'identity' in rec else None}
    ck(fd_count==468 and identity_count==21,'independent FD/identity totals')
    review={'status':'passed_completed_native_prefix_review','result':pin(rp),'implementation':r['implementation'],'reviewer':pin(__file__),'config':pin(WORK/'run_v1/CONFIG_RUN.json'),'source_manifest':cfg['source_manifest'],'build_receipt':cfg['build_receipt'],'analyze':cfg['analyze'],'testgrad':cfg['testgrad'],'completed_programs':27,'untouched_programs':30,'completed_internal_FD_energies':fd_count,'exact_original_EF_identities':identity_count,'maximum_FD_error_kcal_mol_A_decimal':str(max_fd),'induction_summaries':induction_count,'maximum_printed_induction_residual_Debye':max_res,'checks':checks,'rows':rows,'native_derivative_failure_found':False,'full_qualification_passed':False,'remaining':'Drive/stiffness axis gradients, directional paths, operator exports and static-field checks remain unexecuted. Original aggregate persistence failure is preserved.','new_native_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    dest=WORK/'independent_partial_review_v1';dest.mkdir(exist_ok=False)
    (dest/'REVIEW.json').write_text(json.dumps(review,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'review':pin(dest/'REVIEW.json'),'checks':len(checks),'max_FD_error':str(max_fd),'induction_summaries':induction_count,'max_induction':max_res}))

if __name__=='__main__':main()
