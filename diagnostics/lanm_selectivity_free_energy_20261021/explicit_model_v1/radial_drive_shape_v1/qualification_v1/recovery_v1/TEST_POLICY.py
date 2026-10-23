#!/usr/bin/env python3
"""Tiny manufactured parser and mocked staged-execution controls; no native/solves."""
import argparse,copy,importlib.util,json,os,resource,subprocess,tempfile,time
from pathlib import Path
import RUN as R
P=R.P;F=R.F

def main(cfg,out):
    tick,cpu=time.monotonic(),time.process_time();checks=[]
    def check(name,ok):P.require(ok,name);checks.append({'name':name,'passed':True})
    def rejects(name,fn):
        try:fn()
        except (ValueError,AssertionError):check(name,True)
        else:check(name,False)
    # Source/header tests only: the numerical qualifier is never invoked.
    inp=P.read(cfg['inputs']);matrix=P.module(inp['matrix_helper'],'radial_recovery_control_matrix');units={'electric':332.0637133,'debye_per_eA':4.803206802,'elefield':1439.96455,'MV_cm_per_field_au':5142.206758812493}
    with tempfile.TemporaryDirectory(prefix='radial_recovery_controls_') as td:
        d=Path(td);(d/'system.xyz').write_text('4 artificial\n1 M 0 0 0 500\n2 O 0 1 0 36\n3 H 1 1 0 37\n4 H -1 1 0 37\n');(d/'amoeba09.prm').write_text('polarize 500 1.0\npolarize 36 1.0\npolarize 37 1.0\n')
        for sign in (-1,1):
            (d/'system.key').write_text('EXTERNAL-FIELD 0 0 '+format(sign*.001*units['MV_cm_per_field_au'],'.17g')+'\n')
            pins={k:P.pin(d/name) for k,name in P.FILES};job={'id':'fields__water_central_La__analyze__'+str(sign),'stage':'fields','program':'analyze','cell':'water_central_La','sign':sign,'field_au':sign*.001,'fd':False,'input':pins};row={'atoms':4,'types':[500,36,37,37],'native_names':['M','O','H','H'],'input':pins}
            z=float((d/'system.key').read_text().split()[-1])/units['elefield'];tokens=['EXTERNAL','T','0','0',format(z,'.17e')]
            check('valid_signed_header_'+str(sign),F.header(tokens,row,units,job,matrix)['maximum_error']==0.)
            for name,token in [('false_flag',['EXTERNAL','F','0','0',str(z)]),('wrong_sign',['EXTERNAL','T','0','0',str(-z)]),('zero',['EXTERNAL','T','0','0','0']),('transverse',['EXTERNAL','T','1e-15','0',str(z)]),('nonfinite',['EXTERNAL','T','0','0','nan'])]:rejects(name+'_'+str(sign),lambda token=token:F.header(token,row,units,job,matrix))
            bad=copy.deepcopy(job);bad['cell']='aquo9_central_La';rejects('undeclared_cell_'+str(sign),lambda:F.header(tokens,row,units,bad,matrix))
            lines=['TINKER_RESPONSE_MATRIX_V1','COUNTS 1 4 4 12','UNITS A=Angstrom^-3 RHS=e/Angstrom^2 U=eA','CONSTANTS 332.0637133 1 4.803206802',' '.join(tokens),'CHECKS 0 0 0 0 0 0','RESTORED T']
            for i,(name,typ,xyz) in enumerate(zip(row['native_names'],row['types'],[(0,0,0),(0,1,0),(1,1,0),(-1,1,0)]),1):lines.append('SITE '+str(i)+' '+str(i)+' '+str(typ)+' '+name+' 1 '+' '.join(map(str,xyz)))
            for i in range(12):lines.extend(['VECTOR '+str(i+1)+' 0 0 0 0','A '+str(i+1)+' '+' '.join('1' if i==j else '0' for j in range(12))])
            lines.append('END_TINKER_RESPONSE_MATRIX_V1');p=d/'artificial.dat';p.write_text('\n'.join(lines)+'\n')
            proxy=F.scoped(matrix,job);parsed=proxy.parse_export(p,row,units)
            check('full_manufactured_field_parser_'+str(sign),parsed['sequence']==1 and len(parsed['A'])==12 and parsed['external_field']['maximum_error']==0.)
            check('numerical_qualifier_code_unchanged_'+str(sign),proxy.qualify_exports.__code__ is matrix.qualify_exports.__code__)
            rejects('old_zero_parser_still_rejects_'+str(sign),lambda:matrix.parse_export(p,row,units))
            badlines=lines.copy();badlines[6]='RESTORED F';p.write_text('\n'.join(badlines)+'\n');rejects('other_parser_gate_preserved_'+str(sign),lambda:proxy.parse_export(p,row,units))
    # Mock only the immutable context/native/parser interfaces. Real run() supplies order/caps.
    real_context,real_reparse,real_field,real_parse,real_copy=R.context,R.reparse,R.C.field_checks,R.C.parse,R.shutil.copyfile
    orig_run=subprocess.run
    def forbidden(*a,**kw):raise AssertionError('No subprocess/native permitted in controls')
    subprocess.run=forbidden
    try:
        for fail in ('none','field_parse','field_numeric','bypass_native'):
            d=out/fail;d.mkdir();jobs=[];rows={};metadata={};calls=[];stages=[]
            for i in range(42):
                g='saved'+str(i);job={'id':g,'stage':'active','cell':g};jobs.append(job);rows[g]={'status':'complete','literal':i}
            for g in P.DIMERS:
                # Existing active rows referenced by the unchanged field stage.
                oldkey='saved'+str(0 if g.endswith('La') else 1);job=jobs[0 if g.endswith('La') else 1];job['id']='active__'+g+'__analyze';rows[job['id']]=rows.pop(oldkey)
                for sign in (-1,1):
                    job={'id':'fields__'+g+'__analyze__'+str(sign),'stage':'fields','cell':g,'sign':sign};jobs.append(job);rows[job['id']]={'status':'failed','error':R.ERROR,'receipts':[{'returncode':0}]}
                metadata[g]={}
            with tempfile.TemporaryDirectory(prefix='radial_recovery_mock_') as td:
                f=Path(td);pins={}
                for k,name in P.FILES:(f/name).write_text(name);pins[k]=P.pin(f/name)
                for g in P.BYPASS:
                    job={'id':'bypass__'+g+'__analyze','stage':'bypass','cell':g,'input':pins};jobs.append(job);rows[job['id']]={'status':'not_attempted'};metadata[g]={}
                pc={'manifest':{},'source_manifest':{},'build':{},'executables':{'analyze':{'path':'FORBIDDEN_NATIVE'}}};old={'rows':rows,'error':R.ERROR,'cost':{},'operator_checks':{},'directional_checks':{}}
                acct=f/'acct.json';P.save(acct,{})
                def invoke(argv,folder,output,receipts):
                    calls.append(argv);receipts.append({'returncode':0});output.write_text('manufactured output')
                    if fail=='bypass_native':raise ValueError('manufactured native failure')
                class Parent:
                    @staticmethod
                    def invocation(h):return invoke
                R.context=lambda c,runtime:( {'parent_result':{},'parent_sacct':P.pin(acct)},Parent,pc,old,{'programs':jobs,'rows':metadata},{},None)
                def parse_field(rec,*args):
                    stages.append('field_parse')
                    if fail=='field_parse':raise ValueError('manufactured field parse failure')
                    rec.pop('error');rec['status']='complete';return rec
                def field(*args):
                    stages.append('field_numeric')
                    if fail=='field_numeric':raise ValueError('manufactured field numerical failure')
                    return {'passed':True}
                R.reparse=parse_field;R.C.field_checks=field;R.C.parse=lambda *args:None
                prior=os.environ.get('SLURM_JOB_ID');os.environ['SLURM_JOB_ID']='CONTROL_ONLY'
                try:r=R.run(cfg,d)
                finally:
                    if prior is None:del os.environ['SLURM_JOB_ID']
                    else:os.environ['SLURM_JOB_ID']=prior
                check(fail+'_old42_unchanged',r['reused42_rows_exact'])
                if fail=='none':check('success_exact3_total49',len(calls)==3 and r['native_programs_completed']==49 and r['numerical_checks_passed'])
                elif fail.startswith('field'):check(fail+'_zero_native_untouched_bypass',not calls and all(r['rows']['bypass__'+g+'__analyze']=={'status':'not_attempted'} for g in P.BYPASS) and r['global_stop'])
                else:check('first_bypass_failure_stops',len(calls)==1 and r['new_native_program_invocation_attempts']==1 and r['new_native_programs_completed']==0 and r['global_stop'])
    finally:R.context,R.reparse,R.C.field_checks,R.C.parse,R.shutil.copyfile=real_context,real_reparse,real_field,real_parse,real_copy;subprocess.run=orig_run
    result={'status':'passed_radial_field_recovery_controls','implementation':P.pin(__file__),'runner':P.pin(R.__file__),'field_adapter':P.pin(F.__file__),'checks':checks,'count':len(checks),'native_calls':0,'scientific_matrix_qualifications':0,'optimizer_calls':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'numerical_threads':1}}
    P.save(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'checks':len(checks)}))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=False)
    try:main(json.loads(Path(a.config).read_text()),out)
    except Exception as e:P.save(out/'FAILURE.json',{'error':repr(e),'native_calls':0});raise
