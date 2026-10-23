"""Tiny saved-metadata and mocked-binder checks; no actual scientific outcomes."""
import argparse
import copy
from pathlib import Path
import tempfile
from types import SimpleNamespace
from unittest.mock import patch
import BATCH as B


def run(out):
    tick=B.time.monotonic();cpu=B.time.process_time();checks=[]
    def check(name,ok):B.require(ok,name);checks.append(name)
    def rejected(name,fn):
        try:fn()
        except (ValueError,KeyError):checks.append(name)
        else:raise AssertionError(name)
    with tempfile.TemporaryDirectory() as td:
        p=Path(td)
        def doc(name,obj):q=p/name;q.parent.mkdir(parents=True,exist_ok=True);B.save(q,obj);return B.pin(q)
        token=p/'raw';token.write_text('tiny manufactured raw receipt\n');tp=B.pin(token)
        slots=['fixed_slot_'+str(i) for i in range(32)]
        m={'native_sources':tp,'rows':{s:{'role':'ACTIVE','inputs':{},'source_inputs':{}} for s in slots},'execution_sources':{'BIND.py':tp,'NATIVE.py':tp}}
        mp=doc('manifest.json',m)
        cs={'execution_manifest':mp,'certificate_sources':{'BIND.py':tp,'CERTIFY.py':tp}}
        sources={'manifest':mp,'execution_sources':m['execution_sources'],'certificate_source_inputs':doc('certificate_sources.json',cs)}
        sp=doc('sources.json',sources);build=doc('build.json',{'status':'stage2_build_complete','manifest':mp,'source_manifest':tp})
        scout=doc('scout.json',{'status':'admitted_stage2_central_and_bypass_scout','manifest':mp,'build':build,'all_three_bypass_tokens_equal':True,'central_dual_energy_passed':True})
        cfg={'status':'bound_stage2_metadata_batch','implementation':B.pin(B.__file__),'sources':sp,'build':build,'central_scout_review':scout,'phase':'initializers','upstream_results':dict.fromkeys(slots),'existing_attempts':dict.fromkeys(slots),'root_batch_admission':None}
        check('exact32 original slot inventory',B.ids(m)==slots)
        B.base(cfg);checks.append('matching build and admitted scout accepted')
        bad=copy.deepcopy(cfg);bad['upstream_results'].pop(slots[0]);rejected('missing map slot rejected',lambda:B.base(bad))
        bad=copy.deepcopy(cfg);bad['build']=doc('failedbuild.json',{'status':'stage2_build_failed'});rejected('failed build rejected',lambda:B.base(bad))
        bad=copy.deepcopy(cfg);bad['central_scout_review']=doc('badscout.json',{'status':'admitted_stage2_central_and_bypass_scout','manifest':mp,'build':build,'all_three_bypass_tokens_equal':False,'central_dual_energy_passed':True});rejected('failed scout parity rejected',lambda:B.base(bad))
        calls=[]
        def mock_bind(a):
            calls.append(a.state_id);a.out.mkdir(parents=True,exist_ok=False);B.save(a.out/'CONFIG.json',{'state_id':a.state_id});B.save(a.out/'ROOT_RELEASE_UNBOUND.json',{'status':None})
        fake=SimpleNamespace(run=mock_bind)
        cfg['existing_attempts'][slots[0]]=doc('failed_existing.json',{'state_id':slots[0],'manifest':mp,'status':'stage2_native_failed','error':'preserved','cost':{'CPU_seconds':1}})
        cp=doc('cfg_init.json',cfg)
        with patch.object(B,'load',return_value=fake):B.run(B.verify(cp),p/'initializers')
        result=B.read(B.pin(p/'initializers/RESULT.json'))
        check('existing failed attempt never rebound',slots[0] not in calls and len(calls)==31 and result['rows'][slots[0]]['producer_error']=='preserved')
        check('complete denominator retained',len(result['rows'])==32 and result['counts']['original_slots']==32)
        try:B.run(B.verify(cp),p/'initializers')
        except FileExistsError:checks.append('immutable output cannot be reused')
        else:raise AssertionError('output overwrite')
        cfg['phase']='certificates';cfg['existing_attempts']=dict.fromkeys(slots)
        good=doc('good_upstream.json',{'status':'stage2_initialization_complete'});bad=doc('bad_upstream.json',{'status':'stage2_native_failed','error':'native stopped'})
        cfg['upstream_results'][slots[0]]=good;cfg['upstream_results'][slots[1]]=bad
        def mock_admit(item,*args):
            if item==bad:raise ValueError('native stopped')
            return {}
        calls.clear();cp=doc('cfg_cert.json',cfg)
        with patch.object(B,'load',return_value=fake),patch.object(B,'admit_initializer',side_effect=mock_admit):B.run(B.verify(cp),p/'certificates')
        result=B.read(B.pin(p/'certificates/RESULT.json'))
        check('only terminal-admitted initializer reaches binder',calls==[slots[0]])
        check('failed and missing upstream rows retained',result['rows'][slots[1]]['producer_error']=='native stopped' and result['rows'][slots[2]]['status']=='missing_upstream_result')
        cfg['phase']='gradients';cfg['upstream_results'][slots[1]]=None;calls.clear();cp=doc('cfg_proposal.json',cfg)
        with patch.object(B,'load',return_value=fake),patch.object(B,'admit_certificate',return_value={}):B.run(B.verify(cp),p/'proposal')
        proposal_pin=B.pin(p/'proposal/PROPOSED_ADMISSION.json');proposal=B.read(proposal_pin)
        check('no gradient config before root acceptance',calls==[] and not (p/'proposal/bindings').exists())
        root={'status':'accepted_stage2_certificate_batch_admission','proposal':proposal_pin,**{k:proposal[k] for k in ('manifest','build','central_scout_review','accepted_certificates')}}
        rootpin=doc('root.json',root);B.root_acceptance(rootpin,proposal_pin,proposal);checks.append('exact one-map root acceptance accepted')
        badroot=copy.deepcopy(root);badroot['accepted_certificates']={};rejected('changed root accepted map rejected',lambda:B.root_acceptance(doc('badroot.json',badroot),proposal_pin,proposal))
        cfg['root_batch_admission']=rootpin;cp=doc('cfg_gradient.json',cfg)
        with patch.object(B,'load',return_value=fake),patch.object(B,'admit_certificate',return_value={}):B.run(B.verify(cp),p/'gradient')
        result=B.read(B.pin(p/'gradient/RESULT.json'));receipt=B.read(result['rows'][slots[0]]['state_admission'])
        check('root batch expanded into exact state receipt',calls==[slots[0]] and receipt['certificate']==good and receipt['root_batch_admission']==rootpin)
        # Real initializer metadata-admission checks on a tiny fabricated receipt.
        state=slots[0];where=p/'initializer_result';where.mkdir()
        nativecfg={'status':'bound_stage2_native_configuration','operation':'initialize','manifest':mp,'state_id':state,'build':build,'central_scout_review':scout,'implementation':tp}
        ncp=doc('nativecfg.json',nativecfg)
        rel={'status':'released_stage2_native','operation':'initialize','maximum_native_programs':1,'implementation':tp,'configuration':ncp,'state_id':state,'output':str(where)};rp=doc('native_release.json',rel)
        nr={'status':'stage2_initialization_complete','error':None,'implementation':tp,'state_id':state,'manifest':mp,'build':build,'inputs':{},'source_inputs':{},'source_manifest':tp,'configuration':ncp,'release':rp,'counts':{'native_programs':1,'direct_initializations':1,'mutual_applications':0,'induced_solves':0,'gradient_calls':0,'energy_only_calls':0},'PME_checks':{'site_checks':{'passed':True},'summary_comparison':{'passed':True},'operator_norms':{'q_negative_count':0,'q_min':0.}},'state_output':tp,'pme_export':tp,'native_attempts':1,'commands':{'native':{'status':'complete','returncode':0,'stdout':tp,'stderr':tp}},'initialization_receipt':tp}
        npin=doc('initializer_result/RESULT.json',nr)
        B.admit_initializer(npin,state,cfg,sources,m);checks.append('declared initializer counts and lineage admitted without replay')
        mutate=copy.deepcopy(nr);mutate['counts']['induced_solves']=1;rejected('hidden initializer solve rejected',lambda:B.admit_initializer(doc('initializer_result/RESULT.json',mutate),state,cfg,sources,m))
        mutate=copy.deepcopy(nr);mutate['commands']['native']['returncode']=1;rejected('nonzero native return rejected',lambda:B.admit_initializer(doc('initializer_result/RESULT.json',mutate),state,cfg,sources,m))
        mutate=copy.deepcopy(nr);mutate['state_id']=slots[1];rejected('wrong upstream state rejected',lambda:B.admit_initializer(doc('initializer_result/RESULT.json',mutate),state,cfg,sources,m))
    out.mkdir(parents=True,exist_ok=False)
    result={'status':'passed_stage2_metadata_batch_controls','implementation':B.pin(__file__),'source':B.pin(B.__file__),'checks':checks,'count':len(checks),'actual':{'wall_seconds':B.time.monotonic()-tick,'CPU_seconds':B.time.process_time()-cpu,'threads':1,'native_calls':0,'numerical_imports':0,'assemblies':0,'factors':0,'FFT_calls':0,'submissions':0}}
    B.save(out/'RESULT.json',result);print(B.json.dumps(result))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.out)
