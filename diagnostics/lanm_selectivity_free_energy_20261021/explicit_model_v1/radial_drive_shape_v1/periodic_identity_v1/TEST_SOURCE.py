"""Focused source/manufactured identity controls. No physical box readback."""
import ast, copy, hashlib, io, json, math, tempfile, time, traceback
from pathlib import Path
import numpy as np
import IDENTITY as I
import STATE as S
import PREPARE_SOURCE as P

HERE=Path(__file__).resolve().parent
CHECKS=[]
def check(name,condition):
    if not condition:raise AssertionError(name)
    CHECKS.append({'name':name,'passed':True})
def reject(name,fn):
    try:fn()
    except (ValueError,KeyError,TypeError,IndexError):CHECKS.append({'name':name,'passed':True});return
    raise AssertionError('Expected rejection: '+name)
def rows(text):
    r=S.Rows.__new__(S.Rows);r.stream=io.StringIO(text);return r

def header(model):
    c=model['cell'];settings=model['settings'];t=np.diag(64/np.array(c))
    out=[I.SCHEMA,'UNITS xyz=angstrom alpha=angstrom^3 beta=angstrom^-3 h=electron/angstrom^2',
         'RHS LABELS raw_Ed raw_Ep; add_shared_h_to_both; no_solved_dipoles',
         'GRID 64 64 64 5 .4 '+repr(c[0]*c[1]*c[2]),'MODE FFTW TINFOIL MUTUAL ALL_ACTIVE SIMPLE_LOOP',
         'CELL '+' '.join(map(repr,c)),'CUTOFFS 12 9 8.1','EWALD_SPLIT .4 .4','INDUCTION 1e-8 100']
    for a in range(3):out.append('TRANS '+str(a+1)+' '+' '.join(map(repr,t[a].tolist())))
    for m in ('900','901'):out.append('COEFFICIENT '+m+' '+' '.join(map(repr,model['coefficients'][m])))
    for m in ('900','901'):out.append('RADIAL_COEFFICIENT '+m+' '+repr(model['radial'][m]))
    return '\n'.join(out)+'\n'

def tiny_body():
    sites=[dict(atom=1,type=900,xyz=[0.,0.,0.],alpha=1.5,igrid=[0,0,0],theta=np.zeros((3,5,2))),
           dict(atom=2,type=402,xyz=[2.5,0.,0.],alpha=1.,igrid=[0,0,0],theta=np.zeros((3,5,2)))]
    expected=[{k:s[k] for k in ('type','xyz','alpha')} for s in sites];out=[]
    for i,s in enumerate(sites,1):
        out+=['SITE '+f'{i} {s["type"]} {i} '+' '.join(map(str,s['xyz']+[s['alpha'],1.,.4])),
              f'RESPONSE {i} 0 0 0 0',f'DIRECT {i} 1 2 3 -4 -5 -6']
        if s['type']==402:out.append(f'PAIR_GAMMA {i} .4')
        out.append(f'IGRID {i} 0 0 0')
        for slot in range(1,6):
            for axis in range(1,4):out.append(f'THETA {i} {axis} {slot} 0 0')
    out+=['QFAC 1 1 1 .75','END '+''.join(f'{n:12d}' for n in (2,30,1,2))]
    return '\n'.join(out)+'\n',{'sites':sites,'qfac':np.full((1,1,1),.75)},expected

def main(output):
    wall=time.monotonic();cpu=time.process_time();out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    result={'status':'failed_periodic_identity_v2_controls','checks':CHECKS,
            'scope':'Source and tiny manufactured records only; no actual candidate or physical snapshot.',
            'counts':{'native_calls':0,'solver_calls':0,'FFT_replays':0,'factorizations':0,'compilations':0}}
    try:
        declaration=json.loads((HERE/'ROOT_DECLARATION.json').read_text());parents=declaration['parents']
        for p in parents.values():P.checked(p)
        check('All 25 parent pin authorities verified',len(parents)==25)
        src,edits,diff=P.sources(parents)
        for key,original,name in [('writer','state_writer_actual_parent','STATE_WRITER_V2.f90'),('gate','cutoff12_gate_actual_parent','GATE_V2.f90')]:
            restored=src[name]
            for old,new in reversed(edits[key]):restored=restored.replace(new,old,1)
            check(key+' inverse restores exact compiled parent',restored==P.checked(parents[original]))
        writer=src['STATE_WRITER_V2.f90'];gate=src['GATE_V2.f90']
        check('One shared body has both ordered radial slots','do slot=1,2' in writer and "'RADIAL_COEFFICIENT ',slot+899,lrpshape(slot)" in writer)
        check('Both wrappers call identical body with distinct raw RHS',src['STATE_ENTRY_POINTS.f90'].count('call periodic_write_state_v2(')==2 and writer.count("'DIRECT ',ia,direct1(:,ia),direct2(:,ia)")==1)
        check('Writer has no solve/field calls','call ' not in writer.lower())
        check('Gate retains actual intrinsic recovery and cutoff12', 'intrinsic :: maxval' in gate and 'ewaldcut/=12d0' in gate)
        check('Gate requires V2 ticket and V2 schema', 'STAGE2_BOUND_INPUT_V2' in gate and "trim(line)/='PERIODIC_PARITY_STATE_V2'" in gate)
        check('Exact comparison precedes mutual permission', "role=='INITIALIZE') return" in gate and "role=='ACTIVE'.and..not.certified" in gate)
        check('INITIALIZE cannot solve', "role=='INITIALIZE') error stop 'Stage2 unauthorized solve entry'" in gate)
        check('Radial gate compares coefficients, not h/H',"any(lrpshape/=0d0).and..not.lractive" in gate)
        for name in ('STATE.py','IDENTITY.py','PREPARE_SOURCE.py','TEST_SOURCE.py'):
            ast.parse((HERE/name).read_text());check('Python AST '+name,True)
        old=P.checked(parents['state_parser_cell_aware'])
        oldfun=old[old.index('def compact_response'):old.index('def read_current')].replace('def compact_response(', 'def _legacy_compact_response(')
        check('Exact original ordered h/H source preserved',oldfun in (HERE/'STATE.py').read_text())
        model={'cell':[48.,50.,52.],'settings':copy.deepcopy(I.SETTINGS),'coefficients':{'900':[.2,.3,.1,-.1],'901':[0.,0.,0.,0.]},'radial':{'900':.25,'901':-.75}}
        text=header(model);h=S.read_header(rows(text),model)
        check('V2 header keeps BOTH explicit nonzero e',h['radial']==model['radial'])
        variants={
            'V1 cannot fall back':text.replace(I.SCHEMA,'PERIODIC_PARITY_STATE_V1'),
            'Missing e900':text.replace('RADIAL_COEFFICIENT 900 0.25\n',''),
            'Duplicate e900':text.replace('RADIAL_COEFFICIENT 901 -0.75','RADIAL_COEFFICIENT 900 -0.75'),
            'Nonfinite e':text.replace('900 0.25','900 nan'),
            'Out of domain e':text.replace('900 0.25','900 1.01'),
            'Extra radial token':text.replace('900 0.25','900 0.25 0'),
            'Malformed integer radial type':text.replace('RADIAL_COEFFICIENT 900','RADIAL_COEFFICIENT 900.0'),
            'Cutoff9 rejected':text.replace('CUTOFFS 12','CUTOFFS 9'),
            'Missing CELL':text.replace('CELL 48.0 50.0 52.0\n',''),
            'Fractional iterations rejected':text.replace('1e-8 100','1e-8 100.5'),
            'Grid changed':text.replace('GRID 64 64 64','GRID 32 64 64'),
            'Cell volume mismatch':text.replace(repr(48.*50.*52.),'1.0'),
            'Cell transform mismatch':text.replace('TRANS 1 1.3333333333333333','TRANS 1 1.5')}
        for name,value in variants.items():reject(name,lambda value=value:S.read_header(rows(value),model))
        body,snapshot,expected=tiny_body();parsed=S._read_body(rows(body),snapshot,expected,(1,1,1))
        check('Tiny full body preserves unequal raw Ed/Ep',np.array_equal(parsed['direct'][0],[1,2,3,-4,-5,-6]))
        for name,bad in [('Trailing output',body+'GARBAGE\n'),('Bad END radial denominator',body[:-13]+f'{1:12d}\n'),('Missing QFAC',body.replace('QFAC 1 1 1 .75\n','')),('Reordered spline',body.replace('THETA 1 1 1','THETA 1 2 1',1)),('Raw nonfinite',body.replace('DIRECT 1 1','DIRECT 1 nan'))]:
            reject(name,lambda bad=bad:S._read_body(rows(bad),snapshot,expected,(1,1,1)))
        sites=expected;gamma=[0.,.4];pdamp=[1.,1.]
        zero=S.compact_response(sites,[.2,.3,.1,-.1],gamma,pdamp,[48.,50.,52.],0.)
        oldzero=S._legacy_compact_response(sites,[.2,.3,.1,-.1],gamma,pdamp,[48.,50.,52.])
        check('Zero e exactly reuses old h and H',all(np.array_equal(a,b) for a,b in zip(zero,oldzero)))
        r0a=S.compact_response(sites,[0.,0.,0.,0.],gamma,pdamp,[48.,50.,52.],.2)
        r0b=S.compact_response(sites,[0.,0.,0.,0.],gamma,pdamp,[48.,50.,52.],.7)
        check('r0 distinct e has identical h/H',all(np.array_equal(a,b) for a,b in zip(r0a,r0b)))
        check('r0 e-only is ACTIVE',I.validate_role('ACTIVE',{'900':[0.]*4,'901':[0.]*4},{'900':.7,'901':0.},True,True))
        reject('r0 e-only cannot bypass',lambda:I.validate_role('BYPASS_CLONE_ZERO',{'900':[0.]*4,'901':[0.]*4},{'900':.7,'901':0.},False,False))
        shifted=copy.deepcopy(sites);shifted[1]['xyz']=[45.5,0.,0.]
        imaged=S.compact_response(shifted,[0.]*4,gamma,pdamp,[48.,50.,52.],.7)
        othercell=S.compact_response(shifted,[0.]*4,gamma,pdamp,[50.,50.,52.],.7)
        check('Declared cell changes nearest-image radial response',not np.array_equal(imaged[0],othercell[0]))
        contact=copy.deepcopy(sites);contact[1]['xyz']=[0.,0.,0.]
        check('Contact radial h finite zero',np.array_equal(S.compact_response(contact,[0.]*4,gamma,pdamp,[48.,50.,52.],1.)[0],np.zeros((2,3))))
        outside=copy.deepcopy(sites);outside[1]['xyz']=[21.,0.,0.]
        check('Compact endpoint is zero',np.array_equal(S.compact_response(outside,[0.]*4,gamma,pdamp,[48.,50.,52.],1.)[0],np.zeros((2,3))))
        stages={k:copy.deepcopy(I.SETTINGS) for k in ('key','gate','serialized','parser','assembler','certificate')}
        check('Six cutoff12 settings joins pass',I.validate_settings_joins(**stages))
        for stage in stages:
            bad=copy.deepcopy(stages);bad[stage]['ewaldcut_A']=9.
            reject('Cutoff corrupt '+stage,lambda bad=bad:I.validate_settings_joins(**bad))
        # Tiny immutable metadata files bind exactly the chain, including central.
        with tempfile.TemporaryDirectory(prefix='identity-v2-control-') as td:
            td=Path(td)
            def file(name,obj):
                f=td/name;f.write_text(json.dumps(obj,sort_keys=True)+'\n');return P.pin(f)
            sample=file('sample',{'manufactured':True})
            identity={'schema':I.SCHEMA,'state_id':'central',**copy.deepcopy(model),'inputs':{k:sample for k in ('xyz','key','parameters')},
                      'source_manifest':sample,'build':sample,'state_output':sample,'pme_export':sample}
            check('Identity exact self join',I.validate_identity(identity,identity))
            for label,mut in [('Both e identity at r0',lambda v:v['radial'].__setitem__('900',.7)),('Unused metal e identity',lambda v:v['radial'].__setitem__('901',0.)),('Source identity',lambda v:v['source_manifest'].__setitem__('sha256','0'*64))]:
                bad=copy.deepcopy(identity);mut(bad);reject(label,lambda bad=bad:I.validate_identity(bad,identity))
            nom={'U_B':3.,'U_R':4.,'delta_odd':1e-12,'U_nominal':7.,'tau':1e-6*7+1e-12,'e_app':None,'application_error_bound_established':False}
            records={};pins={}
            def put(k,status,**extra):
                records[k]={'status':status,'identity':copy.deepcopy(identity),**extra};pins[k]=file(k,records[k])
            put('initializer','stage2_initialization_complete',state_output=sample,pme_export=sample)
            put('assembly','stage2_assembly_complete',initializer=pins['initializer'],cutoff_A=12.,ewald_inverse_A=.4)
            put('factor','passed_nominal_static_sufficient_numerical_test',assembly=pins['assembly'],nominal=nom)
            put('certificate','stage2_state_nominal_certificate_passed',initializer=pins['initializer'],assembly=pins['assembly'],factor=pins['factor'],expected_state=sample,pme_export=sample,nominal=nom,sufficient_numerical_test_passed=True)
            put('root_admission','admitted_periodic_identity_v2_certificate',certificate=pins['certificate'],initializer=pins['initializer'],fresh_same_state_chain=True)
            check('Fresh central full chain passes',I.validate_bound_chain(identity,pins)['state_id']=='central')
            bad=copy.deepcopy(records);bad['certificate']['initializer']=None
            reject('Central cannot implicitly reuse old certificate',lambda:I._validate_chain_records(identity,bad,pins))
            bad=copy.deepcopy(records);bad['assembly']['cutoff_A']=9.
            reject('Actual assembly cutoff cannot differ',lambda:I._validate_chain_records(identity,bad,pins))
            bad=copy.deepcopy(records);bad['root_admission']['fresh_same_state_chain']=False
            reject('Root must admit explicit fresh chain',lambda:I._validate_chain_records(identity,bad,pins))
            bad=copy.deepcopy(records);bad['factor']['nominal']['tau']*=2
            reject('Nominal shift unchanged',lambda:I._validate_chain_records(identity,bad,pins))
            badpins=copy.deepcopy(pins);badpins['initializer']['sha256']='0'*64
            reject('Raw chain file hash required',lambda:I.validate_bound_chain(identity,badpins))
            # Legacy comparison is diagnostic only, not a synthetic runnable state.
            m0=copy.deepcopy(model);m0['radial']={'900':0.,'901':0.};v2=td/'v2';v2.write_text(header(m0)+body)
            oldlines=[l for l in v2.read_text().splitlines(True) if not l.startswith(('CELL ','CUTOFFS ','EWALD_SPLIT ','INDUCTION ','RADIAL_COEFFICIENT '))]
            oldlines[0]='PERIODIC_PARITY_STATE_V1\n';oldlines[-1]=oldlines[-1][:-13]+'\n';v1=td/'v1';v1.write_text(''.join(oldlines))
            authority={'legacy_schema':'PERIODIC_PARITY_STATE_V1','zero_radial_proven':True,'cell':m0['cell'],'settings':copy.deepcopy(I.SETTINGS),
                       'proof_pins':{'state':P.pin(v1),'key':sample,'composition':sample,'source':sample,'admission':sample}}
            check('Same settings legacy zero byte projection',S.compare_legacy_zero(v1,v2,authority,m0)['expected_state_created'] is False)
            bad=copy.deepcopy(authority);bad['settings']['ewaldcut_A']=9.
            reject('No 9 to12 legacy reuse',lambda:S.compare_legacy_zero(v1,v2,bad,m0))
            bad=copy.deepcopy(authority);del bad['proof_pins']['key']
            reject('No missing legacy authority',lambda:S.compare_legacy_zero(v1,v2,bad,m0))
            bad=copy.deepcopy(authority);bad['zero_radial_proven']=False
            reject('No inferred missing-e zero',lambda:S.compare_legacy_zero(v1,v2,bad,m0))
            v2.write_text(header(model)+body)
            reject('Legacy projection nonzero e forbidden',lambda:S.compare_legacy_zero(v1,v2,authority,model))
        result['status']='passed_periodic_identity_v2_source_controls'
    except Exception as exc:
        result['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    result['check_count']=len(CHECKS)
    result['implementation']=P.pin(__file__);result['source_pins']={n:P.pin(HERE/n) for n in ('STATE.py','IDENTITY.py','PREPARE_SOURCE.py')}
    result['cost']={'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'threads':1}
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'checks':len(CHECKS),'result':P.pin(out/'RESULT.json'),'error':result.get('error'),'cost':result['cost']}))
    if result['status'].startswith('failed'):raise SystemExit(1)
if __name__=='__main__':
    import sys
    main(sys.argv[1])
