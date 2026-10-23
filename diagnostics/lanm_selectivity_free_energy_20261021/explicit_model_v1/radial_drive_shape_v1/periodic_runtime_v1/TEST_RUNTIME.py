"""Finite manufactured/source controls. Never calls native, FFT, factor or solve."""
import ast,copy,hashlib,json,tempfile,time,traceback
from pathlib import Path
import COMMON as C
import POLICY as P
import METADATA as M
import GRADIENT_PARSE as G

HERE=Path(__file__).resolve().parent;CHECKS=[]
def check(name,ok):
    if not ok:raise AssertionError(name)
    CHECKS.append({'name':name,'passed':True})
def reject(name,fn):
    try:fn()
    except (ValueError,TypeError,KeyError,IndexError):CHECKS.append({'name':name,'passed':True});return
    raise AssertionError('Expected rejection '+name)
def function(text,name):
    return next(n for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name==name)
def ast_equal(a,b):return ast.dump(a,include_attributes=False)==ast.dump(b,include_attributes=False)
def pin(p):return C.pin(p)

def gradient(role,radial=(.2,-.3),v2=True):
    def f(x):return format(x,'.17E')
    out=['STAGE2_ONE_ANALYTIC_GRADIENT_V1','UNITS E=kcal/mol G=positive_dE/dr_kcal/mol/angstrom VIRIAL=kcal/mol',
         'TOTAL_ENERGY '+f(-3.2),'POLARIZATION_ENERGY '+f(-3.2),
         'CONTRACTION_CONSTANTS electric dielec '+f(1.)+' '+f(1.),
         'TOTAL_GRADIENT 1 '+' '.join([f(0.)]*3)]
    out += ['TOTAL_VIRIAL '+str(i)+' '+' '.join([f(0.)]*3) for i in range(1,4)]
    if v2:out+=['STATE_SCHEMA PERIODIC_PARITY_STATE_V2','RADIAL_COEFFICIENT 900 '+f(radial[0]),'RADIAL_COEFFICIENT 901 '+f(radial[1])]
    out+=['ROLE '+role,'STATE_ID central','INPUT_BINDING '+'a'*64]
    if role=='BYPASS_OLD_ABSENT':out+=['INTERNAL_COUNTERS unavailable_unchanged_old_archive']
    else:
        out+=['COUNTS induced_solve direct mutual 1 1 4','COUNTS extra_diagnostic_operator_applications=0 internal_FD=0',
              'PCG_PRE_PEEK 3 '+' '.join(map(f,[1e-9,1e-18,2e-18,2e-9]))]
        if role=='ACTIVE':
            out+=['SHARED_H_READBACK_MAX_ERROR '+f(0.),
                  'RAW_DIRECT 1 '+' '.join(map(f,[.4,.5,.6,.7,.8,.9])),
                  'CORRECTED_RHS 1 '+' '.join(map(f,[.5,.7,.9,.8,1.,1.2])),
                  'FINAL_POST_PEEK 1 '+' '.join(map(f,[1.,1.,2.,3.,-1.,-2.,-3.])),
                  'CORRECTION 1 '+' '.join(map(f,[0.,.1,.2,.3]))]
    out+=['END gradient_calls=1 energy_only_calls=0 internal_FD=0']
    return '\n'.join(out)+'\n'

def main(output,native_manifest):
    wall=time.monotonic();cpu=time.process_time();out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    result={'status':'failed_periodic_v2_runtime_source_controls','checks':CHECKS,'error':None,
            'counts':{'native_calls':0,'builds':0,'FFT_replays':0,'assemblies':0,'factorizations':0,'solver_calls':0}}
    try:
        decl=json.loads((HERE/'ROOT_DECLARATION.json').read_text());parents=decl['parents']
        for p in parents.values():C.verify(p)
        check('All declared parent pins verified',len(parents)==39)
        native=json.loads(Path(native_manifest).read_text());inverse=C.read(native['inverse'])
        for p in native['files'].values():C.verify(p)
        def native_text(name):return C.verify(native['files'][name]).read_text()
        restored=native_text('induce_qualification.f')
        for a,b in reversed(inverse['qualification']):
            check('Unique inverse qualification anchor '+str(len(CHECKS)),restored.count(b)==1);restored=restored.replace(b,a,1)
        check('Exact whole radial induce inverse',restored==C.verify(inverse['parent_induce']).read_text())
        restored=native_text('PARITY_V2.f90')
        for a,b in reversed(inverse['parity']):restored=restored.replace(b,a,1)
        check('Exact Stage1 parent inverse',restored==C.verify(parents['stage1_harness']).read_text())
        restored=native_text('PARITY_SETUP.f')
        for a,b in reversed(inverse['setup']):restored=restored.replace(b,a,1)
        parent=C.verify(inverse['parent_induce']).read_text();a=parent.index('      subroutine periodicresponsesetup (handled)\n');b=parent.index('      subroutine periodicresponsebackendguard\n',a)
        check('Radial setup inverse and e parser retained',restored==parent[a:b] and 'call periodicradialsetup' in native_text('PARITY_SETUP.f'))
        gate=native_text('GATE_V2.f90');parity=native_text('PARITY_V2.f90')
        check('Shared V2 writer/gate exact frozen identity sources',all(native['files'][k]['sha256']==p['sha256'] for k,p in C.read(parents['identity_source_manifest'])['outputs'].items()))
        check('Parity full-state comparison before any operator',parity.index('call compare_expected_v2')<parity.index('call ufield0(field1,field2)'))
        check('Both e/coeff persistent identity', 'any(lrpshape/=old_e)' in parity and 'any(lrcoef/=old_coef)' in parity)
        check('Both normal fatal bodies retained',native['normal_fatal_guards_unchanged'] is True)
        check('Initializer still has zero solve/energy/gradient calls',not any('call '+x in native_text('INITIALIZE.f90').lower() for x in ('induce','energy','gradient','ufield')))
        source={f.name:f.read_text() for f in HERE.glob('*.py')}
        for n,s in source.items():ast.parse(s)
        check('All runtime Python ASTs parse',len(source)>=9)
        old_common=C.verify(parents['transaction_wrapper']).read_text()
        check('Invocation/allocation/persistence helpers exact',source['COMMON.py'].startswith(old_common[:old_common.index('def source_context(')]))
        oldcert=C.verify(parents['certificate_producer']).read_text()
        check('Numerical factor routine AST unchanged',ast_equal(function(oldcert,'factor'),function(source['CERTIFY.py'],'factor')))
        oldpar=C.verify(parents['stage1_parser_replay']).read_text()
        for name in ('scatter_matrix','complex_reciprocal','real_action','compare'):
            check('Exact parity numerical helper '+name,ast_equal(function(oldpar,name),function(source['PARITY_REPLAY.py'],name)))
        oldnative=C.verify(parents['native_wrapper']).read_text()
        check('Literal site reader AST unchanged',ast_equal(function(oldnative,'sites'),function(source['METADATA.py'],'sites')))
        oldparse=C.verify(parents['gradient_parser']).read_text();v1=function(source['GRADIENT_PARSE.py'],'_parse_v1');v1.name='parse_gradient'
        check('Explicit old bypass parser unchanged',ast_equal(function(oldparse,'parse_gradient'),v1))
        check('Actual assembler arguments use validated12 settings',"cutoff=model['settings']['ewaldcut_A'],ewald=model['settings']['apewald_inverse_A']" in source['CERTIFY.py'])
        # No call to any runtime invocation function is made by these controls.
        check('Central initialize has no parity/scout dependency',P.phase_needs('initialize','central','ACTIVE')=={'central_parity':False,'central_scout':False,'own_certificate':False})
        check('Central factor bootstrap has no parity dependency',not P.phase_needs('certificate','central','ACTIVE')['central_parity'])
        check('Central parity requires own certificate only',P.phase_needs('parity','central','ACTIVE')=={'central_parity':False,'central_scout':False,'own_certificate':True})
        check('Central gradient needs new parity and own certificate',P.phase_needs('gradient','central','ACTIVE')=={'central_parity':True,'central_scout':False,'own_certificate':True})
        for op in ('initialize','certificate','gradient'):
            need=P.phase_needs(op,'strain_1_plus','ACTIVE');check('Displaced phase needs campaign parity/scout '+op,need['central_parity'] and need['central_scout'])
        reject('No32 implicit displaced parity programs',lambda:P.phase_needs('parity','strain_1_plus','ACTIVE'))
        for role in P.BYPASSES:reject('No bypass initializer '+role,lambda role=role:P.phase_needs('initialize',role,role))
        check('Old baseline gets V1 ticket only',P.ticket('gradient','BYPASS_OLD_ABSENT','BYPASS_OLD_ABSENT','a'*64).splitlines()[0]=='STAGE2_BOUND_INPUT_V1')
        check('Radial initializer gets V2 schema ticket',P.ticket('initialize','ACTIVE','central','a'*64).splitlines()[:3]==['STAGE2_BOUND_INPUT_V2','PERIODIC_PARITY_STATE_V2','INITIALIZE'])
        reject('Ticket newline injection',lambda:P.ticket('gradient','ACTIVE','central\nEND','a'*64))
        text='\n'.join(['ewald','ewald-cutoff 12','vdw-cutoff 9','vdw-taper 8.1','ewald-alpha .4','pewald-alpha .4','pme-grid 64 64 64','pme-order 5','ppme-order 5','fft-package FFTW','polarization MUTUAL','polar-iter 100','polar-eps 1e-8','PME-NORM-EXPORT 1',''])
        check('Explicit native key settings',M.key_settings(text)==C.I.SETTINGS)
        for name,bad in [('cutoff9',text.replace('cutoff 12','cutoff 9')),('boundary keyword',text+'ewald-boundary TINFOIL\n'),('missing eps',text.replace('polar-eps 1e-8\n','')),('duplicate',text+'ewald-cutoff 12\n')]:reject('Bad settings '+name,lambda bad=bad:M.key_settings(bad))
        response='polar-local-response 900 402 0 0 0 0\npolar-local-response 901 402 0 0 0 0\npolar-radial-drive 900 402 .2\npolar-radial-drive 901 402 -.3\n'
        coefs,radial=M.response_records(['',response]);check('Explicit two e-only slots parse',coefs=={'900':[0.]*4,'901':[0.]*4} and radial=={'900':.2,'901':-.3})
        for name,bad in [('missing slot',response.replace('polar-radial-drive 901 402 -.3\n','')),('duplicate',response+'polar-radial-drive 900 402 .2\n'),('nonfinite',response.replace('402 .2','402 nan')),('outdomain',response.replace('402 .2','402 1.1')),('trailing',response.replace('402 .2','402 .2 4'))]:reject('Bad radial '+name,lambda bad=bad:M.response_records([bad]))
        reject('Cross input override forbidden',lambda:M.response_records([response,response]))
        with tempfile.TemporaryDirectory(prefix='periodic-runtime-control-') as td:
            td=Path(td);f=td/'artifact';f.write_text('manufactured');pp=pin(f)
            m={'state_schema':C.I.SCHEMA,'campaign':pp,'native_sources':pp,'model':{'candidate':pp,'composition':pp,'settings':dict(C.I.SETTINGS),'coefficients':coefs,'radial':radial},
               'rows':{k:{'role':'ACTIVE','cell':c,'inputs':{s:pp for s in ('xyz','key','parameters')}} for k,c in [('central',[48.,48.,48.]),('strain_1_plus',[48.,48.01,48.])]}}
            central=P.full_identity(m,pp,'central',pp,pp);displaced=P.full_identity(m,pp,'strain_1_plus',pp,pp)
            admission={'status':'admitted_periodic_v2_central_operator_parity','campaign_identity':P.campaign_identity(m,pp),'central_identity':central,'sampled_operator_parity_passed':True,'result':pp,'certificate_admission':pp}
            check('Campaign parity admits scope independently of displaced coordinates',P.parity_scope(admission,m,pp) and central!=displaced)
            C.I.validate_identity(displaced,P.full_identity(m,pp,'strain_1_plus',pp,pp));check('Displaced own certificate remains exact',True)
            reject('Displaced cannot claim central full state',lambda:C.I.validate_identity(displaced,central))
            bad=copy.deepcopy(admission);bad['campaign_identity']['radial']['901']=0.
            reject('Unused e changes campaign identity',lambda:P.parity_scope(bad,m,pp))
            scout={'status':'admitted_periodic_v2_central_and_bypass_scout','campaign_identity':P.campaign_identity(m,pp),'central_parity_admission':pp,'all_three_bypass_tokens_equal':True,'central_dual_energy_passed':True}
            check('Actual central scout scope joins',P.scout_scope(scout,m,pp,pp))
            bad=copy.deepcopy(scout);bad['central_dual_energy_passed']=False
            reject('Failed central dual gate blocks changed states',lambda:P.scout_scope(bad,m,pp,pp))
            g=td/'gradient';g.write_text(gradient('ACTIVE'));got=G.parse_gradient(g,'ACTIVE',m,n=1)
            check('V2 full manufactured dual readback passes',got['shared_h']['passed'] and got['polarization_contraction']['passed'] and got['radial']==radial)
            check('Unequal solved channels preserved',got['polarization_contraction']['maximum_u_v_difference']==6.)
            for name,body in [('wrong e',gradient('ACTIVE',(.4,-.3))),('missing e',gradient('ACTIVE').replace('RADIAL_COEFFICIENT 901 -2.99999999999999989E-01\n','')),('V1 fallback',gradient('ACTIVE',v2=False))]:
                g.write_text(body);reject('Strict gradient '+name,lambda:G.parse_gradient(g,'ACTIVE',m,n=1))
            g.write_text(gradient('BYPASS_OLD_ABSENT',v2=False));old=G.parse_gradient(g,'BYPASS_OLD_ABSENT',m,n=1)
            g.write_text(gradient('BYPASS_CLONE_ABSENT',(0.,0.)));clone=G.parse_gradient(g,'BYPASS_CLONE_ABSENT',m,n=1)
            check('Old/V2 bypass E/G/V token parity survives metadata',old['E_G_virial_token_sha256']==clone['E_G_virial_token_sha256'])
            reject('No automatic old parser fallback',lambda:G.parse_gradient(g,'BYPASS_OLD_ABSENT',m,n=1))
        result['status']='passed_periodic_v2_runtime_source_controls'
    except Exception as exc:result['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    result.update(check_count=len(CHECKS),implementation=pin(__file__),native_manifest=pin(native_manifest),source_pins={f.name:pin(f) for f in HERE.glob('*.py')},cost={'wall_seconds':time.monotonic()-wall,'CPU_seconds':time.process_time()-cpu,'threads':1})
    C.save(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'checks':len(CHECKS),'error':result['error'],'result':pin(out/'RESULT.json'),'cost':result['cost']}))
    if result['error']:raise SystemExit(1)
if __name__=='__main__':
    import sys
    main(sys.argv[1],sys.argv[2])
