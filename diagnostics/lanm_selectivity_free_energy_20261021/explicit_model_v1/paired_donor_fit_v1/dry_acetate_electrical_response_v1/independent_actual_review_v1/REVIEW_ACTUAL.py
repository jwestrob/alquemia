#!/usr/bin/env python3
"""One saved admission plus independent raw Decimal response/origin replay; no native calls."""
import hashlib,importlib.util,json,math,re,time
from pathlib import Path
from decimal import Decimal as D,getcontext
getcontext().prec=60
ROOT=Path.cwd();REL=Path('explicit_model_v1/paired_donor_fit_v1/dry_acetate_electrical_response_v1')
DB=ROOT/'diagnostics/lanm_selectivity_free_energy_20261021'/REL;WB=ROOT/'workspaces/lanm_selectivity_free_energy_20261021'/REL
OUT=WB/'independent_actual_review_v1';OUT.mkdir(exist_ok=False)
t0=time.monotonic();c0=time.process_time();checks=[];diffs=[]
def pin(p):
 p=Path(p).resolve();return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def read(p):return json.loads(Path(p).read_text())
def ck(ok,name):
 checks.append(dict(check=name,passed=bool(ok)));assert ok,name
def close(x,y,name):
 dx=abs(float(x)-float(y));diffs.append(dx);ck(dx<=2e-10*max(1,abs(float(x))),name)
def half(t):return D('0.5')*D(10)**D(t).as_tuple().exponent
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def verify(p):ck(pin(p['path'])==p,'pin:'+p['path']);return Path(p['path'])
RPATH=WB/'run_v1/execution/RESULT.json';r=read(RPATH);m=read(r['manifest']['path']);s=r['summary']
ck(pin(RPATH)['sha256']=='71d831858cbf9cc1c7faca24bb61a6a1536d025f96e3b56cc0e1cd91848facd4','frozen actual result')
sp=importlib.util.spec_from_file_location('frozen_field40',verify(r['implementation']));mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod)
ck(mod.admit_result(verify(r['manifest']),RPATH)==r,'single full saved producer admission')
K=D(str(m['conversion']['HARTREE']));N=r'[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?'
components=['Bond Stretching','Angle Bending','Stretch-Bend','Out-of-Plane Bend','Torsional Angle','Van der Waals','Atomic Multipoles','Polarization']
def raw(p):
 text=verify(p).read_text();tok=re.findall(r'^\s*Total Potential Energy\s*:\s*('+N+r')\s+Kcal/mole',text,re.M);ck(len(tok)==1,'one total:'+p['path']);energies={}
 for comp in components:
  a=re.findall(r'^\s*'+re.escape(comp)+r'\s+('+N+r')\s+(\d+)\s*$',text,re.M);ck(len(a)==1,'one component:'+comp+':'+p['path']);energies[comp]=a[0][0]
 energies['Total']=tok[0];ck(abs(D(tok[0])-sum((D(energies[c]) for c in components),D(0)))<=sum((half(t) for t in energies.values()),D(0)),'raw component closure:'+p['path']);return energies
raws={}
for jid,rec in r['rows'].items():
 toks=raw(rec['receipts'][0]['output']);ck(toks==rec['parsed']['tokens'],'raw tokens:'+jid);raws[jid]=toks['Total']
ck(len(raws)==40 and len(m['rows'])==10,'40 fields/10 zeros');states={};gates={};qms={}
for gid,row in m['rows'].items():
 zero=raw(row['zero_energy']['native']['analyze_output'])['Total'];ck(zero==row['zero_energy']['tokens']['Total'],'raw zero:'+gid)
 # Independently recover the common quantum origin from raw property vectors and effective nuclear charges.
 prop=verify(row['quantum']['property_file']).read_text();geo=prop.split('$Geometry\n',1)[1].split('$End',1)[0];xyz=re.findall(r'^\s*([A-Za-z]+)\s+('+N+r')\s+('+N+r')\s+('+N+r')\s*$',geo,re.M);ck(len(xyz)==8,'property sites:'+gid)
 ds=prop.split('$SCF_Dipole_Moment\n',1)[1].split('$End',1)[0]
 def vec(label):
  body=ds.split('&'+label,1)[1].split('&',1)[0] if False else re.search(r'&'+label+r'[^\n]*\n(.*?)(?=\n\s*&|$)',ds,re.S).group(1)
  arr=re.findall(r'^\s*([012])\s+('+N+r')\s*$',body,re.M);ck([x[0] for x in arr]==['0','1','2'],'dipole order:'+gid+'/'+label);return [D(x[1]) for x in arr]
 total,nuc,ele=vec('dipoleTotal'),vec('dipoleNucContrib'),vec('dipoleElecContrib');z=[11,6,6,8,8,1,1,1]
 origin=[(sum((D(z[i])*D(xyz[i][k+1]) for i in range(8)),D(0))-nuc[k])/42 for k in range(3)];lab=[total[k]+2*origin[k] for k in range(3)]
 for k in range(3):
  close(total[k],nuc[k]+ele[k],gid+'/quantum component');close(origin[k],row['quantum']['origin_Bohr'][k],gid+'/origin');close(lab[k],row['quantum']['total_lab_au'][k],gid+'/lab')
 qms[gid]=sum((lab[k]*D(str(m['direction'][k])) for k in range(3)),D(0));close(qms[gid],row['quantum_projected_au'],gid+'/QM projection')
 ck(sum(map(D,row['charge_origin']['site_monopoles_e']))==D(2) and row['charge_origin']['origin_A']==[0,0,0],'charged lab-origin:'+gid)
 # Verify native field vector converts back to the declared lab field for all four calls.
 for tag,fh in [('minus_h',-.001),('minus_half',-.0005),('plus_half',.0005),('plus_h',.001)]:
  jid=gid+'__'+tag;text=verify(r['rows'][jid]['input']['key']).read_text();fv=re.findall(r'^EXTERNAL-FIELD\s+(.*)$',text,re.M);ck(len(fv)==1,'one field:'+jid)
  for k,v in enumerate(fv[0].split()):close(D(v)/D(str(m['conversion']['FIELD_FACTOR'])),D(str(fh))*D(str(m['direction'][k])),jid+'/field')
 vals={}
 for tag,h in [('h',D('.001')),('half',D('.0005'))]:
  p,n=raws[gid+'__plus_'+tag],raws[gid+'__minus_'+tag];mu=-(D(p)-D(n))/(2*h*K);chi=-(D(p)+D(n)-2*D(zero))/(h*h*K);mb=(half(p)+half(n))/(2*h*K);cb=(half(p)+half(n)+2*half(zero))/(h*h*K)
  actual=s['states'][gid]['steps'][tag]
  for name,v in [('mu_E_au',mu),('chi_E_a0_cubed',chi),('mu_print_bound_au',mb),('chi_print_bound_a0_cubed',cb),('model_minus_QM_mu_au',mu-qms[gid])]:close(v,actual[name],gid+'/'+tag+'/'+name)
  vals[tag]=dict(mu=mu,chi=chi,mb=mb,cb=cb)
 states[gid]=vals;gates[gid+'/mu_field_step']=abs(vals['half']['mu']-vals['h']['mu'])<=D('.0001');gates[gid+'/chi_field_step']=abs(vals['half']['chi']-vals['h']['chi'])<=D('.5')
def state(metal,q):return next(k for k,v in m['rows'].items() if v['metal']==metal and v['q_A']==q)
paired={};derivatives={}
for q in [0.,-.01,-.005,.005,.01]:
 la,dy=state('La',q),state('Dy',q);paired[str(q)]={}
 for tag in ['h','half']:
  a,b=states[la][tag],states[dy][tag];v={'mu':b['mu']-a['mu'],'chi':b['chi']-a['chi'],'mb':a['mb']+b['mb'],'cb':a['cb']+b['cb']};paired[str(q)][tag]=v
  for old,new in [('mu','mu_E_au'),('chi','chi_E_a0_cubed'),('mb','mu_print_bound_au'),('cb','chi_print_bound_a0_cubed')]:close(v[old],s['paired_Dy_minus_La'][str(q)]['steps'][tag][new],str(q)+'/paired/'+tag+'/'+old)
 gates[f'paired/q={q}/chi_field_step']=abs(paired[str(q)]['half']['chi']-paired[str(q)]['h']['chi'])<=1
for kind,coeff in [('La',{'La':1}),('Dy',{'Dy':1}),('paired',{'La':-1,'Dy':1})]:
 derivatives[kind]={}
 for q in [.01,.005]:
  vals={}
  for tag in ['h','half']:
   v=b=qm=qb=D(0)
   for metal,c in coeff.items():
    p,n=state(metal,q),state(metal,-q);den=2*D(str(q));v+=c*(states[p][tag]['mu']-states[n][tag]['mu'])/den;b+=abs(c)*(states[p][tag]['mb']+states[n][tag]['mb'])/den;qm+=c*(qms[p]-qms[n])/den;qb+=abs(c)*sum((D(str(m['rows'][k]['quantum_projected_print_bound_au'])) for k in [p,n]),D(0))/den
   actual=s['rotation_dipole_derivatives'][kind][str(q)]['field_steps'][tag]
   for name,x in [('dmu_dq_ea0_per_A',v),('print_bound_ea0_per_A',b),('QM_dmu_dq_ea0_per_A',qm),('QM_print_bound_ea0_per_A',qb),('model_minus_QM_ea0_per_A',v-qm),('model_minus_QM_print_bound_ea0_per_A',b+qb)]:close(x,actual[name],kind+'/'+str(q)+'/'+tag+'/'+name)
   vals[tag]={'model':v,'QM':qm,'difference':v-qm,'print_bound':b,'QM_print_bound':qb}
  derivatives[kind][str(q)]=vals;gates[f'{kind}/q={q}/dmu_field_step']=abs(vals['half']['model']-vals['h']['model'])<=D('.0001')*sum(abs(c) for c in coeff.values())/D(str(q))
 for tag in ['h','half']:
  close(derivatives[kind]['0.005'][tag]['model']-derivatives[kind]['0.01'][tag]['model'],s['rotation_dipole_derivatives'][kind]['smaller_minus_larger_q_descriptive'][tag]['dmu_dq_ea0_per_A'],kind+'/q-step descriptive')
ck(gates==s['numerical_checks'] and len(gates)==31 and all(gates.values()),'independent31 numerical gates')
ck(r['central_gate']['passed'] and all(r['central_gate']['checks'].values()),'central8 gate passed before32')
ck(all(v['parsed']['native']['polarization']['printed_iteration_residual'] is None for v in r['rows'].values()),'no invented printed induction residual')
event=read(WB/'run_v1/completion_event.json')['records']['18764'];ck(event['state']=='COMPLETED' and event['exit_code']=='0:0','scheduler completed');cost=dict(r['cost'],scheduler_elapsed_seconds=int(event['elapsed_seconds']),scheduler_allocated_CPUs=int(event['allocated_cpus']),scheduler_allocated_CPU_seconds=int(event['elapsed_seconds'])*int(event['allocated_cpus']))
ck(abs(sum(rec['receipts'][0]['elapsed_seconds'] for rec in r['rows'].values())-cost['native_elapsed_seconds'])<1e-10,'native cost sum')
raw_replay={'states':states,'paired':paired,'derivatives':derivatives,'raw_field_energy_tokens':raws,'numerical_checks':gates}
# Decimal values retained as strings, with original producer floats bound separately.
def serial(x):
 if isinstance(x,D):return str(x)
 if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
 if isinstance(x,list):return [serial(v) for v in x]
 return x
save(OUT/'RAW_REPLAY.json',serial(raw_replay));save(OUT/'CHECKS.json',checks)
review={'status':'passed_dry_acetate_electrical_actual_review','implementation':r['implementation'],'manifest':r['manifest'],'result':pin(RPATH),'release':r['release'],'analyze':r['analyze'],'review_implementation':pin(__file__),'source_review':pin(WB/'independent_review_v1/REVIEW.json'),'raw_replay':pin(OUT/'RAW_REPLAY.json'),'checks':pin(OUT/'CHECKS.json'),'check_count':len(checks),'complete_cells':40,'zero_energies_reused':10,'quantum_origin_replays':10,'central_gate_passed':True,'numerical_checks_passed':True,'numerical_check_count':31,'maximum_independent_float_difference':max(diffs),'new_native_calls':0,'saved_admission_calls':1,'scientific_cost':cost,'review_cost':{'wall_seconds':time.monotonic()-t0,'process_CPU_seconds':time.process_time()-c0},'limitations':['Numerical gates are not physical agreement thresholds.','General dual-field energy derivative; no uind=uinp or uind-only moment assumption.','No new quantum susceptibility; saved zero-field quantum moments only.','POLAR-EPS declared1e-8 but analyze E does not print actual induction residual.','Two geometry steps are descriptive; no causality attribution.']}
save(OUT/'REVIEW.json',review);print(json.dumps({'review':pin(OUT/'REVIEW.json'),'checks':len(checks),'maxdiff':max(diffs),'cost':cost},indent=2))
