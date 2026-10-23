#!/usr/bin/env python3
"""One finite saved-native field parity panel; no native invocation interface."""
import argparse,hashlib,importlib.util,json,math,platform,re,sys,time
from decimal import Decimal
from pathlib import Path
import torch

def pin(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
def verify(p):
    assert pin(p['path'])==p,('changed input',p);return Path(p['path'])
def read(p):return json.loads(verify(p).read_text())
def load(p):
    spec=importlib.util.spec_from_file_location('cached_field_kernel',verify(p));mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod);return mod
def half(s):return float(Decimal(5).scaleb(Decimal(s.replace('D','E')).as_tuple().exponent-1))
def total_print(m):
    s=verify(m['analyze_output']).read_text();a=re.findall(r'Total Potential Energy\s*:\s*([-+\d.EeDd]+)',s);assert len(a)==1 and float(a[0])==m['energy_kcal_mol'];return half(a[0])
def ar(dim,*xs):return 128*dim*torch.finfo(torch.float64).eps*max(1.,*[abs(float(x)) for x in xs])
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--inputs',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(1);torch.set_num_interop_threads(1);tick=time.monotonic();cpu=time.process_time();cfg=json.loads(a.inputs.read_text());r={'status':'failed','configuration':pin(a.inputs),'implementation':pin(__file__),'kernel':cfg['kernel'],'rows':{},'susceptibility':{},'kernel_evaluations':0,'native_calls':0,'quantum_calls':0,'error':None,'physical_accuracy_qualified':False,'full_coordinate_forces_qualified':False}
    T=lambda x:torch.tensor(x,dtype=torch.float64)
    try:
        mod=load(cfg['kernel']);read(cfg['kernel_ready']);data=read(cfg['dataset']);compat=read(cfg['compatibility']);bridge=read(cfg['shared_bridge']);assert compat['dataset']==cfg['dataset'] and bridge['water_dataset']==cfg['dataset']
        units=read(data['baseline'])['units'];D=units['debye_per_eA'];K=units['electric'];B=units['bohr_A'];fields=read(data['fields']);exports=read(data['inputs']['export_result']);assert len(fields)==cfg['maximum_kernel_evaluations']==16
        perstate={}
        for id,row in fields.items():
            f=read(row['baseline_zero_operator']);op=f['native_operator'];geom=f['geometry'];mid=f['id'];record=exports['rows'][mid];q=record['admission']['qualification'] if record.get('inherited') else record['numerical'];snap=q['exports'][1];dim=len(op['Ed_e_per_A2']);n=dim//3;bounds=snap['numerical_bounds'];lmin=snap['spectrum']['A'][0];ar0=bounds['arithmetic'];an=math.sqrt(dim)*ar0
            Araw=T(op['A0_Aminus3']);A,sym=mod.symmetrize_native_operator(Araw,max_absolute_asymmetry=ar0,admission_authority={'fixture':row['baseline_zero_operator'],'native_qualification':data['inputs']['export_review']})
            ed,ep=T(op['Ed_e_per_A2']),T(op['Ep_e_per_A2']);assert torch.equal(ed,ep),'This finite panel needs existing equal-channel vacuum state'
            base=mod.Baseline(A,ed,ep,tuple(op['site_atom_indices']),K,T(op['nonpolar_energy_kcal']),T(op['permanent_mu_eA']),authority=f['baseline_id'])
            geometry=mod.Geometry(geom['metal_atom'],tuple(geom['atomic_numbers']),tuple(mod.Ligand.from_dict(x) for x in geom['ligands']))
            native=row['model_finite_field_record'];moment=native['moments'];Fau=T(native['job']['field_au']);F=mod.field_au_to_native(Fau,units)
            assert cfg['compatibility'] and compat['field_baseline_joins'][id]['same_zero_field_physics_and_geometry']
            assert float(Fau[2])==row['reference_cell']['field_au'] and float(Fau[0])==float(Fau[1])==0
            assert r['kernel_evaluations']<16;r['kernel_evaluations']+=1
            out=mod.evaluate(T(geom['xyz_A']),geometry,mod.Parameters(),base,F)
            tokens=[s.split() for s in verify(moment['induced_output']).read_text().splitlines() if s.strip()];assert int(tokens[0][0])==n and len(tokens)==n+1
            parsed=[];printbounds=[]
            for i,t in enumerate(tokens[1:]):
                assert int(t[0])==i+1;parsed.append([float(v) for v in t[-3:]]);printbounds.append([half(v) for v in t[-3:]])
            assert parsed==moment['induced_Debye'];native_u=T(parsed);pb=T(printbounds)
            raw_res=float(torch.linalg.vector_norm(Araw@out['u']-ed-F.repeat(n)))
            oldq=moment['induction']['maximum_printed_RMS_residual_Debye']
            native_bound=bounds['post_peek_amplification']*math.sqrt(n)*(oldq+.5e-10)/D+an
            solution_bound=(native_bound+raw_res+an)/lmin
            pred=torch.empty((n,3),dtype=torch.float64)
            for pos,atom in enumerate(op['site_atom_indices']):pred[atom]=D*out['u'][3*pos:3*pos+3]
            checks={}
            def check(name,error,limit):
                checks[name]={'error':float(error),'limit':float(limit),'passed':float(error)<=float(limit)}
                if not checks[name]['passed']:raise ValueError((id,name,checks[name]))
            check('native_physical_u_each_site_Debye',torch.max(torch.abs(pred-native_u)),torch.max(pb)+D*solution_bound)
            permanent=base.permanent_mu_eA*D;physical=permanent+pred.sum(0);native_physical=permanent+native_u.sum(0);conjugate=out['mu_eA']*D
            moment_bound=float(pb.sum(0).max())+D*math.sqrt(n)*solution_bound+ar(dim,*native_physical)
            check('physical_total_moment_from_u',torch.max(torch.abs(physical-native_physical)),moment_bound)
            check('conjugate_vs_physical_equal_RHS_only',torch.max(torch.abs(conjugate-physical)),ar(dim,*physical))
            check('saved_precise_total_moment',torch.max(torch.abs(conjugate-T(moment['total_precise_lab_Debye']))),moment_bound+ar(dim,*conjugate))
            energy_prints=total_print(f['baseline_observables'])+q['polarization_energy']['print_half_unit']+total_print(moment)
            epfield=ep+F.repeat(n);energy_bound=energy_prints+.5*K*float(torch.linalg.vector_norm(epfield))*solution_bound+ar(dim,out['energy'],moment['energy_kcal_mol'],op['nonpolar_energy_kcal'])
            check('native_field_total_energy_kcal_mol',abs(float(out['energy'])-moment['energy_kcal_mol']),energy_bound)
            chi=mod.response_to_source_au(out['mu_eA'],out['chi_A3'],units)['chi_source_au'][:,2]
            r['rows'][id]={'field_au':Fau.tolist(),'field_native_e_per_A2':F.tolist(),'fixture':row['baseline_zero_operator'],'native_uind':moment['induced_output'],'native_output':moment['analyze_output'],'symmetry_receipt':sym,'checks':checks,'native_physical_moment_Debye':native_physical.tolist(),'predicted_physical_moment_Debye':physical.tolist(),'predicted_conjugate_moment_Debye':conjugate.tolist(),'native_v_nonzero_field_available':False,'u_v_equality_scope':'exact equal baseline RHS plus same uniform field; not an unequal-channel qualification','solution_error_bound_eA':solution_bound,'CG_print_half_unit_Debye':.5e-10,'per_site_print_bounds_Debye':printbounds,'energy_total_printing_bound':energy_prints,'native_energy_kcal_mol':moment['energy_kcal_mol'],'predicted_energy_kcal_mol':float(out['energy']),'chi_source_au_column_z':chi.tolist()}
            perstate.setdefault(mid,{})[native['job']['id'].split('__')[-1]]={'native_u':native_u,'print':pb,'bound':solution_bound,'chi':chi,'dim':dim,'n':n,'f':float(Fau[2])}
        for mid,v in perstate.items():
            r['susceptibility'][mid]={}
            for tag in ('h','half'):
                plus,minus=v['plus_'+tag],v['minus_'+tag];h=plus['f'];assert h==-minus['f'] and h in (.001,.0005)
                observed=(plus['native_u']-minus['native_u']).sum(0)/(2*h*D*B)
                printbound=float((plus['print']+minus['print']).sum(0).max())/(2*h*D*B)
                limit=printbound+math.sqrt(plus['n'])*(plus['bound']+minus['bound'])/(2*h*B)+ar(plus['dim'],*plus['chi'])
                err=float(torch.max(torch.abs(plus['chi']-observed)));r['susceptibility'][mid][tag]={'native_finite_difference_column_source_au':observed.tolist(),'predicted_column_source_au':plus['chi'].tolist(),'error':err,'limit':limit,'passed':err<=limit}
                assert err<=limit,(mid,tag,r['susceptibility'][mid][tag])
        assert len(r['rows'])==16 and len(r['susceptibility'])==4;r['status']='passed_saved_native_field_parity'
    except Exception as exc:r['error']=repr(exc)
    r['cost']={'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1,'host':platform.node()};(a.out/'RESULT.json').write_text(json.dumps(r,indent=2,allow_nan=False)+'\n');print(json.dumps({'status':r['status'],'completed':len(r['rows']),'kernel_evaluations':r['kernel_evaluations'],'error':r['error'],'cost':r['cost']}));return 0 if r['error'] is None else 1
if __name__=='__main__':raise SystemExit(main())
