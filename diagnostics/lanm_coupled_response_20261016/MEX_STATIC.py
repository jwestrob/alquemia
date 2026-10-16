#!/usr/bin/env python3
"""Static2049-atom Mex additive model:2local MACE+2fixed-FF-charge ALPB calls.

CLI --config CONFIG --output FRESH; one GPU allocation, root submits.
Config pins: implementation,plan,fortran,field_result12570,solvent_force_result12566,
mace_force_result12500,native_helper(LBT3_FD.py),classical{La:pin,Dy:pin}.
No new SCC/MM evaluations, coordinate motion or optimization.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import socket
import subprocess
import sys
import time
import traceback

ROOT = Path('/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs')
HA = 627.509474
BOHR = .529177210903
EV = 23.06054783061903
GROUPS = ('retained_bonded','MM_Coulomb','MM_LJ','QM_MM_LJ')


def require(ok,message):
    if not ok: raise ValueError(message)


def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return {'path':str(p),'sha256':h.hexdigest()}


def verify(p):
    require(pin(p['path'])['sha256']==p['sha256'],'Changed pin: '+p['path'])
    return Path(p['path'])


def read(p): return json.loads(verify(p).read_text())


def write(p,d):
    with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')


def module(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def inputs(config,out):
    cfg=json.loads(config.read_text());verify(cfg['plan'])
    require(pin(__file__)['sha256']==cfg['implementation']['sha256'],'Implementation differs')
    field=read(cfg['field_result']);solvent=read(cfg['solvent_force_result']);mace=read(cfg['mace_force_result'])
    require(field['status']=='complete' and field['slurm_job_id']=='12570' and field['attempted']==field['accepted']==4,
            'Four accepted local field calls required')
    require(solvent['status']=='passed' and solvent['slurm_job_id']=='12566' and len(solvent['native_calls'])==26
            and all(c['pass'] for c in solvent['checks']),'Native fixed-charge force qualification differs')
    require(mace['status']=='complete' and mace['slurm_job_id']=='12500' and mace['declared_cells']==4
            and mace['analytic_forces_requested'] and mace['model_parameters_unchanged'],'Native MACE force precedent differs')
    fcfg=read(field['config']);scfg=read(solvent['config']);mcfg=read(mace['config'])
    require(pin(sys.executable)['sha256']==mcfg['python']['sha256'],'Pinned MACE Python differs')
    require(mace['model']==mcfg['model'] and mace['versions']==mcfg['expected_versions'],'Native MACE model/software differs')
    local=module('mex_static_field',verify(field['implementation']))
    world=local.inputs(verify(field['config']))
    np,mapper=world[5],world[8]
    prep=read(field['preparation_result']);source=read(prep['source_inputs'])
    require(prep['status']=='prepared' and prep['slurm_job_id']=='12569' and prep['molecular_calls']==0,
            'Prepared physical source unavailable')
    export=read(source['export']);native=read(export['atoms'])
    require(export['atom_count']==len(native)==2045 and abs(export['charge_e']+10)<=1e-10,'Native physical FF inventory')
    native_by={a['id']:a for a in native};require(len(native_by)==2045,'Duplicate FF source IDs')
    for p in [export['system'],export['bonds'],*export['forcefields']]:verify(p)
    by={(r['task']['metal'],r['task']['kind']):r for r in field['rows']}
    require(len(by)==4 and set(by)=={(m,k) for m in ('La','Dy') for k in ('isolated','embedded')},'Field cell denominator')
    for row in by.values():
        require(row['status']=='accepted','Unavailable field row')
        saved=read(row['result']);require(saved=={k:v for k,v in row.items() if k!='result'},'Saved field row differs')
        for p in row['artifacts'].values():verify(p)
        execution=read(row['execution']);d=Path(row['execution']['path']).parent
        replay=local.admit(d,row['task'],execution,world)
        for key in ('energy_hartree','charges_e','qm_gradient','pointcharge_gradient','mapped_physical_gradient'):
            require(replay[key]==row[key],'Saved field raw replay differs: '+key)
    cells={};signature=None
    for metal in ('La','Dy'):
        mp=read(prep['mappings'][metal]);physical=mp['physical_atoms'];ids=[a['id'] for a in physical]
        require(len(ids)==len(set(ids))==2049 and len(mp['core'])==209 and len(mp['caps'])==4 and len(mp['field'])==1840,
                'Physical/QM/field/cap inventory differs')
        target_rows=[a for a in mp['core'] if a['kind']=='metal'];require(len(target_rows)==1,'One target required')
        target=mapper.target_source_id(target_rows[0]);mi=target_rows[0]['qm_index']
        ions={a['id']:a for a in source['spectator_occupancy']}
        require(len(ions)==3 and all(a['element']=='Nd' and a['charge_e']==3 for a in ions.values()),'Spectator state differs')
        require(set(ids)-set(native_by)==set(ions)|{target},'Physical FF/ion source join differs')
        q=[]
        for a in physical:
            if a['id'] in native_by:
                n=native_by[a['id']];require(n['element']==a['element'],'Physical FF element differs');q.append(n['charge_e'])
            else:
                expected=metal if a['id']==target else 'Nd'
                require(a['element']==expected,'Physical ion identity differs');q.append(3.)
        q=np.asarray(q,dtype=float);x=np.asarray([a['xyz_A'] for a in physical],dtype=float)
        require(np.isfinite(q).all() and np.isfinite(x).all() and abs(math.fsum(q)-2)<=1e-5,'Fixed FF charge closure/coordinates')
        sig=(ids,q.tolist(),x.tolist())
        if signature is None:signature=sig
        else:require(sig==signature,'Paired physical geometry/charge vector differs')
        symbols,qp=local.xyz(source['configurations']['A']['endpoints'][metal]['xyz'])
        require(mi==0 and symbols[mi]==metal and symbols.count(metal)==1,'Local target/source order')
        classical=read(cfg['classical'][metal]);ledger=read(classical['ledger'])
        require(classical['inputs']==prep['source_inputs']==ledger['inputs'] and classical['metal']==metal
                and classical['quantity']=='gradient' and classical['units']=='kcal/mol/angstrom'
                and classical['QM_MM_Coulomb_included'] is False and classical['C4_induction_included'] is False,
                'Classical source/state/term exclusions differ')
        require(ledger['boundary_policy']=='fixed_cap_radial_bonds_retained_substitutable_angular_terms_omitted_v1'
                and ledger['native_system']==export['system'],'Classical boundary/native source differs')
        verify(ledger['implementation']);verify(ledger['reused_builder'])
        for p in ledger['artifacts'].values():verify(p)
        for v in ledger['ion_parameters'].values():verify(v['file'])
        ci={sid:i for i,sid in enumerate(classical['source_ids'])}
        require(len(ci)==2049 and set(ci)==set(ids),'Classical physical IDs differ')
        order=[ci[sid] for sid in ids];cg=np.asarray(classical['gradient_kcal_mol_A'],dtype=float)[order]
        require(np.allclose(np.asarray(classical['coordinates_A'])[order],x,atol=1e-8,rtol=0)
                and cg.shape==(2049,3) and np.isfinite(cg).all(),'Classical geometry/gradient differs')
        require(set(classical['components_kcal_mol'])==set(GROUPS)
                and set(classical['component_gradients_kcal_mol_A'])==set(GROUPS),'Classical component inventory')
        require(abs(math.fsum(classical['components_kcal_mol'].values())-classical['total_kcal_mol'])<=1e-8
                and np.allclose(sum(np.asarray(classical['component_gradients_kcal_mol_A'][k])[order] for k in GROUPS),cg,atol=1e-10,rtol=1e-12),
                'Saved classical component sum differs')
        a,b=by[(metal,'embedded')],by[(metal,'isolated')]
        delta_e=math.fsum([a['energy_hartree'],-b['energy_hartree']])
        delta_g=np.asarray(a['mapped_physical_gradient']['gradient_kcal_mol_A'])-np.asarray(b['mapped_physical_gradient']['gradient_kcal_mol_A'])
        recorded=next(r for r in field['field_differences'] if r['metal']==metal)
        require(recorded['status']=='available' and recorded['energy_hartree']==delta_e
                and np.array_equal(recorded['gradient_kcal_mol_A'],delta_g),'Saved field difference reconstruction')
        cells[metal]={'mapping':mp,'mapping_pin':prep['mappings'][metal],'physical':physical,'ids':ids,'x':x,'q':q,
                      'symbols':symbols,'qm_positions':np.asarray(qp),'target_index':mi,'target_id':target,
                      'xyz':source['configurations']['A']['endpoints'][metal]['xyz'],
                      'classical':classical,'classical_gradient':cg,'classical_order':order,
                      'field_energy_hartree':delta_e,'field_gradient':delta_g,
                      'field_rows':{k:by[(metal,k)]['result'] for k in ('embedded','isolated')}}
    reference=verify(scfg['fortran']).read_text()
    expected=reference.replace('dq = 5.0e-9_wp ! half the archived eight-decimal charge unit',
           'dq = 0.0_wp ! exact declared FF model charges, serialized with round-trip precision').replace(
           "if (n /= 256) error stop 'qualified fixture scope requires 256 atoms'",
           "if (n /= 2049) error stop 'declared Mex physical scope requires 2049 atoms'")
    require(verify(cfg['fortran']).read_text()==expected,'Unexpected native wrapper functional change')
    for p in scfg['build_inputs']:verify(p)
    previous=read(scfg['previous_build']);require(previous['status']=='built','Native static build unavailable')
    for p in previous['linked_libraries']:verify(p)
    helper=module('mex_static_native',verify(cfg['native_helper']))
    require(helper.EV_TO_KCAL==EV and helper.HA_TO_KCAL==HA,'Native helper unit convention differs')
    write(out/'PRECHECK.json',{'status':'passed','source_inputs':prep['source_inputs'],'physical_atoms':2049,
          'QM_atoms':209,'caps':4,'field_sites':1840,'fixed_charge_e':2,'charge_rule':'unshifted ff19SB/TIP3P plus formal3 ions',
          'MACE_calls':0,'ALPB_calls':0,'field_result':cfg['field_result'],'saved_classical':cfg['classical']})
    return cfg,cells,scfg,mcfg,helper,mapper,np


def mace_cell(native,c,metal):
    torch=native.torch;symbols=c['symbols'];positions=c['qm_positions'];mult=1 if metal=='La' else 6
    require(len(symbols)==209 and positions.shape==(209,3) and symbols[0]==metal,'Local MACE inventory')
    require((sum(native.Z[s] for s in symbols)-mult+1)%2==0,'Local physical spin parity')
    atoms=native.Atoms(symbols,positions=positions,pbc=False);atoms.info.update(charge=0,spin=mult)
    batch=native.calc._atoms_to_batch(atoms)
    state=native.input_batch(native.calc,atoms,0,mult,batch,0);tick=time.monotonic()
    with torch.enable_grad():
        native.calls+=1
        value=native.model(batch.to_dict(),training=False,compute_force=True,compute_virials=False,compute_stress=False,
                          compute_displacement=False,compute_hessian=False,compute_edge_forces=False,compute_atomic_stresses=False)
    energy=float(value['energy'].detach().cpu().item());force=value['forces'].detach().cpu()
    require(math.isfinite(energy) and tuple(force.shape)==(209,3) and bool(torch.isfinite(force).all()),'Invalid local native MACE result')
    require([p._version for p in native.model.parameters()]==native.parameter_versions,'Native model mutated')
    return {'energy_eV':energy,'forces_eV_A':force.tolist(),'input_state_check':state,'xyz':c['xyz'],
            'physical_multiplicity':mult,'charge':0,'wall_seconds':time.monotonic()-tick,'force_convention':'negative vacuum energy derivative'}


def solvent_cell(binary,c,d,env,np):
    d.mkdir();xp=d/'physical.xyz';qp=d/'charges'
    xp.write_text('2049\nExact repaired physical Mex; fixed unshifted FF plus formal3 charges\n'+''.join(
        '%s %.17g %.17g %.17g\n'%(a['element'],*v) for a,v in zip(c['physical'],c['x'])))
    qp.write_text(''.join('%.17g\n'%q for q in c['q']))
    require(np.array_equal(np.loadtxt(qp),c['q']),'FF charge serialization changed declared vector')
    command=[str(binary),str(xp),str(qp)];start=time.monotonic()
    write(d/'START.json',{'command':command,'xyz':pin(xp),'charges':pin(qp),'physical_mapping':c['mapping_pin']})
    with (d/'stdout.txt').open('x') as so,(d/'stderr.txt').open('x') as se:
        proc=subprocess.run(command,cwd=d,env=env,stdout=so,stderr=se)
    execution={'command':command,'returncode':proc.returncode,'wall_seconds':time.monotonic()-start}
    write(d/'EXECUTION.json',execution);require(proc.returncode==0,'Native whole-solvent execution failed')
    text=(d/'stdout.txt').read_text()
    def values(label,n):
        found=re.findall(r'^'+label+r'\s+(.+)$',text,re.M);require(len(found)==1,'Missing/duplicate '+label)
        v=[float(x) for x in found[0].split()];require(len(v)==n and np.isfinite(v).all(),'Invalid '+label);return v
    parts=dict(zip(('Gelec','Ghb','Gsasa','Gshift'),values('COMPONENTS',4)))
    charge=values('CHARGE_SUM',1)[0];require(abs(charge-2)<=1e-5,'Native global charge closure')
    require(values('ROUNDING_BOUNDS',2)==[0.,0.],'Unexpected archived-GFN charge uncertainty assumption')
    g=np.loadtxt(d/'gradient.dat');require(g.shape==(2049,3) and np.isfinite(g).all(),'Invalid full physical solvent gradient')
    result={'components_hartree':parts,'total_hartree':math.fsum(parts.values()),'charge_e':charge,
            'gradient_hartree_per_bohr':g.tolist(),'gradient':pin(d/'gradient.dat'),'xyz':pin(xp),'charges':pin(qp),
            'fixed_charge_rule':'unmodified ff19SB/TIP3P physical charges plus four formal3 ions; no cap population projection',
            'execution':pin(d/'EXECUTION.json'),'cavity_global':values('CAVITY_GLOBAL',4)}
    write(d/'RESULT.json',result);return result


def run(config,out,result):
    cfg,cells,scfg,mcfg,helper,mapper,np=inputs(config,out)
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',
             NUMEXPR_NUM_THREADS='1',LD_LIBRARY_PATH=scfg['library_directory'])
    binary=out/'mex_fixed_solvent';command=scfg['compiler_command']+[str(verify(cfg['fortran'])),'-o',str(binary)]+scfg['link_arguments']
    with (out/'compile.out').open('x') as so,(out/'compile.err').open('x') as se:
        proc=subprocess.run(command,cwd=out,env=env,stdout=so,stderr=se)
    result['build']={'command':command,'returncode':proc.returncode,'fortran':cfg['fortran'],'parent':cfg['solvent_force_result']}
    require(proc.returncode==0,'Thin native wrapper compile failed');result['build']['binary']=pin(binary)
    native=helper.NativeMACE(mcfg);rows=[];result['rows']=rows
    for metal in ('La','Dy'):
        c=cells[metal];d=out/metal;d.mkdir();row={'metal':metal,'status':'unavailable'};rows.append(row)
        try:
            result['new_MACE_calls_attempted']+=1
            ml=mace_cell(native,c,metal);write(d/'MACE.json',ml);row['MACE']=pin(d/'MACE.json')
            result['new_ALPB_calls_attempted']+=1
            sol=solvent_cell(binary,c,d/'solvent',env,np);row['ALPB']=pin(d/'solvent/RESULT.json')
            qm_g=-np.asarray(ml['forces_eV_A'])*EV*BOHR/HA
            mapped=mapper.scatter(c['mapping'],qm_g,np.zeros((1840,3)))
            component_g={'MACE_vacuum':np.asarray(mapped['gradient_kcal_mol_A']),
                 'GFN2_field_difference':c['field_gradient'],'MM_ledger':c['classical_gradient'],
                 'fixed_FF_charge_ALPB':np.asarray(sol['gradient_hartree_per_bohr'])*HA/BOHR}
            energies={'MACE_vacuum':ml['energy_eV']*EV,'GFN2_field_difference':c['field_energy_hartree']*HA,
                      'MM_ledger':c['classical']['total_kcal_mol'],'fixed_FF_charge_ALPB':sol['total_hartree']*HA}
            total_g=sum(component_g.values());require(total_g.shape==(2049,3) and np.isfinite(total_g).all(),'Invalid total gradient')
            artifact={'metal':metal,'physical_atoms':c['physical'],'source_ids':c['ids'],'physical_mapping':c['mapping_pin'],
                 'source_xyz':c['xyz'],'component_energies_kcal_mol':energies,'total_energy_kcal_mol':math.fsum(energies.values()),
                 'component_gradients_kcal_mol_A':{k:v.tolist() for k,v in component_g.items()},
                 'gradient_kcal_mol_A':total_g.tolist(),'force_kcal_mol_A':(-total_g).tolist(),
                 'translation_gradient_sum_kcal_mol_A':total_g.sum(axis=0).tolist(),
                 'MACE_cap_mapping':mapped,'MACE':row['MACE'],'ALPB':row['ALPB'],'GFN2_reused':c['field_rows'],
                 'classical_reused':cfg['classical'][metal],'total_derivative_FD_qualified':False,
                 'solvent_charge_response_included':False,'optimizer_steps':0,'affinity':None}
            write(d/'ASSEMBLED.json',artifact);row.update(status='complete',assembled=pin(d/'ASSEMBLED.json'),
                  component_energies_kcal_mol=energies,total_energy_kcal_mol=artifact['total_energy_kcal_mol'])
        except Exception as exc:row.update(failure=repr(exc),traceback=traceback.format_exc())
        write(d/'RESULT.json',row);row['result']=pin(d/'RESULT.json')
    result.update(status='complete' if all(r['status']=='complete' for r in rows) else 'unavailable',
                  native_model_calls=native.calls,native_device=native.device,model=mcfg['model'],versions=native.versions)
    if result['status']=='complete':
        result['raw_origin_Dy_minus_La_kcal_mol']=rows[1]['total_energy_kcal_mol']-rows[0]['total_energy_kcal_mol']
        result['component_Dy_minus_La_kcal_mol']={k:rows[1]['component_energies_kcal_mol'][k]-rows[0]['component_energies_kcal_mol'][k] for k in rows[0]['component_energies_kcal_mol']}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    cfg=json.loads(a.config.read_text());cpus=int(os.environ['SLURM_CPUS_PER_TASK'])
    result={'protocol_id':'mex_static_local_MACE_field_difference_MM_fixedFF_ALPB_v1','status':'running',
            'config':pin(a.config),'implementation':pin(__file__),'plan':cfg['plan'],'slurm_job_id':os.environ['SLURM_JOB_ID'],
            'host':socket.gethostname(),'allocated_cpus':cpus,'declared_metals':['La','Dy'],
            'new_MACE_calls_attempted':0,'new_ALPB_calls_attempted':0,'new_SCC_calls':0,'new_MM_calls':0,
            'reused_GFN2_cells':4,'reused_MM_cells':2,'optimizer_steps':0,
            'scope':'EF3 target at fixed three-Nd spectators; static conditional model assembly; no total derivative qualification, protein affinity or classifier'}
    write(out/'START.json',result)
    try:run(a.config.resolve(),out,result)
    except Exception as exc:result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc())
    result['elapsed_seconds']=time.monotonic()-start;result['application_allocated_CPU_seconds']=result['elapsed_seconds']*cpus
    result['application_allocated_GPU_seconds']=result['elapsed_seconds']
    write(out/'RESULT.json',result);print(json.dumps({'status':result['status'],'result':str(out/'RESULT.json')}),flush=True)
    return result['status']=='complete'


if __name__=='__main__':
    if not main():raise SystemExit(1)
