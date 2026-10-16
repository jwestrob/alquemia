#!/usr/bin/env python3
"""Four native fixed-FF ALPB calls; reuse MACE12500 and Mex12685. Root submits.
CLI: --config CONFIG --output FRESH. No imported molecular engine or new SCC/MM.
V2 repairs only the native-key/deposited-value source-map join after zero-call12711.
"""
import argparse
import hashlib
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
import xml.etree.ElementTree as ET

SOURCES=('7CCO','7CCN')
METALS=('La','Dy')
NAMES=('MACE_vacuum','GFN2_field_difference','MM_ledger','fixed_FF_charge_ALPB')
PARTS=('Gelec','Ghb','Gsasa','Gshift')
EV=23.06054783061903
HA=627.509474
TERMINI='free N-terminal NH3+ / C-terminal COO-; Asp3/7/13 and Glu11/14 deprotonated'


def require(ok,message):
    if not ok:raise ValueError(message)


def pin(path):
    p=Path(path).resolve();h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return {'path':str(p),'sha256':h.hexdigest()}


def verify(p):
    require(pin(p['path'])['sha256']==p['sha256'],'Changed pin: '+p['path']);return Path(p['path'])


def read(p):return json.loads(verify(p).read_text())


def write(p,value):
    with Path(p).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def pins_in(value):
    if isinstance(value,dict):
        if 'path' in value and 'sha256' in value:verify(value)
        else:
            for v in value.values():pins_in(v)
    elif isinstance(value,list):
        for v in value:pins_in(v)


def vectors(value,n):
    require(len(value)==n and all(len(v)==3 and all(math.isfinite(x) for x in v) for v in value),'Invalid coordinates/gradient')
    return value


def coordinates_equal(a,b):
    return len(a)==len(b) and all(abs(x-y)<=1e-12 for u,v in zip(a,b) for x,y in zip(u,v))


def xyz(p):
    lines=verify(p).read_text().splitlines();rows=[v.split() for v in lines[2:] if v.strip()]
    require(int(lines[0])==len(rows) and all(len(r)==4 for r in rows),'XYZ atom layout')
    return [r[0] for r in rows],vectors([[float(x) for x in r[1:]] for r in rows],len(rows))


def state(value,n,q,m,index):
    require(all(value[k]==v for k,v in {'atoms':n,'charge':q,'spin_multiplicity':m,
        'selected_metal_index':index,'graph_count':1,'head':'omol','charge_prediction_claimed':False}.items()),'Native MACE state differs')


def native_components(value):
    return {'MACE_vacuum':value['MACE_eV']*EV,'GFN2_field_difference':value['field_Eh']*HA,
            'MM_ledger':value['MM_kcal_mol'],'fixed_FF_charge_ALPB':value['ALPB_Eh']*HA}


def difference(a,b):
    return {'MACE_vacuum':(a['MACE_eV']-b['MACE_eV'])*EV,
            'GFN2_field_difference':(a['field_Eh']-b['field_Eh'])*HA,
            'MM_ledger':a['MM_kcal_mol']-b['MM_kcal_mol'],
            'fixed_FF_charge_ALPB':(a['ALPB_Eh']-b['ALPB_Eh'])*HA}


def mex_replay(result,cfg):
    require(result['status']=='complete' and result['slurm_job_id']=='12685'
        and result['new_MACE_calls_attempted']==result['new_ALPB_calls_attempted']==2
        and result['new_SCC_calls']==result['new_MM_calls']==0,'Admitted Mex parent differs')
    require(len(result['rows'])==2 and {r['metal'] for r in result['rows']}==set(METALS),'Mex denominator')
    values={};records={};paired=None
    for r in result['rows']:
        require(r['status']=='complete','Unavailable Mex cell');a=read(r['assembled']);metal=r['metal'];mult=1 if metal=='La' else 6
        require(read(r['result'])=={k:v for k,v in r.items() if k!='result'},'Mex saved row differs')
        pins_in(a);mp=read(a['physical_mapping']);physical=a['physical_atoms'];ids=[x['id'] for x in physical]
        require(len(ids)==len(set(ids))==2049 and ids==a['source_ids'] and physical==mp['physical_atoms']
            and len(mp['core'])==209 and len(mp['caps'])==4 and len(mp['field'])==1840,'Mex physical/map binding differs')
        require(a['MACE']==r['MACE'] and a['ALPB']==r['ALPB'] and a['metal']==mp['metal']==metal,'Mex component/state binding')
        ml=read(a['MACE']);sol=read(a['ALPB']);mm=read(a['classical_reused'])
        require(a['classical_reused']==cfg['classical'][metal] and ml['xyz']==a['source_xyz']==mp['xyz']
            and ml['charge']==0 and ml['physical_multiplicity']==mult,'Mex local source/charge/spin differs')
        state(ml['input_state_check'],209,0,mult,0);vectors(ml['forces_eV_A'],209)
        g={k:read(p) for k,p in a['GFN2_reused'].items()};require(set(g)=={'embedded','isolated'},'Mex field pair')
        for kind,row in g.items():
            pins_in(row['artifacts']);pins_in(row['task']);verify(row['execution']);t=row['task']
            require(row['status']=='accepted' and t['kind']==kind and t['metal']==metal and t['charge']==0
                and t['uhf']==0 and t['multiplicity']==1 and t['physical_multiplicity']==mult
                and t['mapping']==a['physical_mapping'] and t['source_xyz']==a['source_xyz'],'Mex saved GFN state/source differs')
        require(mm['metal']==metal and mm['inputs']==mp['inputs'] and not mm['QM_MM_Coulomb_included']
            and not mm['C4_induction_included'],'Mex MM ledger differs');pins_in(read(mm['ledger']))
        require(abs(math.fsum(mm['components_kcal_mol'].values())-mm['total_kcal_mol'])<=1e-8,'Mex MM component sum')
        for p in (sol['xyz'],sol['charges'],sol['gradient'],sol['execution']):verify(p)
        require(set(sol['components_hartree'])==set(PARTS) and math.fsum(sol['components_hartree'].values())==sol['total_hartree']
            and abs(sol['charge_e']-2)<=1e-5,'Mex ALPB components/charge differ')
        inp=read(mp['inputs']);export=read(inp['export']);ff={x['id']:x for x in read(export['atoms'])}
        ions={x['id'] for x in inp['spectator_occupancy']};target=next(x['id'] for x in physical if x['element']==metal)
        require(len(ff)==2045 and len(ions)==3 and set(ids)-set(ff)==ions|{target},'Mex native FF/ion join')
        q=[ff[s]['charge_e'] if s in ff else 3. for s in ids]
        require([float(v) for v in verify(sol['charges']).read_text().split()]==q,'Mex fixed FF charge vector differs')
        sy,x=xyz(sol['xyz']);require(sy==[p['element'] for p in physical] and coordinates_equal(x,[p['xyz_A'] for p in physical]),'Mex solvent geometry differs')
        sig=(ids,q,x)
        if paired is None:paired=sig
        else:require(sig==paired,'Mex paired physical geometry/charge differs')
        v={'MACE_eV':ml['energy_eV'],'field_Eh':math.fsum([g['embedded']['energy_hartree'],-g['isolated']['energy_hartree']]),
           'MM_kcal_mol':mm['total_kcal_mol'],'ALPB_Eh':sol['total_hartree'],'ALPB_parts_Eh':sol['components_hartree']}
        comp=native_components(v);require(all(math.isfinite(x) for x in comp.values()) and comp==a['component_energies_kcal_mol']==r['component_energies_kcal_mol']
            and math.fsum(comp.values())==a['total_energy_kcal_mol']==r['total_energy_kcal_mol'],'Mex raw energy reconstruction differs')
        values[metal]=v;records[metal]={'assembled':r['assembled'],'raw_MACE':a['MACE'],'raw_ALPB':a['ALPB'],
            'raw_GFN2':a['GFN2_reused'],'raw_MM':a['classical_reused'],'mapping':a['physical_mapping'],'exact_energy_replay':True}
    return values,records


def inputs(cfg,out):
    verify(cfg['plan']);verify(cfg['implementation']);require(pin(__file__)['sha256']==cfg['implementation']['sha256'],'Implementation differs')
    verify(cfg['recovery_note']);previous=read(cfg['previous_failure']);oldcfg=read(previous['config'])
    require(previous['status']=='execution_failure' and previous['slurm_job_id']=='12711'
        and previous['new_ALPB_calls_attempted']==previous['new_MACE_calls']==previous['new_GFN2_calls']==previous['new_MM_calls']==previous['new_DFT_calls']==0
        and len(previous['rows'])==4 and all(not r['attempted'] for r in previous['rows']),'Expected unchanged zero-call predecessor')
    verify(previous['implementation']);verify(oldcfg['implementation'])
    require({k:cfg[k] for k in oldcfg if k!='implementation'}=={k:v for k,v in oldcfg.items() if k!='implementation'},'Scientific configuration changed during source-map recovery')
    require(cfg['sources']==list(SOURCES) and cfg['metals']==list(METALS) and cfg['declared_ALPB_calls']==4,'Finite scope differs')
    prep=read(cfg['preparation']);mace=read(cfg['mace_result']);mex=read(cfg['mex_result'])
    require(prep['status']=='complete' and prep['slurm_job_id']=='12469' and len(prep['rows'])==2,'Preparation parent differs')
    require(mace['status']=='complete' and mace['slurm_job_id']=='12500' and mace['declared_cells']==len(mace['rows'])==4
        and mace['model_parameters_unchanged'] and mace['analytic_forces_requested'],'Four native MACE records required')
    pcfg=read(prep['config']);mcfg=read(mace['config']);xcfg=read(mex['config'])
    verify(prep['implementation']);verify(prep['plan']);verify(pcfg['plan'])
    require(prep['plan']==pcfg['plan'],'Preparation plan differs')
    for parent,c in ((mex,xcfg),):
        verify(parent['implementation']);verify(c['implementation']);verify(parent['plan']);verify(c['plan'])
        require(parent['implementation']['sha256']==c['implementation']['sha256'],'Parent source aliases differ')
    verify(mcfg['implementation']);verify(mcfg['plan']);verify(mcfg['model']['checkpoint'])
    require(pin(sys.executable)['sha256']==mcfg['python']['sha256'],'Pinned Python differs')
    require(mcfg['preparation']==cfg['preparation'] and xcfg['mace_force_result']==cfg['mace_result']
        and mace['model']==mex['model']==mcfg['model'] and mace['versions']==mex['versions']==mcfg['expected_versions'],'Shared MACE recipe differs')
    solvent=read(xcfg['solvent_force_result']);scfg=read(solvent['config']);build=read(scfg['previous_build'])
    require(solvent['status']=='passed' and solvent['slurm_job_id']=='12566' and len(solvent['native_calls'])==26
        and all(c['pass'] for c in solvent['checks']) and build['status']=='built','Native force/build precedent differs')
    for p in scfg['build_inputs']+build['linked_libraries']:verify(p)
    reference=verify(scfg['fortran']).read_text();source=verify(xcfg['fortran']).read_text()
    require(source==reference.replace('dq = 5.0e-9_wp ! half the archived eight-decimal charge unit',
        'dq = 0.0_wp ! exact declared FF model charges, serialized with round-trip precision').replace(
        "if (n /= 256) error stop 'qualified fixture scope requires 256 atoms'",
        "if (n /= 2049) error stop 'declared Mex physical scope requires 2049 atoms'"),'Mex native functional changed')
    marker="if (n /= 2049) error stop 'declared Mex physical scope requires 2049 atoms'"
    require(source.count(marker)==1 and mex['build']['fortran']==xcfg['fortran'],'Native wrapper lineage differs')
    native_source=source.replace(marker,"if (n /= 256) error stop 'declared LBT physical scope requires 256 atoms'")
    mexvalues,mexrecords=mex_replay(mex,xcfg);cells=[];signatures=[];proof=[]
    by={r['source_id']:r for r in prep['rows']};mr={(r['source_id'],r['metal']):r for r in mace['rows']}
    require(set(by)==set(SOURCES) and set(mr)=={(s,m) for s in SOURCES for m in METALS},'Source/metal denominator differs')
    for sid in SOURCES:
        pr=by[sid];p=read(pr['preparation']);physical=read(p['physical_atoms']);ff=read(p['native_peptide_atoms']);bonds=read(p['covalent_bonds'])
        ids=[a['id'] for a in physical];x=vectors([a['xyz_A'] for a in physical],256)
        require(pr['status']==p['status']=='prepared' and p['source_id']==sid and p['model']==1
            and p['source_conditioned_metal']==('La' if sid=='7CCO' else 'Lu') and p['sequence']=='FIDTNNDGWIEGDELLA'
            and p['atoms']==256 and p['protein_atoms']==255 and p['protein_charge']==-5 and p['charge']==-2
            and p['terminal_chemistry']==TERMINI and p['caps']==p['waters']==p['spectators']==[],'LBT chemistry/partition differs')
        require(len(ids)==len(set(ids))==256 and len(ff)==255 and ids[-1]=='metal' and physical[-1]['element']=='M'
            and p['deposited_H_coordinates_preserved'] and p['deposited_heavy_coordinates_preserved']
            and p['target_has_no_covalent_graph_edges'] and p['peptide_bond_count']==16,'LBT physical inventory differs')
        require({a['id'] for a in p['terminal_additions']}=={'A/1//H2','A/1//H3','A/17//OXT'},'Terminal additions differ')
        for key in ('source_pins','forcefield','completed_peptide','hydrogen_audit','native_peptide_system'):pins_in(p[key])
        require([a['id'] for a in ff]==ids[:-1] and all({k:a[k] for k in ('id','element','resname','xyz_A')}
            =={k:b[k] for k in ('id','element','resname','xyz_A')} for a,b in zip(ff,physical)),'Native FF/physical atom join differs')
        byid={a['id']:a for a in physical}
        require(all(sid in byid and byid[sid]['deposited_source']==a and a['element']==byid[sid]['element']
            and coordinates_equal([a['xyz_A']],[byid[sid]['xyz_A']]) for sid,a in p['source_to_native'].items()),'Deposited source mapping differs')
        tree=ET.parse(verify(p['native_peptide_system']));nb=tree.findall("./Forces/Force[@type='NonbondedForce']/Particles/Particle")
        require(len(tree.findall('./Particles/Particle'))==len(nb)==255,'Native FF XML inventory')
        q=[a['charge_e'] for a in ff]+[3.]
        require(all(math.isfinite(v) for v in q) and [float(a.attrib['q']) for a in nb]==q[:-1]
            and abs(math.fsum(q[:-1])+5)<=1e-6,'Native FF XML charge replay/closure')
        graph=[]
        for b in bonds:
            i,j=b['native_particle_indices'];require(0<=i<255 and 0<=j<255 and ids[i]==b['atom_a_id'] and ids[j]==b['atom_b_id'],'Covalent graph identity')
            graph.append({k:v for k,v in b.items() if k!='distance_A'})
        signatures.append({'identity':[(a['id'],a['element'],a.get('resname')) for a in physical],
            'q':q,'system_sha256':p['native_peptide_system']['sha256'],'graph':graph,'sequence':p['sequence'],'termini':TERMINI,'partition':'whole256'})
        proof.append({'source_id':sid,'preparation':pr['preparation'],'physical':p['physical_atoms'],
            'FF_atoms':p['native_peptide_atoms'],'FF_system':p['native_peptide_system'],'charge_sum_e':math.fsum(q)})
        paired=None
        for metal in METALS:
            t=p['endpoints']['whole'][metal];r=mr[(sid,metal)];mult=1 if metal=='La' else 6
            require(t in pr['tasks'] and all(r[k]==v for k,v in t.items()) and r['status']=='complete' and r['call_attempted']
                and any(all(task.get(k)==v for k,v in t.items()) and task['cache_key']==r['cache_key'] for task in mcfg['tasks']), 'Saved MACE row/task binding')
            require(t['atom_count']==256 and t['charge']==-2 and t['metal_index']==255 and t['multiplicity']==t['spin_multiplicity']==mult
                and t['oxidation_state']==3 and t['source_cap_map'] is None and t['physical_source']==p['physical_atoms'],'LBT electronic state differs')
            state(r['input_state_check'],256,-2,mult,255);vectors(r['forces_eV_A'],256);require(math.isfinite(r['energy_eV']),'Nonfinite MACE energy')
            sy,pos=xyz(t['xyz']);require(sy==[a['element'] for a in physical[:-1]]+[metal] and coordinates_equal(pos,x),'MACE/physical geometry differs')
            if paired is None:paired=pos
            else:require(pos==paired,'LBT paired geometry differs')
            cells.append({'source_id':sid,'metal':metal,'task':t,'physical':physical,'x':pos,'symbols':sy,'q':q,'MACE_eV':r['energy_eV']})
    pool_ok=signatures[0]==signatures[1]
    preflight={'status':'passed','declared_cells':4,'molecular_calls':0,'sources':proof,'Mex_replay':mexrecords,
        'common_two_source_pool_identity_proved':pool_ok,'pool_rule':'same physical membership/FF q/topology/terminal/electronic model; distinct source coordinates',
        'MACE_result':cfg['mace_result'],'Mex_result':cfg['mex_result'],'native_functional_parent':xcfg['solvent_force_result'],
        'native_wrapper_parent':xcfg['fortran'],'LBT_partition_zeros':['GFN2_field_difference','MM_ledger']}
    write(out/'PREFLIGHT.json',preflight)
    return cells,mexvalues,scfg,build,native_source,pool_ok


def solvent(binary,c,d,env):
    d.mkdir();xp=d/'physical.xyz';qp=d/'charges'
    xp.write_text('256\nWhole LBT3 fixed native FF plus formal3 target\n'+''.join('%s %.17g %.17g %.17g\n'%(s,*v) for s,v in zip(c['symbols'],c['x'])))
    qp.write_text(''.join('%.17g\n'%v for v in c['q']))
    require(xyz(pin(xp))==(c['symbols'],c['x']) and [float(v) for v in qp.read_text().split()]==c['q'],'Serialization changed declared input')
    command=[str(verify(binary)),str(xp),str(qp)];start=time.monotonic()
    write(d/'START.json',{'command':command,'xyz':pin(xp),'charges':pin(qp),'physical_source':c['task']['physical_source']})
    with (d/'stdout.txt').open('x') as so,(d/'stderr.txt').open('x') as se:proc=subprocess.run(command,cwd=d,env=env,stdout=so,stderr=se)
    write(d/'EXECUTION.json',{'command':command,'returncode':proc.returncode,'elapsed_seconds':time.monotonic()-start})
    require(proc.returncode==0,'Native ALPB call failed');text=(d/'stdout.txt').read_text()
    def values(label,n):
        matches=re.findall(r'^'+label+r'\s+(.+)$',text,re.M);require(len(matches)==1,'Native output field '+label)
        v=[float(x) for x in matches[0].split()];require(len(v)==n and all(math.isfinite(x) for x in v),'Invalid native field '+label);return v
    parts=dict(zip(PARTS,values('COMPONENTS',4)));charge=values('CHARGE_SUM',1)[0]
    require(abs(charge+2)<=1e-5 and values('ROUNDING_BOUNDS',2)==[0.,0.],'Native charge/rounding bookkeeping differs')
    vectors([[float(v) for v in line.split()] for line in (d/'gradient.dat').read_text().splitlines()],256)
    r={'components_hartree':parts,'total_hartree':math.fsum(parts.values()),'charge_e':charge,'rounding_bounds_hartree':[0.,0.],
       'cavity_global':values('CAVITY_GLOBAL',4),'gradient_units':'Hartree/Bohr; archived only',
       'artifacts':{p.name:pin(p) for p in d.iterdir() if p.is_file()}}
    write(d/'RESULT.json',r);return r


def comparison(rows,mex,pool_ok):
    mexdiff=difference(mex['Dy'],mex['La']);by={(r['source_id'],r['metal']):r for r in rows};reports=[]
    def exchange(label,la,dy):
        d=difference(dy['native_energies'],la['native_energies']);components={k:math.fsum([mexdiff[k],-d[k]]) for k in NAMES}
        alpb={k:((mex['Dy']['ALPB_parts_Eh'][k]-mex['La']['ALPB_parts_Eh'][k])-
            (dy['native_energies']['ALPB_parts_Eh'][k]-la['native_energies']['ALPB_parts_Eh'][k]))*HA for k in PARTS}
        return {'scope':label,'status':'available','LBT_Dy_minus_La_components_kcal_mol':d,'Mex_Dy_minus_La_components_kcal_mol':mexdiff,
            'R_components_kcal_mol':components,'R_kcal_mol':math.fsum(components.values()),'R_ALPB_subcomponents_kcal_mol':alpb,
            'LBT_selected_sources':{'La':la['source_id'],'Dy':dy['source_id']}}
    for sid in SOURCES:
        pair=[by[(sid,m)] for m in METALS]
        reports.append(exchange(sid,*pair) if all(r['status']=='complete' for r in pair) else {'scope':sid,'status':'unavailable'})
    pool={'status':'unavailable','identity_proved':pool_ok,'declared_sources':list(SOURCES),'requires_all_four_cells':True}
    if pool_ok and all(r['status']=='complete' for r in rows):
        selected={m:min((by[(s,m)] for s in SOURCES),key=lambda r:r['total_energy_kcal_mol']) for m in METALS}
        pool.update(exchange('common_two_source_pool',selected['La'],selected['Dy']))
        pool['exact_ties']={m:[s for s in SOURCES if by[(s,m)]['total_energy_kcal_mol']==selected[m]['total_energy_kcal_mol']] for m in METALS}
        pool['interpretation']='finite common source menu; not equilibrium populations'
    return {'source_specific':reports,'common_two_source_pool':pool,'sign':'positive R favors Mex-La plus LBT-Dy',
            'exchange':'Mex-La(+2) + LBT-Dy(-2) -> Mex-Dy(+2) + LBT-La(-2)'}


def run(cfg,out,result):
    cells,mex,scfg,build,source,pool_ok=inputs(cfg,out);result['preflight']=pin(out/'PREFLIGHT.json')
    result.update(reused_LBT_MACE_cells=4,reused_Mex_static_cells=2)
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',LD_LIBRARY_PATH=scfg['library_directory'])
    for key in ('LD_PRELOAD','PYTHONHOME','PYTHONPATH'):env.pop(key,None)
    result['native_environment']={k:env.get(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','LD_LIBRARY_PATH','LD_PRELOAD','PYTHONHOME','PYTHONPATH')}
    fp=out/'LBT_FIXED_SOLVENT.f90';fp.write_text(source);binary=out/'lbt_fixed_solvent'
    command=scfg['compiler_command']+[str(fp),'-o',str(binary)]+scfg['link_arguments']
    with (out/'compile.out').open('x') as so,(out/'compile.err').open('x') as se:proc=subprocess.run(command,cwd=out,env=env,stdout=so,stderr=se)
    result['build']={'command':command,'returncode':proc.returncode,'fortran':pin(fp),'build_inputs':scfg['build_inputs'],'runtime_libraries':build['linked_libraries']}
    require(proc.returncode==0,'Native wrapper build failed');bp=pin(binary);result['build']['binary']=bp
    for c,r in zip(cells,result['rows']):
        try:
            for p in build['linked_libraries']:verify(p)
            result['new_ALPB_calls_attempted']+=1;r['attempted']=True;d=out/(c['source_id']+'_'+c['metal'])
            sol=solvent(bp,c,d,env);v={'MACE_eV':c['MACE_eV'],'field_Eh':0.,'MM_kcal_mol':0.,'ALPB_Eh':sol['total_hartree'],'ALPB_parts_Eh':sol['components_hartree']}
            energies=native_components(v);r.update(status='complete',ALPB=pin(d/'RESULT.json'),MACE_result=cfg['mace_result'],MACE_task_id=c['task']['task_id'],
                native_energies=v,component_energies_kcal_mol=energies,total_energy_kcal_mol=math.fsum(energies.values()),partition_defined_zeros=['GFN2_field_difference','MM_ledger'])
        except Exception as exc:r.update(failure=repr(exc),traceback=traceback.format_exc())
    result['comparison']=comparison(result['rows'],mex,pool_ok)
    result['status']='complete' if all(r['status']=='complete' for r in result['rows']) else 'unavailable'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Slurm required');require(int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==1,'Declared one CPU required')
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);cfg=json.loads(a.config.read_text());start=time.monotonic()
    result={'status':'running','protocol_id':'lbt_whole_MACE_native_fixedFF_ALPB_Mex_balanced_exchange_v1',
        'implementation':pin(__file__),'plan':cfg['plan'],'config':pin(a.config),'previous_failure':cfg['previous_failure'],
        'recovery_note':cfg['recovery_note'],'slurm_job_id':os.environ['SLURM_JOB_ID'],
        'host':socket.gethostname(),'declared_cells':4,'new_ALPB_calls_attempted':0,'new_MACE_calls':0,'new_GFN2_calls':0,'new_MM_calls':0,'new_DFT_calls':0,
        'optimizer_steps':0,'reused_LBT_MACE_cells':0,'reused_Mex_static_cells':0,'affinity':None,
        'rows':[{'source_id':s,'metal':m,'status':'unavailable','attempted':False} for s in SOURCES for m in METALS]}
    try:run(cfg,out,result)
    except Exception as exc:result.update(status='execution_failure',failure=repr(exc),traceback=traceback.format_exc())
    result['elapsed_seconds']=time.monotonic()-start;result['application_allocated_CPU_seconds']=result['elapsed_seconds']
    write(out/'RESULT.json',result)
    lines=['# Matched fixed-FF solvent static comparison','',f"Status: {result['status']}. New ALPB attempts: {result['new_ALPB_calls_attempted']}/4.",'']
    if 'comparison' in result:
        lines+=['| LBT scope | Status | R, kcal/mol |','|---|---|---:|']
        for r in result['comparison']['source_specific']+[result['comparison']['common_two_source_pool']]:
            lines.append(f"| {r.get('scope','common_two_source_pool')} | {r['status']} | {r.get('R_kcal_mol','unavailable')} |")
        lines+=['','| LBT scope | MACE | Field difference | MM ledger | Fixed-FF ALPB |','|---|---:|---:|---:|---:|']
        for r in result['comparison']['source_specific']+[result['comparison']['common_two_source_pool']]:
            if r['status']=='available':lines.append('| '+r['scope']+' | '+' | '.join(str(r['R_components_kcal_mol'][k]) for k in NAMES)+' |')
    lines+=['','Positive R favors Mex–La + LBT–Dy. Exact component contributions and all four cells are in RESULT.json.',
        'Conditional fixed EF3 with three Nd spectators versus whole peptide; different Q partitions, source hydration, terminal and electronic-state limitations remain.',
        'The common pool is a finite two-source menu, not populations. No protein Kd or classifier claim.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n');return 0 if result['status']=='complete' else 1


if __name__=='__main__':raise SystemExit(main())
