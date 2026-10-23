#!/usr/bin/env python3
"""Stage the declared 102-state finite panel; no native calls or model solves."""
import argparse,json,math,time
from pathlib import Path
from PREPARE import pin,read,verify,save,IDS

PAIR={'v0':.15,'vr':-.08,'va':.12,'hr':.003,'hf':-.002,'krr':.004,'kff':.003}
TRIPLE={'v0':.20,'va':-.13,'vf':.09,'h_neighbor':.0015,'k_cross':.001}
STEPS={'h':1e-4,'half':5e-5}
def dirs(g):
    n=len(g['xyz_A']);m=g['metal_atom'];lig=g['ligands'][0];o=lig['donor_atoms'][0];frame=lig['frame_atoms_by_donor'][str(o)][0]
    v={}
    for name,i,j in [('metal_x',m,0),('donor_y',o,1),('frame_z',frame,2)]:
        q=[[0.,0.,0.] for _ in range(n)];q[i][j]=1.;v[name]=q
    q=[math.sin(j+1) for j in range(3*n)];norm=math.sqrt(math.fsum(x*x for x in q));v['whole_sine']=[[q[3*i+j]/norm for j in range(3)] for i in range(n)]
    return v
def render(original,xyz):
    rows=original.decode().splitlines();out=[rows[0]]
    for s,coords in zip(rows[1:],xyz):
        t=s.split();out.append(f'{int(t[0]):6d} {t[1]:6s} '+ ' '.join(f'{x:.17g}' for x in coords)+' '+ ' '.join(t[5:]))
    assert len(out)==len(xyz)+1
    return ('\n'.join(out)+'\n').encode()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--centers',type=Path,required=True);ap.add_argument('--phase',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();tick=time.monotonic();cpu=time.process_time();c=read(pin(a.centers));phase=a.phase.read_text();assert '102' in phase and '1e-4' in phase and '5e-5' in phase
    assert tuple(c['entries'])==IDS;a.out.mkdir(parents=True,exist_ok=False);states={};directions={};scout=[]
    for cid,center in c['entries'].items():
        g=center['geometry'];directions[cid]=dirs(g);xyz0=verify(center['inputs']['xyz']).read_bytes();key=verify(center['inputs']['key']).read_bytes();prm=verify(center['inputs']['parameters']).read_bytes()
        for mode in ('zero','combined'):
            slots=[('center',None,None,0)]
            slots += [(f'{name}__{tag}__{sign}',name,tag,s) for name in directions[cid] for tag in STEPS for sign,s in [('minus',-1),('plus',1)]]
            for suffix,direction,tag,sign in slots:
                id=f'{cid}__{mode}__{suffix}';folder=a.out/'states'/id;folder.mkdir(parents=True)
                xyz=g['xyz_A'] if direction is None else [[x+sign*STEPS[tag]*q for x,q in zip(row,drow)] for row,drow in zip(g['xyz_A'],directions[cid][direction])]
                data=xyz0 if direction is None else render(xyz0,xyz)
                (folder/'system.xyz').write_bytes(data);(folder/'system.key').write_bytes(key);(folder/'amoeba09.prm').write_bytes(prm)
                parsed=[[float(x) for x in line.split()[2:5]] for line in data.decode().splitlines()[1:]];assert parsed==xyz
                kinds=sorted({l['kind'] for l in g['ligands']});tk=sorted({tuple(sorted((l['kind'],r['kind']))) for i,l in enumerate(g['ligands']) for r in g['ligands'][i+1:]})
                states[id]={'id':id,'case_id':cid,'coefficient_mode':mode,'direction':direction,'step_tag':tag,'step_A':None if tag is None else STEPS[tag],'sign':sign,'fixture':center['fixture'],'geometry':dict(g,xyz_A=xyz),'center_geometry':g,'parameters':{'pairs':{} if mode=='zero' else {k:PAIR for k in kinds},'triples':[] if mode=='zero' else [{'kinds':list(k),'coefficients':TRIPLE} for k in tk]},'inputs':{k:pin(folder/v) for k,v in [('xyz','system.xyz'),('key','system.key'),('parameters','amoeba09.prm')]},'source_center_inputs':center['inputs'],'native_calls_planned':1}
                if cid==IDS[0] and (suffix=='center' or mode=='combined' and direction=='metal_x'):scout.append(id)
    assert len(states)==102 and len(scout)==6
    manifest={'status':'prepared_declared_102_state_native_force_panel','implementation':pin(__file__),'phase':pin(a.phase),'centers':pin(a.centers),'shared_baseline':c['configuration']['shared_baseline'],'kernel':c['configuration']['kernel'],'states':states,'directions':directions,'scout_ids':scout,'remaining_ids':[id for id in states if id not in scout],'maximum_native_programs':102,'scout_native_programs':6,'remaining_native_programs':96,'fd_criteria':{'absolute_kcal_mol_A':5e-5,'relative_to_absolute_analytic':1e-6,'steps_A':STEPS,'expected_full_comparisons':48},'runtime_build':None,'runtime_release':None,'new_native_calls':0,'operator_solves':0,'cost':{'wall_seconds':time.monotonic()-tick,'CPU_seconds':time.process_time()-cpu,'threads':1}}
    save(a.out/'MANIFEST.json',manifest);print(json.dumps({'status':manifest['status'],'states':len(states),'scout':len(scout),'cost':manifest['cost']}))
if __name__=='__main__':main()
