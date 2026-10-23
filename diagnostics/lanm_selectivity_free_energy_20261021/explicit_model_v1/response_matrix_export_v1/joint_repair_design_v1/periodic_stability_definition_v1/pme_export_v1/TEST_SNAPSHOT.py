"""Finite artificial PME exporter/replay checks; no native invocation."""
import copy
import difflib
import hashlib
import io
import json
from pathlib import Path
import sys
import time
import numpy as np
import SNAPSHOT as S

BASE=Path(__file__).resolve().parent


def artifact(s):
    stats=S.replay(s);shape=s['shape'];p=s['p'];sites=s['sites'];M=np.prod(shape)
    lines=['PME_NORM_EXPORT_V1',f'GRID 1 {shape[0]} {shape[1]} {shape[2]} {p} {len(sites)} {len(sites)} {M}',
           'FFT FFTPACK',f'CHUNKS 1 1 1 {shape[0]} {shape[1]} {shape[2]} 0 1 2','EWALD_VOLUME .4 1000']
    for j,row in enumerate(s['trans'],1): lines.append('TRANS '+str(j)+' '+' '.join(map(repr,row.tolist())))
    lines.append('Q_STATS '+' '.join(str(stats[k]) for k in ['q_min','q_max','q_abs_max','q_zero','q_conjugate_difference_max']))
    lines.append('NORMS '+' '.join(str(stats[k]) for k in ['c1','cinf','cF','U_R_nominal','evenness_perturbation_nominal']))
    for site in sites:
        lines.append('SITE '+' '.join(map(str,[site['atom'],site['type'],*site['igrid'],site['alpha'],*site['xyz']])))
        for j in range(p):
            for b in range(3):lines.append('THETA '+' '.join(map(str,[site['atom'],b+1,j+1,*site['theta'][b,j]])))
    for z in range(shape[2]):
        for y in range(shape[1]):
            for x in range(shape[0]):lines.append(f'QFAC {x+1} {y+1} {z+1} {s["qfac"][x,y,z]!r}')
    lines.append(f'END {len(sites)} {3*p*len(sites)} {M} {len(sites)*p**3}')
    return '\n'.join(lines)+'\n'


def run():
    checks=[]
    def check(name,value):
        checks.append({'name':name,'passed':bool(value)})
        if not value:raise AssertionError(name)
    def close(name,a,b):check(name,np.allclose(a,b,rtol=1e-12,atol=1e-12))
    theta=np.array([[[.7,-1.],[.3,1.]],[[.6,-1.],[.4,1.]],[[.8,-1.],[.2,1.]]])
    sites=[dict(atom=i+1,type=i+5,alpha=[.4,.7][i],igrid=[-2+i,1,3-i],xyz=[0.,0.,0.],theta=theta.copy()) for i in range(2)]
    q=np.zeros((4,4,4))
    for index in np.ndindex(q.shape):
        z=sum(min(k,4-k)**2 for k in index)
        q[index]=1/(1+z) if z else 0.
    s=dict(shape=q.shape,p=2,sites=sites,qfac=q,trans=np.array([[.4,.02,0.],[0.,.5,.03],[.01,0.,.6]]))
    streamed=S.replay(s)
    B=np.zeros((64,6))
    for grid,i,b in S.entries(s):B[np.ravel_multi_index(grid,q.shape,order='F'),3*i:3*i+3]+=b
    close('norm1',streamed['c1'],np.linalg.norm(B,1))
    close('norm_inf',streamed['cinf'],np.linalg.norm(B,np.inf))
    close('frobenius_squared',streamed['cF'],np.square(B).sum())
    FB=np.stack([np.fft.fftn(B[:,j].reshape(q.shape,order='F')).ravel(order='F') for j in range(6)],axis=1)
    R=FB.conj().T@(q.ravel(order='F')[:,None]*FB)
    close('even_q_real_channel',R.imag,0.)
    check('nominal_upper_bounds_tiny_operator',np.linalg.eigvalsh(R.real)[-1]<=streamed['U_R_nominal']+1e-12)
    check('no_PSD_admission_from_file',streamed['PSD_admitted'] is False)
    doubled=copy.deepcopy(s)
    for site in doubled['sites']:site['alpha']*=4
    close('sqrt_alpha_normalization',S.replay(doubled)['U_R_nominal'],4*streamed['U_R_nominal'])
    wrapped=copy.deepcopy(s)
    for site in wrapped['sites']:site['igrid']=[x+8 for x in site['igrid']]
    close('wrapped_stencil_identity',S.replay(wrapped)['cF'],streamed['cF'])
    odd=copy.deepcopy(s);odd['qfac'][1,0,0]+=.2
    close('conjugate_difference',S.replay(odd)['q_conjugate_difference_max'],.2)
    check('asymmetry_retained',S.replay(odd)['evenness_perturbation_nominal']>0)
    neg=copy.deepcopy(s);neg['qfac'][0,0,0]=-.1
    check('negative_multiplier_reported',S.replay(neg)['q_negative_count']==1)
    text=artifact(s);parsed=S.read(io.StringIO(text));again=S.replay(parsed)
    for key in ['c1','cinf','cF','U_R_nominal','evenness_perturbation_nominal']:
        close('saved_replay_'+key,again[key],streamed[key])
    badcases={'missing_END':text[:text.rfind('END ')],
              'trailing':text+'EXTRA\n','duplicate_q':text.replace('QFAC 2 1 1','QFAC 1 1 1',1),
              'wrong_sequence':text.replace('GRID 1 ','GRID 2 ',1),
              'missing_theta':text.replace('THETA 1 1 1','THETA 1 1 2',1),
              'aliased_stencil':text.replace('GRID 1 4 4 4 2','GRID 1 4 4 4 5',1)}
    for name,bad in badcases.items():
        try:S.read(io.StringIO(bad))
        except ValueError:check('reject_'+name,True)
        else:check('reject_'+name,False)
    inputs=json.loads((BASE/'PATCH_INPUTS.json').read_text())
    original=Path(inputs['parent_source']['path']).read_text()
    check('parent_source_pin',hashlib.sha256(original.encode()).hexdigest()==inputs['parent_source']['sha256'])
    helper=(BASE/'PME_NORM_EXPORT.f').read_text()
    anchor=inputs['hook_anchor'];replacement=anchor.replace('      call fftback','      call pmenormexport\n      call fftback',1)
    check('unique_initialized_hook',original.count(anchor)==1 and original.index(anchor)>original.index('      qfac(1,1,1) = 0.0d0'))
    changed=original.replace(anchor,replacement)+helper
    patch=''.join(difflib.unified_diff(original.splitlines(True),changed.splitlines(True),fromfile='a/source/induce.f',tofile='b/source/induce.f'))
    check('isolated_patch_literal',patch==(BASE/'PME_NORM_EXPORT.patch').read_text())
    check('absent_return_before_arrays',helper.index('if (nrequest.eq.0) return')<helper.index('allocated(qfac)'))
    check('new_file_collision_guard',"status='new'" in helper and 'if (ios.ne.0) goto 900' in helper)
    check('no_qgrid_access','qgrid' not in helper.lower())
    check('first_state_only','if (sequence.ne.1) return' in helper)
    check('fixed_form_lines',all(len(line)<=72 for line in helper.splitlines() if line and line[0].lower()!='c'))
    check('stability_guards_preserved',original.count('call periodicresponsestability')==changed.count('call periodicresponsestability') and
          original.count('call periodicresponsebackendguard')==changed.count('call periodicresponsebackendguard'))
    return checks


if __name__=='__main__':
    out=Path(sys.argv[1]);wall,cpu=time.perf_counter(),time.process_time()
    if out.exists():raise ValueError('Immutable output exists')
    result={'status':'failed','native_calls':0,'builds':0,'physical_box_calls':0,
            'backend_factorizations':0,'implementation_sha256':hashlib.sha256(Path(S.__file__).read_bytes()).hexdigest(),
            'controls_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    try:result['checks']=run();result['check_count']=len(result['checks']);result['status']='passed'
    except Exception as e:result['error']=repr(e);raise
    finally:
        result['cost']={'wall_seconds':time.perf_counter()-wall,'process_CPU_seconds':time.process_time()-cpu,'numerical_threads':1}
        out.parent.mkdir(parents=True,exist_ok=True)
        with out.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
