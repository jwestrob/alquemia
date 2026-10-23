"""Pure prospective key transformation; no files, candidate loading or native calls."""
import hashlib
import math

SETTINGS='''
# Operator-only derivative input: no energy/model/stability admission.
ewald
ewald-cutoff 9.0
vdw-cutoff 9.0
vdw-taper 8.1
ewald-alpha 0.4
pewald-alpha 0.4
pme-grid 64 64 64
pme-order 5
ppme-order 5
fft-package FFTW
polarization MUTUAL
polar-iter 100
PME-NORM-EXPORT 1
'''


def derive(composed_key):
    if not isinstance(composed_key,str) or not composed_key.endswith('\n'):
        raise ValueError('Exact newline-terminated composed key required')
    banned={'ewald','ewald-cutoff','ewald-boundary','ewald-alpha','pewald-alpha','dewald','pme-grid',
            'pme-order','ppme-order','fft-package','vdw-cutoff','vdw-taper','vdw-correction',
            'polarization','polar-iter','pme-norm-export','response-matrix-export','polar-donor-response',
            'external-field','polarizable','solvate','a-axis','b-axis','c-axis','alpha','beta','gamma'}
    kept=[];removed=[];seen=set()
    for index,line in enumerate(composed_key.splitlines(keepends=True)):
        words=line.split()
        if not words or words[0].startswith('#'):
            kept.append(line);continue
        key=words[0].lower()
        if key in banned:raise ValueError('Unexpected existing control '+key)
        if key!='polar-local-response':
            kept.append(line);continue
        if len(words)!=7:raise ValueError('Response record shape')
        metal,oxygen=map(int,words[1:3]);a,b,c,d=map(float,words[3:])
        if metal not in (900,901) or metal in seen or oxygen!=402:
            raise ValueError('Response namespace/duplicate')
        if not all(math.isfinite(v) for v in (a,b,c,d)) or max(abs(a),abs(b))>1 or abs(c)+abs(d)>math.log(4)/2+1e-15:
            raise ValueError('Original response domain')
        seen.add(metal);removed.append(dict(line_index=index,literal_line=line))
    if seen!={900,901}:raise ValueError('Both exact composed response records required')
    derivative=''.join(kept)+SETTINGS
    reconstructed=list(kept)
    for item in removed:reconstructed.insert(item['line_index'],item['literal_line'])
    if ''.join(reconstructed)!=composed_key:raise ValueError('Literal reconstruction failure')
    sha=lambda s:hashlib.sha256(s.encode()).hexdigest()
    return derivative,dict(status='operator_only_derivative_not_energy_model',original_key_sha256=sha(composed_key),
                           derivative_key_sha256=sha(derivative),removed=removed,appended_settings_sha256=sha(SETTINGS),
                           exact_reconstruction=True,native_calls=0,energy_admitted=False,full_stability_admitted=False)
