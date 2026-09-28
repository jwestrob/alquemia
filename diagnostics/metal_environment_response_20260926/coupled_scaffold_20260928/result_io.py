"""Research-only result writer: native JSON conversion before exclusive creation."""
import json
from pathlib import Path

def _native(value):
    if type(value).__module__.split('.')[0] == 'numpy':
        import numpy as np
        if isinstance(value,np.ndarray):return value.tolist()
        if isinstance(value,np.generic):
            result=value.item()
            if not isinstance(result,np.generic):return result
    raise TypeError(f'Object of type {type(value).__name__} is not JSON serializable')

def write_new(path,value):
    text=json.dumps(value,indent=2,sort_keys=True,allow_nan=False,default=_native)+'\n'
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x') as f:f.write(text)

def write_configuration(path,*,energies,forces,**metadata):
    """Call immediately after a configuration, before computing aggregate gates.

    Forces must be the actual evaluated arrays; no missing-value substitution.
    This function performs no molecular calculation or scientific admission.
    """
    if energies is None or forces is None:raise ValueError('actual energy and force receipt required')
    write_new(path,dict(energies=energies,forces=forces,**metadata))
