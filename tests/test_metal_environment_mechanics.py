"""Parameter accounting on exact real 1H4I artifacts; no molecular evaluations."""
import json
import sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from metal_environment_mechanics import build, metal_parameters
BASE=ROOT/'workspaces/metal_environment_response_20260926'
METAL=Path('/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/dat/leap/parm/frcmod.ions234lm_126_tip3p')

@pytest.fixture(scope='module')
def ledger():
    return build(BASE/'preparation/scout_v3/INPUTS.json',BASE/'hybrid_preparation_v1/protein_v2/RESULT.json',BASE/'hybrid_preparation_v1/pqq_types_v1/PQQ_CROSS_LJ.json',METAL)

def test_real_partition(ledger):
    r,p,e,b,c=ledger
    assert r['counts']['physical_atoms']==9141
    assert r['counts']['real_QM']==51
    assert len(c)==3
    assert r['counts']['MM']==9090
    assert r['MM_charge_e']==pytest.approx(-6,abs=1e-12)
    ids={a['id']:a for a in p}
    for cap in c:
        omitted=ids[cap['omitted_source_id']]
        assert omitted['region']=='MM' and omitted['charge_e']==0
        assert omitted['epsilon_kJ_mol']>0
        assert cap['id'] not in ids

def test_exceptions_after_redistribution(ledger):
    _,p,ex,_,_=ledger
    changed=[]
    for e in ex:
        a,b=e['atoms']
        if e['region']=='MM':
            assert e['chargeprod_e2']==pytest.approx(p[a]['charge_e']*p[b]['charge_e']*e['coulomb_scale'],abs=1e-14)
            if abs(e['chargeprod_e2']-e['original_chargeprod_e2'])>1e-12:changed.append(e)
        else:assert e['chargeprod_e2'] is None
    assert changed

def test_bonded_and_no_qm_charge(ledger):
    r,p,_,b,_=ledger
    assert all(any(p[i]['region']=='MM' for i in t['atoms']) for t in b)
    assert all(a['charge_e'] is None for a in p if a['region']=='QM')
    assert r['qualification']['molecular_energy_force_calls']==0
    assert not r['qualification']['full_hybrid_force_checks']

def test_exact_metal_family(tmp_path):
    params=metal_parameters(METAL)
    assert params['Ca']['Rmin_half_A']==1.649
    assert params['La']['epsilon_kcal_mol']==0.15060822
    copy=tmp_path/'frcmod.ionslm_1264_opc'
    copy.write_text(METAL.read_text())
    with pytest.raises(ValueError,match='pure 12-6'):metal_parameters(copy)
