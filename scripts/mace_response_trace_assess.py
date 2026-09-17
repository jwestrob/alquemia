"""Summarize saved charge-update traces; no model evaluations."""
import argparse
from pathlib import Path
import time
import numpy as np
from affordable_common import InvalidArtifact, read_json, record, verify, write_new


def assess(medium, large, output):
    start=time.monotonic()
    result={'status':'complete','implementation':record(__file__),'models':{},'new_model_calls':0}
    for label,path in (('medium',medium),('large',large)):
        c=read_json(path)
        if c['status']!='complete' or not all(x['pass'] for x in c['trace_equivalence']):
            raise InvalidArtifact('trace replay did not reproduce its reference')
        model={'collection':record(path),'pairs':{},
               'max_energy_replay_error_eV':max(abs(x['energy_delta_eV']) for x in c['trace_equivalence']),
               'max_force_replay_error_eV_A':max(x['force_delta_eV_A'] for x in c['trace_equivalence']),
               'max_density_replay_error':max(x['density_delta'] for x in c['trace_equivalence'])}
        for pair in ('qm33','qm36','full'):
            data={}
            for metal in ('La','Ca'):
                name=f'1h4i_{pair}_{metal}'+('_primary' if pair=='full' else '')
                trace=read_json(verify(c['rows'][name]['charge_trace']))
                with np.load(verify(trace['arrays'])) as arrays:
                    data[metal]=(trace,{k:arrays[k].copy() for k in arrays.files})
            stages=[]
            for stage in range(data['La'][0]['recursion_steps']+1):
                key=f'stage_{stage}_restored_density_internal'
                charge={metal:data[metal][1][key][:,:,0].sum(axis=1) for metal in ('La','Ca')}
                delta=charge['Ca']-charge['La']
                correction={metal:data[metal][1][f'stage_{stage}_restoration_correction_e'].sum(axis=1) for metal in ('La','Ca')}
                change=correction['Ca']-correction['La']
                row={'stage':stage,'Ca_minus_La_charge_sum_e':float(delta.sum()),
                     'Ca_minus_La_charge_L1_e':float(np.abs(delta).sum()),
                     'Ca_minus_La_restoration_correction_L1_e':float(np.abs(change).sum()),
                     'Ca_minus_La_max_atomic_charge_change_e':float(np.abs(delta).max()),
                     'minimum_weight_cancellation_ratio':min(x for metal in ('La','Ca') for x in data[metal][0]['stages'][stage]['weight_sum_cancellation_ratio'])}
                if stage:
                    row['field_input_range']={}
                    for metal in ('La','Ca'):
                        field=data[metal][1][f'update_{stage-1}_potential_features']
                        row['field_input_range'][metal]={'min':float(field.min()),'max':float(field.max()),'RMS':float(np.sqrt(np.square(field).mean()))}
                stages.append(row)
            model['pairs'][pair]={'atoms':len(delta),'stages':stages}
        result['models'][label]=model
    result['wall_seconds']=time.monotonic()-start
    write_new(output,result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('medium','large','output'):p.add_argument('--'+k,required=True)
    a=p.parse_args();assess(a.medium,a.large,a.output)
