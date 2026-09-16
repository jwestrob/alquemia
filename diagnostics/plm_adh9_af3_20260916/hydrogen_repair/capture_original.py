#!/usr/bin/env python3
"""Observe the unchanged frozen hydrogen minimizer; preserve its failure."""
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
import openmm as mm
from openmm import unit
import repair_hydrogens as recovery

ROOT=recovery.ROOT
OUT=ROOT/'workspaces/plm_adh9_af3_20260916/hydrogen_repair/original_capture'
CAL=ROOT/'diagnostics/pqq_pmdh_fixed_core_calibration_20260914'

def main():
    OUT.mkdir(exist_ok=False,parents=True)
    spec=importlib.util.spec_from_file_location('capture_frozen_protonator',CAL/'protonate_standard_only.py')
    helper=importlib.util.module_from_spec(spec);sys.modules[spec.name]=helper;spec.loader.exec_module(helper)
    original=mm.LocalEnergyMinimizer.minimize;calls=[]
    def observe(context,tolerance=10,maxIterations=0,reporter=None):
        system=context.getSystem();n=len(calls)
        hydrogen=np.asarray([system.getParticleMass(i).value_in_unit(unit.dalton)>0 for i in range(system.getNumParticles())])
        xml=OUT/f'system_{n}.xml';xml.write_text(mm.XmlSerializer.serialize(system))
        before,bm=recovery.state(context,hydrogen);np.save(OUT/f'initial_positions_{n}_nm.npy',before)
        trace=recovery.Reporter();original(context,tolerance,maxIterations,trace)
        after,am=recovery.state(context,hydrogen);np.save(OUT/f'final_positions_{n}_nm.npy',after)
        calls.append({'system':recovery.record(xml),'tolerance':tolerance,'iterations':maxIterations,
                      'initial':bm,'final':am,'positions_unchanged':bool(np.array_equal(before,after)),
                      'reporter_iterations':len(trace.entries),'trace':trace.entries})
    mm.LocalEnergyMinimizer.minimize=observe
    try:
        helper.protonate_standard_only(recovery.OLD/(recovery.STEM+'_normalized.pdb'),OUT/'protonated.pdb',pins_path=CAL/'implementation_pins.json')
    finally:mm.LocalEnergyMinimizer.minimize=original
    recovery.write(OUT/'observation.json',{'script':recovery.record(__file__),'frozen_wrapper':recovery.record(CAL/'protonate_standard_only.py'),'calls':calls})
    print(json.dumps([{k:v for k,v in x.items() if k!='trace'} for x in calls],indent=2))

if __name__=='__main__':main()
