#!/usr/bin/env python3
"""Plot actual transfer-report rows, with no fitting or molecular calls; Slurm only."""
import argparse
import collections
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import traceback

PYTHON = {'path':'/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python3.11',
          'sha256':'53244dfc692ccbfd6a065f053dfa2026e2f3e4755beef84d600c8f7dbc9b5b27'}
MPL_METADATA = {'path':'/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/lib/python3.11/site-packages/matplotlib-3.10.9.dist-info/METADATA',
                'sha256':'a37f68bc5fec598f62d7641cdde0b10e137514400edd912cca88836177798559'}
REPORTER_SHA = '5c526e22a182ae24ef6c0f6cfefdd506963765b4eafbcefab8fc990df0004f99'
STYLE_PARENT_SHA = '589ed7e616aa3fadbdc4e6f021b3c2b526f3062533b6e1a4415d8a157a62683e'
POOLS = ('common3','common6')
FLAGS = ('lost_correct','new_wrong','lost_score_coverage')
STYLES = {'Ca':('#367F82','o'),'La':('#8A6195','^'),'2Ca_v1':('#BC8530','s'),None:('#687581','D')}

def require(value, message):
    if not value: raise ValueError(message)

def pin(path):
    p = Path(path).resolve()
    return {'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}

def check(record):
    require(pin(record['path'])==record,'Changed pinned artifact: '+record['path'])
    return Path(record['path'])

def read(path): return json.loads(Path(path).read_text())

def write(path, value):
    with Path(path).open('x') as f: json.dump(value,f,indent=2,allow_nan=False); f.write('\n')

def tsv(path):
    with Path(path).open(newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

def number(value):
    if value=='': return None
    x=float(value); require(math.isfinite(x),'Nonfinite exported score'); return x

def boolean(value):
    require(value in ('True','False'),'Unexpected exported boolean')
    return value=='True'

def data(directory, policy):
    names=('RESULT.json','EXECUTION.json','INPUTS.json','MEMBERS.tsv','TRANSITIONS.tsv')
    inputs={name:pin(directory/name) for name in names}
    result=read(directory/'RESULT.json'); execution=read(directory/'EXECUTION.json'); source=read(directory/'INPUTS.json')
    require(execution['status']=='complete' and execution['implementation']==source['implementation'] and source['implementation']['sha256']==REPORTER_SHA,'Expected completed pinned saved reporter')
    check(source['implementation'])
    require(result['status'] in ('complete','partial') and result['declared_contexts']==257 and result['declared_memberships']==260,'Reporter population differs')
    require(result['source_collection']==source['collection'] and result['source_execution']==source['execution'],'Reporter input joins differ')
    require(result['new_calibration'] is False and result['new_class_decisions'] is False and result['new_molecular_calls']==0 and result['consumed_correlated_development_evidence'],'Unexpected reporter scope')
    cal=result['recorded_calibrations']['calibrations']
    require(all(cal[p][policy]['status']=='available' for p in POOLS),'Fixed canonical calibration unavailable')
    bands=cal['common3'][policy]['bands']
    require(bands==cal['common6'][policy]['bands'] and bands['Ca_max']<bands['La_min'],'Canonical bands changed between displayed pools')
    members=[r for r in tsv(directory/'MEMBERS.tsv') if r['policy']==policy and r['pool'] in POOLS]
    by_record=collections.defaultdict(dict)
    for r in members:
        require(r['pool'] not in by_record[r['record_id']],'Duplicate exported membership/pool')
        by_record[r['record_id']][r['pool']]=r
    require(len(members)==520 and len(by_record)==260 and all(set(v)==set(POOLS) for v in by_record.values()),'All 260 memberships and both pools required')
    transitions=[r for r in tsv(directory/'TRANSITIONS.tsv') if r['policy']==policy and r['before_pool']=='common3' and r['after_pool']=='common6']
    tr={r['record_id']:r for r in transitions}
    require(len(transitions)==len(tr)==260 and set(tr)==set(by_record),'Full recorded transition denominator differs')
    contexts={}
    for rid, pair in by_record.items():
        before,after=(pair[p] for p in POOLS); transition=tr[rid]
        identity=('context_id','source_case_id','case_id','protein','biological_group','source_conditioning_metal','expected_class')
        require(all(before[k]==after[k]==transition[k] for k in identity),'Exported source/context identity changed')
        require(transition['before_R_kcal_mol']==before['R_kcal_mol'] and transition['after_R_kcal_mol']==after['R_kcal_mol'],'Recorded scores/transitions differ')
        row={k:before[k] for k in identity}; row['conditioning']=before['source_conditioning_metal'] or None
        row['scores']={p:number(pair[p]['R_kcal_mol']) for p in POOLS}
        row['outcomes']={p:pair[p]['outcome'] for p in POOLS}
        row['decisions']={p:pair[p]['decision'] for p in POOLS}
        row['pool_status']={p:pair[p]['pool_status'] for p in POOLS}
        row['selected']={p:{z:pair[p][z+'_candidate'] for z in ('Ca','La')} for p in POOLS}
        row['regression_flags']={k:boolean(transition[k]) for k in FLAGS}
        require(all((row['scores'][p] is None)==(row['outcomes'][p]=='unavailable') for p in POOLS),'Unavailable score/outcome mismatch')
        cid=row['context_id']
        if cid not in contexts: contexts[cid]={'scientific':row,'record_ids':[],'cohorts':set()}
        require(contexts[cid]['scientific']==row,'Aliased memberships differ; refusing selective deduplication')
        contexts[cid]['record_ids'].append(rid); contexts[cid]['cohorts'].update(json.loads(before['cohorts']))
    require(len(contexts)==257 and sum(len(v['record_ids']) for v in contexts.values())==260,'Expected 257 distinct contexts / 260 memberships')
    rows=[]
    for cid,v in contexts.items():
        row={**v['scientific'],'record_ids':sorted(v['record_ids']),'cohorts':sorted(v['cohorts'])}
        row['display_scores']={p:None if row['scores'][p] is None else row['scores'][p]-bands['Ca_max'] for p in POOLS}
        rows.append(row)
    for p in POOLS:
        require(sum(r['scores'][p] is not None for r in rows)==result['unique_context_coverage'][p][policy],'Unique coverage differs from saved reporter')
    return rows,bands,inputs,result

def source_order(row):
    hit=re.search(r'sample[-_](\d+)',row['source_case_id'])
    return (list(STYLES).index(row['conditioning']) if row['conditioning'] in STYLES else 99,
            int(hit[1]) if hit else -1,row['source_case_id'],row['context_id'])

def source_label(row):
    hit=re.search(r'sample[-_](\d+)',row['source_case_id']); condition={'2Ca_v1':'2Ca',None:'unspec.'}.get(row['conditioning'],row['conditioning'])
    label=f'{condition} s{hit[1]}' if hit else f'{condition} {row["source_case_id"]}'
    return label if len(label)<=25 else label[:22]+'…'

def render(rows,bands,out,policy):
    import matplotlib
    matplotlib.use('Agg')
    require(matplotlib.__version__=='3.10.9','Matplotlib version differs')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib import font_manager
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'svg.fonttype':'none','pdf.fonttype':42,
        'ps.fonttype':42,'axes.linewidth':.65,'axes.spines.top':False,'axes.spines.right':False,
        'savefig.facecolor':'white','svg.hashsalt':'uniform-amide-transfer-common3-common6-v1'})
    grouped=collections.defaultdict(list)
    for row in rows: grouped[row['biological_group']].append(row)
    groups=sorted(grouped); columns=4; nr=math.ceil(len(groups)/columns)
    h=max(2.6,.17*max(len(v) for v in grouped.values())+.9)
    fig=plt.figure(figsize=(18,nr*h+2),facecolor='white')
    grid=fig.add_gridspec(nr,columns,left=.055,right=.985,bottom=.075,top=.895,hspace=.52,wspace=.45)
    gap=bands['La_min']-bands['Ca_max']
    values=[v for row in rows for v in row['display_scores'].values() if v is not None]+[0,gap]
    span=max(values)-min(values); pad=max(2,.04*span); limits=(min(values)-pad,max(values)+pad)
    unknown={r['conditioning'] for r in rows}-set(STYLES)
    require(not unknown,'Unknown conditioning labels require an explicit legend: '+str(unknown))
    for index,group in enumerate(groups):
        inner=grid[index//columns,index%columns].subgridspec(1,2,width_ratios=(9,1),wspace=.05)
        ax=fig.add_subplot(inner[0,0]); missing=fig.add_subplot(inner[0,1],sharey=ax)
        members=sorted(grouped[group],key=source_order); n=len(members)
        ax.axvspan(limits[0],0,color='#EDF5F4',zorder=0); ax.axvspan(0,gap,color='#F3F3F3',zorder=0); ax.axvspan(gap,limits[1],color='#F8F0EB',zorder=0)
        for boundary in (0,gap): ax.axvline(boundary,color='#6A6A6A',linewidth=.7,linestyle='--',zorder=1)
        duplicates=collections.Counter(source_label(r) for r in members); seen=collections.Counter(); labels=[]
        for y,row in enumerate(members):
            label=source_label(row); seen[label]+=1
            if duplicates[label]>1: label+=f' [{seen[label]}]'
            row['figure_group']=group; row['figure_row_label']=label; labels.append(label)
            regression=any(row['regression_flags'].values()); color,marker=STYLES[row['conditioning']]
            vals=[row['display_scores'][p] for p in POOLS]
            if all(v is not None for v in vals):
                line,=ax.plot(vals,[y,y],color='#B64336' if regression else '#B6BCC2',linewidth=1.15 if regression else .7,zorder=2)
                line.set_gid(row['context_id']+'_transition')
            for k,(pool,value) in enumerate(zip(POOLS,vals)):
                if value is None:
                    artist=missing.scatter(k,y,marker='x',s=18,color='#B64336' if regression else '#687581',linewidths=1.0)
                else:
                    artist=ax.scatter(value,y,marker=marker,s=21 if k==0 else 13,facecolors='white' if k==0 else color,edgecolors=color,linewidths=.8,zorder=3+k)
                artist.set_gid(row['context_id']+'_'+pool)
            if regression: missing.text(1.8,y,'!',ha='center',va='center',fontsize=8,color='#B64336',weight='bold')
        classes='/'.join(sorted({r['expected_class'] for r in members}))
        names=sorted({r['protein'] for r in members}); title=group.upper()
        crystals=[name for name in names if '-pqq-' not in name]
        if crystals: title+=' / '+', '.join(crystals)
        ax.set_title(f'{title}  ·  known {classes}  ·  n={n}',loc='left',fontsize=9,weight='bold',pad=9)
        ax.set_xlim(*limits); ax.set_ylim(n-.5,-.5); ax.set_yticks(range(n),labels,fontsize=6.6)
        ax.tick_params(axis='y',length=0,pad=3); ax.tick_params(axis='x',labelsize=6.6,length=2)
        ax.locator_params(axis='x',nbins=4); ax.grid(axis='x',color='#D8DCE0',linewidth=.4,alpha=.6)
        for tick,row in zip(ax.get_yticklabels(),members):
            if any(row['regression_flags'].values()): tick.set_color('#B64336')
        ax.spines['left'].set_visible(False); ax.spines['bottom'].set_color('#999999')
        missing.set_xlim(-.55,2.2); missing.set_xticks([0,1],['3','6'],fontsize=6.5)
        missing.tick_params(axis='x',length=0,pad=2); missing.tick_params(axis='y',left=False,labelleft=False)
        missing.set_title('NA',fontsize=7,pad=9)
        for spine in missing.spines.values(): spine.set_visible(False)
    fig.text(.025,.972,'Uniform amide search: every current source context',fontsize=16,weight='bold')
    fig.text(.025,.950,f'Common3 → common6 · {policy} selection · 257 distinct contexts from 260 memberships · fixed canonical bands',fontsize=10,color='#444444')
    handles=[Line2D([],[],marker='o',color='#555555',linestyle='none',markerfacecolor='white',label='Common3'),Line2D([],[],marker='o',color='#555555',linestyle='none',label='Common6')]
    for condition,(color,marker) in STYLES.items():
        if any(r['conditioning']==condition for r in rows): handles.append(Line2D([],[],marker=marker,color=color,linestyle='none',label={'2Ca_v1':'2Ca conditioning',None:'Conditioning unspecified'}.get(condition,str(condition)+' conditioning')))
    handles.extend([Line2D([],[],marker='x',color='#687581',linestyle='none',label='Unavailable (NA gutter)'),Line2D([],[],color='#B64336',linewidth=1.2,label='Recorded regression / lost coverage')])
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.018,.938),ncol=4,frameon=False,fontsize=8.5,columnspacing=1.8)
    fig.text(.5,.050,'Recorded paired score R − fixed canonical Ca_max (kcal/mol); same display offset and x scale in every panel',ha='center',fontsize=10)
    fig.text(.025,.032,'Shading: Ca-supported | inconclusive | La-supported. Dashed lines are the saved canonical boundaries; no thresholds are refitted.',fontsize=9,color='#444444')
    fig.text(.025,.015,'Consumed development evidence. Structures, cohort overlaps and triples are correlated; no independent biological validation or production promotion.',fontsize=8.8,color='#444444')
    stem=out/('COMMON3_COMMON6_'+policy)
    fig.savefig(stem.with_suffix('.svg'),metadata={'Date':None,'Description':'All 257 distinct source/context pairs, including unavailable results; no averages or fitted thresholds.'})
    fig.savefig(stem.with_suffix('.pdf'),metadata={'Title':'Uniform amide transfer: individual source contexts','Creator':'Matplotlib; source and data pins in FIGURE_DATA.json','CreationDate':None,'ModDate':None})
    fig.savefig(stem.with_suffix('.png'),dpi=220); plt.close(fig)
    return {'matplotlib_version':matplotlib.__version__,'matplotlib_source':pin(matplotlib.__file__),
        'font':pin(font_manager.findfont('DejaVu Sans')),'groups':groups,'shared_display_xlim':limits,
        'outputs':[pin(stem.with_suffix(ext)) for ext in ('.svg','.pdf','.png')]}

def main():
    p=argparse.ArgumentParser(); p.add_argument('--report-dir',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
    p.add_argument('--script-sha',required=True); p.add_argument('--policy',choices=('operational','mathematical'),default='operational'); a=p.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'All reading/aggregation/rendering requires actual Slurm')
    require(pin(__file__)['sha256']==a.script_sha and pin(sys.executable)==PYTHON,'Implementation/interpreter pin differs'); check(MPL_METADATA)
    placement=subprocess.check_output(['scontrol','show','job','-o',os.environ['SLURM_JOB_ID']],text=True)
    fields=dict(v.split('=',1) for v in placement.split() if '=' in v)
    require(fields.get('JobState')=='RUNNING' and fields.get('NumNodes')=='1' and fields.get('NumCPUs')=='1','One running Slurm CPU/node required')
    out=a.output.resolve(); out.mkdir(parents=True,exist_ok=False); start=time.monotonic()
    receipt={'status':'failed','job_id':os.environ['SLURM_JOB_ID'],'placement':placement,'implementation':pin(__file__),'python':PYTHON,'matplotlib_metadata':MPL_METADATA,'new_molecular_calls':0}
    try:
        rows,bands,inputs,result=data(a.report_dir.resolve(),a.policy); rendering=render(rows,bands,out,a.policy)
        figure={'inputs':inputs,'reporter_sha256':REPORTER_SHA,'style_parent_sha256':STYLE_PARENT_SHA,'implementation':pin(__file__),
            'policy':a.policy,'pools':POOLS,'declared_memberships':260,'distinct_contexts':257,
            'deduplicated_memberships':3,'aliases':[{'context_id':r['context_id'],'record_ids':r['record_ids']} for r in rows if len(r['record_ids'])>1],
            'recorded_bands':bands,'display_offset_R_kcal_mol':bands['Ca_max'],'source_report_status':result['status'],
            'calibration_refit':False,'new_decisions':False,'group_averages':False,'rows':rows,'rendering':rendering}
        write(out/'FIGURE_DATA.json',figure)
        with (out/'FIGURE_DATA.tsv').open('x',newline='') as f:
            w=csv.writer(f,delimiter='\t'); w.writerow(['context_id','protein_group','source','conditioning','figure_label','record_ids','common3_R','common6_R','common3_outcome','common6_outcome',*FLAGS])
            for r in rows: w.writerow([r['context_id'],r['biological_group'],r['source_case_id'],r['conditioning'],r['figure_row_label'],json.dumps(r['record_ids']),r['scores']['common3'],r['scores']['common6'],r['outcomes']['common3'],r['outcomes']['common6'],*[r['regression_flags'][k] for k in FLAGS]])
        receipt.update(status='complete',figure_data=pin(out/'FIGURE_DATA.json'),outputs=rendering['outputs'],source_report_status=result['status'])
    except BaseException as exc:
        receipt.update(reason=str(exc),traceback=traceback.format_exc()); raise
    finally:
        receipt['elapsed_seconds']=time.monotonic()-start; write(out/'EXECUTION.json',receipt)

if __name__=='__main__': main()
