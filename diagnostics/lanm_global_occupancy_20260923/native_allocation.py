"""Derive four equal MPI workers from one actual, exclusive Slurm allocation."""
import os
import re
import subprocess
import argparse

def layout(cpus,memory_mib,slots,workers=4):
    if workers<1 or cpus<workers or cpus%workers:raise ValueError('equal workers require a divisible CPU count')
    ranks=cpus//workers
    if ranks>slots:raise ValueError('an MPI worker would exceed declared Slurm task slots')
    # ORCA MaxCore is per MPI rank in MB, Slurm memory is MiB. Reserve 25% for
    # allocations outside MaxCore and round downward, never upward.
    maxcore=int(memory_mib*1.048576*.75/cpus/100)*100
    if maxcore<6000:raise ValueError('insufficient RAM for all CPUs at measured whole-source requirement')
    return ranks,memory_mib,maxcore

def whole_node_memory(job_record,node_record):
    # --mem=0 may omit SLURM_MEM_PER_NODE entirely. Read the authoritative
    # scheduler request and registered node memory instead of that optional env.
    requested=re.search(r'(?:^|\s)MinMemoryNode=(\S+)',job_record)
    if requested is None or requested.group(1)!='0':
        raise ValueError('whole-node memory was not granted by the scheduler')
    registered=re.search(r'(?:^|\s)RealMemory=(\d+)',node_record)
    if registered is None:raise ValueError('missing registered node memory')
    return int(registered.group(1))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    job=os.environ['SLURM_JOB_ID']
    raw=subprocess.check_output(['scontrol','show','job',job,'-o'],text=True)
    def field(key):
        match=re.search(r'(?:^|\s)'+key+r'=(\S+)',raw)
        if match is None:raise ValueError('missing allocation field '+key)
        return match.group(1)
    if field('NumNodes')!='1' or field('OverSubscribe')!='NO':
        raise ValueError('requires one exclusive node')
    cpus=int(os.environ['SLURM_CPUS_ON_NODE'])
    if cpus!=int(field('NumCPUs')) or len(os.sched_getaffinity(0))<cpus:
        raise ValueError('scheduler CPU allocation and process affinity differ')
    node=subprocess.check_output(['scontrol','show','node',field('NodeList'),'-o'],text=True)
    memory=whole_node_memory(raw,node)
    print(*layout(cpus,memory,int(os.environ['SLURM_NTASKS']),args.workers))

if __name__=='__main__':main()
