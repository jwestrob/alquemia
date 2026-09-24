"""Derive four equal MPI workers from one actual, exclusive Slurm allocation."""
import os
import re
import subprocess

def layout(cpus,memory_mib,slots):
    if cpus<4 or cpus%4:raise ValueError('four equal workers require a CPU count divisible by four')
    ranks=cpus//4
    if ranks>slots:raise ValueError('an MPI worker would exceed declared Slurm task slots')
    # ORCA MaxCore is per MPI rank in MB, Slurm memory is MiB. Reserve 25% for
    # allocations outside MaxCore and round downward, never upward.
    maxcore=int(memory_mib*1.048576*.75/cpus/100)*100
    if maxcore<6000:raise ValueError('insufficient RAM for all CPUs at measured whole-source requirement')
    return ranks,memory_mib,maxcore

def main():
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
    memory=int(os.environ['SLURM_MEM_PER_NODE'])
    if memory==0:
        node=subprocess.check_output(['scontrol','show','node',field('NodeList'),'-o'],text=True)
        memory=int(re.search(r'(?:^|\s)RealMemory=(\d+)',node).group(1))
    print(*layout(cpus,memory,int(os.environ['SLURM_NTASKS'])))

if __name__=='__main__':main()
