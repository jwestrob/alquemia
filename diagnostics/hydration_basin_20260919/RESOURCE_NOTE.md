# Execution resource update, before the first computation

The requested standard RTX A5000 node had all GPUs occupied. The H200 node
node-224-2t-8gpu-1 had available GPUs. On2026-09-19, the agent changed only its
own pending job1202050's node constraint and memory request via scontrol:
ReqNodeList=node-224-2t-8gpu-1, MinMemoryNode=200000MiB. OneGPU/16CPUs retained.
Priority, other jobs, scientific manifest, inputs and implementation are unchanged.
The200000MiB request is the previously accepted195GiB H200 scheduler allocation;
it supersedes the initial64474MiB standard-host request in AGREEMENT.md.
No running job was canceled or interrupted; no duplicate computation submitted.
