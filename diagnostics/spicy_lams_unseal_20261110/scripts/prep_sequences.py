"""Sequence-side preparation for all 616 workbook orthologs (sequences are metadata, never sealed):
pairwise global identities and ESM-2 650M mean-pooled embeddings (cached weights, CPU, offline)."""
import os, sys, time
os.environ["HF_HUB_OFFLINE"] = "1"; os.environ["TRANSFORMERS_OFFLINE"] = "1"
import numpy as np, pandas as pd
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from common import *  # noqa
from Bio import Align
from Bio.Align import substitution_matrices
o = pd.read_csv(ROOT / "data/ORTHOLOGS.tsv", sep="\t")
acc = o.accession.values; seqs = o.mature_sequence.values; n = len(acc)
pd.DataFrame(dict(accession=acc, sequence=seqs, family=o.family, genus=o.genus, reserved=o.reserved_panel)).to_csv(ROOT / "data/all616_sequences.tsv", sep="\t", index=False)
t0 = time.time()
al = Align.PairwiseAligner(); al.mode = "global"; al.substitution_matrix = substitution_matrices.load("BLOSUM62")
al.open_gap_score = -10; al.extend_gap_score = -0.5
I = np.eye(n)
for i in range(n):
    for j in range(i + 1, n):
        a = al.align(seqs[i], seqs[j])[0]
        s1, s2 = a[0], a[1]
        m = sum(1 for x, y in zip(s1, s2) if x == y and x != "-")
        I[i, j] = I[j, i] = m / len(s1)
    if i % 100 == 0: print(f"identity row {i}/{n} ({time.time()-t0:.0f} s)", flush=True)
np.save(ROOT / "data/identity_616.npy", I)
print(f"identity matrix done in {time.time()-t0:.0f} s; median off-diagonal {np.median(I[np.triu_indices(n,1)]):.3f}")
import torch
from transformers import AutoTokenizer, EsmModel
torch.set_num_threads(6)
tok = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D"); model = EsmModel.from_pretrained("facebook/esm2_t33_650M_UR50D").eval()
E = np.zeros((n, 1280), dtype=np.float32); t0 = time.time()
with torch.no_grad():
    for s in range(0, n, 16):
        b = list(seqs[s:s + 16]); x = tok(b, return_tensors="pt", padding=True)
        h = model(**x).last_hidden_state
        for k, sq in enumerate(b):
            E[s + k] = h[k, 1:len(sq) + 1].mean(0).numpy()
        if s % 160 == 0: print(f"embedding {s}/{n} ({time.time()-t0:.0f} s)", flush=True)
np.save(ROOT / "data/esm2_650M_mean_616.npy", E); print("embeddings done", E.shape)
