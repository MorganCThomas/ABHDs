"""ProtBert mean-pooled embeddings for every sequence in the UniProt FASTA.

Resumable: per-accession tensors in data/derived/protbert/*.pt are reused; the model is only
loaded if some are missing. Always rewrites the consolidated data/derived/protbert.npz
(accessions, X) that the notebooks read. Needs the `embed` extra (torch, transformers).

Caveat: sequences longer than 510 residues are truncated; the mean includes [CLS]/[SEP].
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C


def parse_fasta(path):
    seqs, acc, parts = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if acc:
                seqs[acc] = "".join(parts)
            acc, parts = line[1:].split("|")[1], []
        else:
            parts.append(line)
    if acc:
        seqs[acc] = "".join(parts)
    return seqs


def mean_pool(hidden, mask):
    m = mask.unsqueeze(-1).expand(hidden.size()).float()
    return (hidden * m).sum(1) / m.sum(1).clamp(min=1e-9)


def embed_missing(seqs, todo):
    from transformers import BertModel, BertTokenizer
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Embedding {len(todo)} sequences on {device}")
    tok = BertTokenizer.from_pretrained(C.PROTBERT_MODEL, do_lower_case=False)
    model = BertModel.from_pretrained(C.PROTBERT_MODEL).to(device).eval()
    for i in range(0, len(todo), C.PROTBERT_BATCH):
        batch = todo[i:i + C.PROTBERT_BATCH]
        txt = [" ".join(re.sub(r"[UZOB]", "X", seqs[a])) for a in batch]
        enc = tok(txt, return_tensors="pt", padding=True, truncation=True, max_length=C.PROTBERT_MAX_LEN)
        enc = {k: v.to(device) for k, v in enc.items()}
        with torch.no_grad():
            out = model(**enc)
        for a, e in zip(batch, mean_pool(out.last_hidden_state, enc["attention_mask"])):
            torch.save(e.cpu(), C.PROTBERT_DIR / f"{a}.pt")
        print(f"  {min(i + C.PROTBERT_BATCH, len(todo))}/{len(todo)}")


def main():
    C.PROTBERT_DIR.mkdir(parents=True, exist_ok=True)
    seqs = parse_fasta(C.SEQUENCES_FASTA)
    todo = [a for a in (e["accession"] for e in json.load(open(C.ENTRIES_JSON))) if not (C.PROTBERT_DIR / f"{a}.pt").exists()]
    if todo:
        embed_missing(seqs, todo)
    accs = [e["accession"] for e in json.load(open(C.ENTRIES_JSON))]   # entries order, not FASTA order
    X = np.stack([torch.load(C.PROTBERT_DIR / f"{a}.pt", weights_only=True).numpy() for a in accs])
    np.savez(C.PROTBERT_NPZ, accessions=np.array(accs), X=X)
    print(f"Wrote {C.PROTBERT_NPZ} {X.shape}")


if __name__ == "__main__":
    main()
