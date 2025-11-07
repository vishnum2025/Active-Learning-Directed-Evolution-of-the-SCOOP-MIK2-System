#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Active learning over 13-mer SCOOP peptides (SxS motif enforced)
- Trains a small DNN ENSEMBLE on your labeled data (sequence, fitness)
- Scores the whole unlabeled domain
- Selects a diverse batch with UCB or Thompson sampling
Outputs:
  results_dir/round1_proposals.csv
  results_dir/round1_logo.png
"""

import os, argparse, time, random
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader
import matplotlib.pyplot as plt

try:
    import logomaker as lm
    HAVE_LOGO = True
except Exception:
    HAVE_LOGO = False

AA20 = list("ACDEFGHIKLMNPQRSTVWY")
AA2I = {a:i for i,a in enumerate(AA20)}

def enforce_and_filter(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["sequence"] = (
        df["sequence"].astype(str).str.upper()
        .str.replace(r"[^A-Z]", "", regex=True).str.strip()
    )
    # 13 aa and contains SxS motif anywhere
    df = df[(df["sequence"].str.len()==13) & (df["sequence"].str.contains("S.S", regex=True))]
    return df.drop_duplicates(subset=["sequence"]).reset_index(drop=True)

def onehot(seqs: list[str]) -> torch.Tensor:
    L = len(seqs[0])
    X = torch.zeros((len(seqs), L*20), dtype=torch.float32)
    for r,s in enumerate(seqs):
        for j,a in enumerate(s):
            i = AA2I.get(a, None)
            if i is not None:
                X[r, j*20 + i] = 1.0
    return X

class MLP(nn.Module):
    def __init__(self, d_in, hidden=(256,128)):
        super().__init__()
        layers = []
        h_prev = d_in
        for h in hidden:
            layers += [nn.Linear(h_prev, h), nn.LeakyReLU(0.1)]
            h_prev = h
        layers += [nn.Linear(h_prev, 1)]
        self.net = nn.Sequential(*layers)
    def forward(self, x): return self.net(x)

def train_one(model, X, y, epochs=400, lr=1e-3, wd=1e-4, bs=64, seed=0):
    torch.manual_seed(seed); np.random.seed(seed); random.seed(seed)
    ds = TensorDataset(X, y.view(-1,1))
    dl = DataLoader(ds, batch_size=bs, shuffle=True, drop_last=False)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=wd)
    loss_fn = nn.MSELoss()
    model.train()
    for _ in range(epochs):
        for xb, yb in dl:
            opt.zero_grad()
            pred = model(xb)
            loss = loss_fn(pred, yb)
            loss.backward()
            opt.step()
    model.eval()

def ensemble_predict(models, X):
    with torch.no_grad():
        preds = [m(X).squeeze(-1) for m in models]
        P = torch.stack(preds, dim=0)  # [M, N]
        mean = P.mean(dim=0)
        std  = P.std(dim=0, unbiased=False)
    return mean, std

def select_batch(mean, std, k=96, mode="ucb", beta=2.0, existing_idx=set(), min_hamm=3, seqs=None):
    """Greedy selection with novelty filter (Hamming distance) to avoid near-duplicates."""
    n = mean.numel()
    if mode.lower()=="ucb":
        score = mean + beta*std
    elif mode.lower()=="ts":
        score = torch.normal(mean, std.clamp_min(1e-4))
    else:
        score = mean
    order = torch.argsort(score, descending=True).tolist()
    chosen = []
    def hamm(a,b):
        return sum(c1!=c2 for c1,c2 in zip(a,b))
    for idx in order:
        if idx in existing_idx: 
            continue
        if seqs and chosen:
            ok = True
            for j in chosen:
                if hamm(seqs[idx], seqs[j]) < min_hamm:
                    ok = False; break
            if not ok: 
                continue
        chosen.append(idx)
        if len(chosen) >= k: break
    return chosen

def make_logo(seqs, out_png, title):
    if not HAVE_LOGO or not seqs: return
    L = len(seqs[0])
    counts = {aa:[0]*L for aa in AA20}
    for s in seqs:
        for i,a in enumerate(s):
            if a in counts: counts[a][i]+=1
    df = pd.DataFrame(counts)
    probs = df.div(df.sum(axis=1), axis=0).fillna(0.0)
    plt.figure(figsize=(0.6*L,3))
    lm.Logo(probs)
    plt.title(title); plt.tight_layout(); plt.savefig(out_png, dpi=220); plt.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="/Users/vishnu/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/alde_athaliana")
    ap.add_argument("--subset_csv", default="alde_training_round0.csv")
    ap.add_argument("--domain_csv", default="alde_domain.csv")
    ap.add_argument("--outdir", default="results_ts/scoop_round1_fixed")
    ap.add_argument("--batch", type=int, default=96)
    ap.add_argument("--mode", default="ts", choices=["ucb","ts"])
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    torch.manual_seed(args.seed); np.random.seed(args.seed); random.seed(args.seed)

    OUT = os.path.join(args.base, args.outdir); os.makedirs(OUT, exist_ok=True)

    # --- labeled ---
    dfL = pd.read_csv(os.path.join(args.base, args.subset_csv))
    if not {"sequence","fitness"}.issubset(dfL.columns):
        raise SystemExit("Training CSV must have columns: sequence, fitness")
    dfL = enforce_and_filter(dfL)
    if len(dfL) < 3:
        raise SystemExit(f"Need ≥3 labeled peptides, found {len(dfL)} after filtering.")

    y = dfL["fitness"].astype(float).to_numpy()
    y = (y - y.min())/(y.max() - y.min() + 1e-12)
    X = onehot(dfL["sequence"].tolist())
    y = torch.tensor(y, dtype=torch.float32)

    # --- domain ---
    dfD = pd.read_csv(os.path.join(args.base, args.domain_csv))
    if "sequence" not in dfD.columns:
        raise SystemExit("Domain CSV must contain 'sequence'")
    dfD = enforce_and_filter(dfD)
    # remove labeled from domain
    labeled_set = set(dfL["sequence"].tolist())
    dfD = dfD[~dfD["sequence"].isin(labeled_set)].reset_index(drop=True)

    domain_seqs = dfD["sequence"].tolist()
    XD = onehot(domain_seqs)

    # --- ensemble training ---
    d_in = X.shape[1]
    M = 8
    models = []
    for m in range(M):
        net = MLP(d_in, hidden=(256,128))
        train_one(net, X, y, epochs=500, lr=1e-3, wd=1e-4, bs=64, seed=args.seed+m)
        models.append(net)

    # --- score domain ---
    mean, std = ensemble_predict(models, XD)

    # --- select diverse top-K ---
    chosen_idx = select_batch(mean, std,
                              k=args.batch,
                              mode=args.mode,
                              beta=2.0,
                              existing_idx=set(),
                              min_hamm=3,
                              seqs=domain_seqs)

    proposals = pd.DataFrame({
        "domain_index": chosen_idx,
        "sequence": [domain_seqs[i] for i in chosen_idx],
        "pred_mean": mean[chosen_idx].numpy(),
        "pred_std":  std[chosen_idx].numpy()
    })
    out_csv = os.path.join(OUT, "round1_proposals.csv")
    proposals.to_csv(out_csv, index=False)

    if len(proposals) > 0:
        make_logo(proposals["sequence"].tolist(),
                  os.path.join(OUT, "round1_logo.png"),
                  title="AL Round 1 Proposals")

    print(f"[OK] Labeled: {len(dfL)} | Domain: {len(dfD)} | Proposed: {len(proposals)}")
    print(f"[OK] Saved -> {out_csv}")
    print("Top 5:\n", proposals.head(5).to_string(index=False))

if __name__ == "__main__":
    main()