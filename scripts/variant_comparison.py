#!/usr/bin/env python3
import os, re, glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT = "/Users/vishnu/Desktop/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/CxC_TxT" 
VAR_DIRS = ["cxc_results","cxs_results","sxc_results","sxt_results","txs_results","txt_results","txc_results","cxt_results"]
OUTDIR = os.path.join(ROOT, "merged")
os.makedirs(OUTDIR, exist_ok=True)

def read_with_variant(path, variant):
    df = pd.read_csv(path)
    df.insert(0, "variant", variant)
    return df

tops, alls, rec_cons, pep_cons = [], [], [], []
for vd in VAR_DIRS:
    vname = vd.replace("_results","")
    base = os.path.join(ROOT, vd)
    tops.append(read_with_variant(os.path.join(base, "AF3_summary_top.csv"), vname))
    alls.append(read_with_variant(os.path.join(base, "AF3_summary_all.csv"), vname))
    rec = read_with_variant(os.path.join(base, "consensus_interface_receptor_af3.csv"), vname)
    pep = read_with_variant(os.path.join(base, "consensus_interface_peptide_af3.csv"), vname)
    rec_cons.append(rec); pep_cons.append(pep)

top_all = pd.concat(tops, ignore_index=True)
all_all = pd.concat(alls, ignore_index=True)
rec_all = pd.concat(rec_cons, ignore_index=True)
pep_all = pd.concat(pep_cons, ignore_index=True)

# Save merged tables
top_path = os.path.join(OUTDIR, "AF3_summary_top_merged.csv")
all_path = os.path.join(OUTDIR, "AF3_summary_all_merged.csv")
rec_path = os.path.join(OUTDIR, "consensus_receptor_by_variant.csv")
pep_path = os.path.join(OUTDIR, "consensus_peptide_by_variant.csv")
top_all.to_csv(top_path, index=False)
all_all.to_csv(all_path, index=False)
rec_all.to_csv(rec_path, index=False)
pep_all.to_csv(pep_path, index=False)

print("Wrote:")
print(" ", top_path)
print(" ", all_path)
print(" ", rec_path)
print(" ", pep_path)

# Leaderboard (mean across the best model per design in each variant folder)
score_cols = [c for c in ["rankingScore","ipTM","iptm","interface_pLDDT_mean","n_contacts"] if c in top_all.columns]
leader = (top_all.groupby("variant", as_index=False)[score_cols]
          .mean(numeric_only=True)
          .sort_values(["interface_pLDDT_mean",
                        "ipTM" if "ipTM" in score_cols else ("iptm" if "iptm" in score_cols else score_cols[0]),
                        "n_contacts" if "n_contacts" in score_cols else score_cols[-1]],
                       ascending=[False, False, False]))
leader_path = os.path.join(OUTDIR, "AF3_variant_leaderboard.csv")
leader.to_csv(leader_path, index=False)
print(" ", leader_path)

# Also collect per-variant receptor residue hits (top N per variant)
rec_hits = (rec_all.sort_values(["variant","frequency"], ascending=[True, False])
                 .groupby("variant", as_index=False)
                 .head(20))
rec_hits.to_csv(os.path.join(OUTDIR, "receptor_top20_residues_per_variant.csv"), index=False)

# --- Quick plots ---
# 1) interface_pLDDT_mean bar
plt.figure(figsize=(7,4))
x = leader["variant"]
y = leader["interface_pLDDT_mean"]
plt.bar(x, y)
plt.ylabel("Interface pLDDT (mean)")
plt.title("AF3: Interface confidence by variant (best-model means)")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "bar_interface_pLDDT_mean.png"), dpi=200)

# 2) n_contacts bar (if present)
if "n_contacts" in leader.columns:
    plt.figure(figsize=(7,4))
    plt.bar(leader["variant"], leader["n_contacts"])
    plt.ylabel("# contacts")
    plt.title("AF3: Contact count by variant (best-model means)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "bar_n_contacts.png"), dpi=200)

print("Plots saved in", OUTDIR)