#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, json
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, ttest_ind, spearmanr



BASE = "/Users/vishnu/v_files/GT/sem1/Research_sem1"
AF3_DIR = f"{BASE}/scoop_mik2_af_pipeline/scoop_cluster_representatives/af3_cluster_outputs"
OUT_DIR = f"{BASE}/scoop_mik2_af_pipeline/scoop_cluster_representatives/af3_analysis"
os.makedirs(OUT_DIR, exist_ok=True)

FOLDERS = os.listdir(AF3_DIR)

records = []

for folder in FOLDERS:
    lvl2 = os.path.join(AF3_DIR, folder, folder)
    if not os.path.isdir(lvl2):
        continue

    conf = os.path.join(lvl2, f"{folder}_summary_confidences.json")
    rank = os.path.join(lvl2, f"{folder}_ranking_scores.csv")

    if not os.path.exists(conf) and not os.path.exists(rank):
        continue

    iptm = ptm = ranking = np.nan
    if os.path.exists(conf):
        try:
            with open(conf) as f:
                data = json.load(f)
            iptm = data.get("iptm", np.nan)
            ptm = data.get("ptm", np.nan)
            ranking = data.get("ranking_score", np.nan)
        except Exception as e:
            print(f"Failed reading {conf}: {e}")

    elif os.path.exists(rank):
        try:
            df_temp = pd.read_csv(rank)
            # attempt to find columns like iptm, ptm, ranking_score
            for c in df_temp.columns:
                lc = c.lower()
                if "iptm" in lc and pd.notna(df_temp[c].iloc[0]):
                    iptm = float(df_temp[c].iloc[0])
                elif "ptm" in lc and pd.notna(df_temp[c].iloc[0]):
                    ptm = float(df_temp[c].iloc[0])
                elif "ranking" in lc:
                    ranking = float(df_temp[c].iloc[0])
        except Exception as e:
            print(f"Failed reading {rank}: {e}")

    score_sum = np.nan if np.isnan(iptm) or np.isnan(ptm) else iptm + ptm

    parts = folder.split("__")
    if len(parts) != 2:
        continue
    receptor_species = parts[0].split("_")[0]
    peptide_species = parts[1].split("_")[0]

    records.append({
        "complex": folder,
        "receptor_species": receptor_species,
        "peptide_species": peptide_species,
        "iptm": iptm,
        "ptm": ptm,
        "iptm+ptm": score_sum,
        "ranking_score": ranking
    })

# --- Save ---
if not records:
    raise RuntimeError("No valid AF3 scores found!")

df = pd.DataFrame(records)
df["pair_type"] = np.where(df["receptor_species"] == df["peptide_species"],
                           "within_species", "cross_species")
out_csv = f"{OUT_DIR}/af3_cluster_scores.csv"
df.to_csv(out_csv, index=False)
print(f"Saved → {out_csv} ({len(df)} complexes)")

print("\nPreview:")
print(df.head())

# Optional summary stats
print("\nMean IPTM+PTM by pair type:")
print(df.groupby("pair_type")["iptm+ptm"].mean())


# ======================================================
#  Statistical comparison and visualization
# ======================================================

df_valid = df.dropna(subset=["iptm"])

within = df_valid.loc[df_valid["pair_type"] == "within_species", "iptm"]
cross  = df_valid.loc[df_valid["pair_type"] == "cross_species", "iptm"]

from scipy.stats import mannwhitneyu, ttest_ind
U, p_mwu = mannwhitneyu(within, cross, alternative="greater")
tval, p_t = ttest_ind(within, cross, equal_var=False, alternative="greater")

print("\n--- IPTM-only Hypothesis Testing ---")
print(f"Within n = {len(within)}, mean = {within.mean():.3f}")
print(f"Cross  n = {len(cross)}, mean = {cross.mean():.3f}")
print(f"Mann–Whitney U p = {p_mwu:.3e}")
print(f"Welch t-test p   = {p_t:.3e}")

# --- Hierarchical clustered heatmap ---
import seaborn as sns, matplotlib.pyplot as plt
pivot_iptm = df_valid.pivot_table(index="receptor_species",
                                  columns="peptide_species",
                                  values="iptm",
                                  aggfunc="mean")

sns.clustermap(pivot_iptm, cmap="mako", annot=True, fmt=".2f",
               cbar_kws={"label": "Mean IPTM (interface confidence)"})
plt.suptitle("Cross-species MIK2–SCOOP IPTM (AF3)", y=1.02)
plt.savefig(f"{OUT_DIR}/cross_species_heatmap_IPTM_clustered.png", dpi=300, bbox_inches="tight")
plt.close()



# ======================================================
# Z-score normalization per receptor species
# ======================================================
import pandas as pd, seaborn as sns, matplotlib.pyplot as plt, numpy as np

pivot_iptm = df_valid.pivot_table(
    index="receptor_species",
    columns="peptide_species",
    values="iptm",
    aggfunc="mean"
)

# Z-score by row (receptor baseline normalization)
zscore = pivot_iptm.sub(pivot_iptm.mean(axis=1), axis=0).div(pivot_iptm.std(axis=1), axis=0)

plt.figure(figsize=(7,6))
sns.heatmap(zscore, cmap="coolwarm", center=0, annot=True, fmt=".2f",
            cbar_kws={"label": "Z-score (deviation from receptor mean)"})
plt.title("Z-score-normalized cross-species MIK2–SCOOP IPTM (AF3)")
plt.xlabel("Peptide species")
plt.ylabel("Receptor species")
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/cross_species_heatmap_IPTM_zscore.png", dpi=300)
plt.close()
print("Saved z-score plot")


# ========================================================================================================================================================================

# Path to your IQ-TREE distance file
dist_path = "/Users/vishnu/v_files/GT/sem1/Research_sem1/phylogenetic_analysis/ML/MIK2_ECD.fasta.mldist"

phylo_dist = pd.read_csv(dist_path, sep=r"\s+", skiprows=1, header=None)
taxa = phylo_dist.iloc[:, 0].tolist()
mat = phylo_dist.iloc[:, 1:]
mat.index = taxa
mat.columns = taxa[:mat.shape[1]]
phylo_dist = mat

# ===== Map short species codes to full sequence IDs =====
species_map = {
    "C": "C_violacea_Clevi0007s1839___MIK2",
    "D": "D_strictus_Distr0009s35100___LIKE",
    "E": "E_salsugineum_Thhalv10028381m"
}

# ===== Correlate IPTM vs Phylo distances =====
iptm_vals, dist_vals = [], []
for r in pivot_iptm.index:
    for p in pivot_iptm.columns:
        if r != p and r in species_map and p in species_map:
            r_full, p_full = species_map[r], species_map[p]
            if r_full in phylo_dist.index and p_full in phylo_dist.columns:
                iptm_vals.append(pivot_iptm.loc[r, p])
                dist_vals.append(phylo_dist.loc[r_full, p_full])

rho, pval = spearmanr(iptm_vals, dist_vals)
print(f"Spearman correlation (IPTM vs phylogenetic distance): rho={rho:.3f}, p={pval:.3e}")

# ===== Plot =====
sns.regplot(x=dist_vals, y=iptm_vals, scatter_kws={"s":70, "alpha":0.7}, color="k")
plt.xlabel("Phylogenetic distance (ML tree branch length)")
plt.ylabel("Mean IPTM (interface confidence)")
plt.title("Correlation of binding confidence with evolutionary divergence")
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/iptm_vs_phylo_distance_real.png", dpi=300)