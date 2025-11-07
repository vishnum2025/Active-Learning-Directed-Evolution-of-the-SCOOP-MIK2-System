#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Analyze ALDE peptide origins vs AF3-predicted fitness.
------------------------------------------------------
This script:
 1. Reconstructs approximate "mutational" (near-native) vs "random" (de novo) labels
 2. Merges them with AF3 iPTM predictions
 3. Plots distribution of predicted iPTM scores per origin category

Outputs:
  - analysis_round1/iptm_origin_comparison.csv
  - analysis_round1/iptm_origin_comparison.png
"""

import os, itertools
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------- PATHS ----------------
BASE_PATH = "/Users/vishnu/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/alde_athaliana"
DOMAIN_FILE = f"{BASE_PATH}/alde_domain.csv"
TRAIN_FILE = f"{BASE_PATH}/alde_training_round0.csv"
AF3_SUMMARY = f"{BASE_PATH}/analysis_round1/round1_af3_summary.csv"
OUT_DIR = f"{BASE_PATH}/analysis_round1"
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------- CONSTANTS ----------------
AA20 = list("ACDEFGHIKLMNPQRSTVWY")
L = 13
MOTIF_POS = (4, 6)  # 0-based indices for SxS motif

# ---------------- STEP 1: LOAD INPUTS ----------------
print("[INFO] Loading data...")
domain = pd.read_csv(DOMAIN_FILE)
train = pd.read_csv(TRAIN_FILE)
af3 = pd.read_csv(AF3_SUMMARY)

seed_peptides = set(train["sequence"].str.strip().str.upper())
domain["sequence"] = domain["sequence"].str.strip().str.upper()

# ---------------- STEP 2: RECONSTRUCT MUTATIONAL SET ----------------
print("[INFO] Reconstructing approximate mutational neighborhood...")

mutational_set = set()

for base in seed_peptides:
    if len(base) != L:
        continue
    base = list(base)
    # enforce motif
    base[MOTIF_POS[0]], base[MOTIF_POS[1]] = "S", "S"

    # Single mutants
    for pos in range(L):
        if pos in MOTIF_POS:
            continue
        for aa in AA20:
            mut = base.copy()
            mut[pos] = aa
            mutational_set.add("".join(mut))

    # Double mutants (optional, uncomment if you want broader coverage)
    # for pos1, pos2 in itertools.combinations([i for i in range(L) if i not in MOTIF_POS], 2):
    #     for aa1, aa2 in itertools.product(AA20, repeat=2):
    #         mut = base.copy()
    #         mut[pos1], mut[pos2] = aa1, aa2
    #         mutational_set.add("".join(mut))

print(f"[INFO] Reconstructed ~{len(mutational_set):,} mutational variants (approx.)")

# ---------------- STEP 3: TAG ORIGIN ----------------
domain["origin"] = domain["sequence"].apply(
    lambda s: "mutational" if s in mutational_set else "random"
)

# ---------------- STEP 4: MERGE WITH AF3 PREDICTIONS ----------------
merged = af3.merge(domain, on="sequence", how="left")
merged = merged.dropna(subset=["iptm"])
merged.to_csv(os.path.join(OUT_DIR, "iptm_origin_comparison.csv"), index=False)
print(f"[OK] Saved merged summary → {OUT_DIR}/iptm_origin_comparison.csv")

# ---------------- STEP 5: PLOT DISTRIBUTIONS (KDE SCALED TO COUNTS) ----------------
plt.style.use("ggplot")
sns.set_context("talk", font_scale=1.2)
sns.set_style("whitegrid", {"axes.facecolor": "#F5F5F5"})

# Count weights per sample — ensures KDE integrates to total count, not 1.0
merged["weight"] = (
    merged.groupby("origin")["origin"].transform(lambda x: 1 / len(x))
)

plt.figure(figsize=(9, 6))
ax = sns.kdeplot(
    data=merged,
    x="iptm",
    hue="origin",
    fill=True,
    common_norm=False,
    alpha=0.35,
    linewidth=3.0,
    bw_adjust=0.6,
    palette={"mutational": "#2530C4", "random": "#0CC3D0"},
    multiple="layer",
    # convert normalized KDE to counts:
    weights=merged["origin"].map(
        merged["origin"].value_counts()
    ).astype(float),
)

# ---- Styling ----
plt.xlabel("iPTM", fontsize=15)
plt.ylabel("Count", fontsize=15)
plt.title("Distribution of Predicted Peptides by Origin", fontsize=17, weight="semibold")
plt.legend(title="Domain Source", labels=["Mutational (≈3.9k)", "Random (≈8k)"], fontsize=13, title_fontsize=13)
plt.xlim(0.4, 1.0)

# Scale y-axis up a bit for readability
ymax = ax.get_ylim()[1]
plt.ylim(0, ymax * 1.25)

# Polishing
plt.gcf().patch.set_facecolor("#F5F5F5")
ax.set_facecolor("#F5F5F5")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()

out_plot = os.path.join(OUT_DIR, "iptm_origin_comparison_kde_counts.png")
plt.savefig(out_plot, dpi=400, bbox_inches="tight", facecolor="#F5F5F5")
plt.close()

print(f"[OK] Saved count-scaled KDE plot → {out_plot}")