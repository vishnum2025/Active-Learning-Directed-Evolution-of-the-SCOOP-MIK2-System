#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Clean analysis of AF3 SCOOP–MIK2 results:
 - Summarizes iPTM/PTM/ranking_score from all modeled complexes
 - Compares distributions between input (training) and predicted (ALDE Round 1)
 - Saves clean minimal histograms + summary CSV in analysis folder
"""

import os, json, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ---------------- PATHS ----------------
BASE = "/Users/vishnu/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/alde_athaliana"
AF3_DIR = f"{BASE}/round1_scoop_outputs"
TRAIN_FILE = f"{BASE}/alde_training_round0.csv"
OUT_DIR = f"{BASE}/analysis_round1"
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------- HELPERS ----------------
def try_load_json(path):
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:
        return None

def extract_sequence(folder_name):
    """Extract peptide sequence from folder name."""
    if "__" in folder_name:
        return folder_name.split("__", 1)[-1]
    m = re.search(r"[A-Z]{10,}", folder_name)
    return m.group(0) if m else folder_name

def summarize_folder(folder_path):
    """Extract iPTM/PTM/ranking_score from *_summary_confidences.json or ranking_scores.csv."""
    result = {"iptm": np.nan, "ptm": np.nan, "ranking_score": np.nan}

    for f in os.listdir(folder_path):
        if f.endswith("_summary_confidences.json"):
            data = try_load_json(os.path.join(folder_path, f))
            if data:
                result["iptm"] = data.get("iptm", np.nan)
                result["ptm"] = data.get("ptm", np.nan)
                result["ranking_score"] = data.get("ranking_score", np.nan)
                return result

    for f in os.listdir(folder_path):
        if f.endswith("_ranking_scores.csv"):
            try:
                df = pd.read_csv(os.path.join(folder_path, f))
                for col in ["iptm", "ptm", "ranking_score"]:
                    if col in df.columns:
                        result[col] = df[col].mean()
                return result
            except Exception:
                pass
    return result

# ---------------- MAIN ----------------
records = []

for name in sorted(os.listdir(AF3_DIR)):
    outer = os.path.join(AF3_DIR, name)
    if not os.path.isdir(outer):
        continue

    # Handle nested dir structure
    inner = os.path.join(outer, name)
    target_dir = inner if os.path.isdir(inner) else outer

    seq = extract_sequence(name).strip().upper()
    scores = summarize_folder(target_dir)
    records.append({
        "sequence": seq,
        "iptm": scores["iptm"],
        "ptm": scores["ptm"],
        "ranking_score": scores["ranking_score"]
    })

df_pred = pd.DataFrame(records)
df_pred["iptm+ptm"] = df_pred[["iptm", "ptm"]].sum(axis=1, min_count=1)
df_pred = df_pred.sort_values("iptm", ascending=False).reset_index(drop=True)

# Load training data
df_train = pd.read_csv(TRAIN_FILE)
df_train = df_train.rename(columns={"fitness": "iptm"})  # unify naming for easy plotting

# Save combined summary CSV
summary_path = os.path.join(OUT_DIR, "round1_af3_summary.csv")
df_pred.to_csv(summary_path, index=False)
print(f"[OK] Saved AF3 summary → {summary_path}")

# -------------------------------- PLOTTING  --------------------------------

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches

plt.style.use("ggplot")
sns.set_context("talk")
sns.set_style("whitegrid", {"axes.facecolor": "#F5F5F5"})

bins = np.linspace(0.4, 1.0, 15)
highlight_threshold = 0.8

# Shared y-axis limit
ymax = max(
    np.histogram(df_train["iptm"].dropna(), bins=bins)[0].max(),
    np.histogram(df_pred["iptm"].dropna(), bins=bins)[0].max()
) * 1.1

fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True, gridspec_kw={"hspace": 0.35})

def plot_hist_with_hatch(ax, data, title, base_color):
    counts, edges = np.histogram(data, bins=bins)

    for i in range(len(edges) - 1):
        bin_center = (edges[i] + edges[i + 1]) / 2
        # Add hatch only for iPTM ≥ threshold
        hatch = "///" if bin_center >= highlight_threshold else None
        ax.bar(
            edges[i], counts[i],
            width=(edges[1] - edges[0]),
            align="edge",
            color=base_color,
            edgecolor="#FFFFFF",
            alpha=0.9,
            hatch=hatch
        )

    ax.set_title(title, fontsize=15, weight="semibold")
    ax.set_ylabel("Count", fontsize=13)
    ax.set_ylim(0, ymax)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# Plot Training (Input)
plot_hist_with_hatch(
    axes[0],
    df_train["iptm"].dropna(),
    "Distribution of Input (Training) Peptides – Fitness (iPTM)",
    base_color="#6BAE44"   # original green
)

# Plot Predicted (ALDE)
plot_hist_with_hatch(
    axes[1],
    df_pred["iptm"].dropna(),
    "Distribution of Predicted (ALDE) Peptides – Fitness (iPTM)",
    base_color="#4B72B0"   # original blue
)
axes[1].set_xlabel("iPTM", fontsize=13)

# ---- Add Legend for Hatching ----
import matplotlib.patches as mpatches

# Create hatch legend patch
hatch_patch = mpatches.Patch(
    facecolor="none", 
    edgecolor="black", 
    hatch="///", 
    label="iPTM ≥ 0.8"
)

# ---- Add legend to top-right INSIDE each subplot (below title) ----
for ax in axes:
    ax.legend(
        handles=[hatch_patch],
        title="Highlighted Region",
        loc="upper right",
        bbox_to_anchor=(1, 0.93),   # Slightly below title area
        fontsize=11,
        title_fontsize=11,
        frameon=True
    )  # Adjust for padding against edges

# ---- Final polish ----
fig.patch.set_facecolor("#F5F5F5")
for ax in axes:
    ax.set_facecolor("#F5F5F5")
    ax.tick_params(axis="both", labelsize=12)
    ax.set_xlim(0.4, 1.0)

plt.tight_layout()
out_path = os.path.join(OUT_DIR, "iptm_comparison_hist_hatched.png")
plt.savefig(out_path, dpi=400, bbox_inches="tight")
plt.close()

print(f"[OK] Saved hatched-highlight plot → {out_path}")