# #!/usr/bin/env python3

# ============================================================================
# Alignment
# ============================================================================

import matplotlib.pyplot as plt
import numpy as np

WT = "PVRSSQSSQAGGR"

PEPTIDES = [
    "AVDHSSSPRHGQR",
    "PVRGSHSSHLPVR",
    "RVGRSKSLRGPQK",
    "GVDHSYSQRAGQR",
    "AVDHSYSTRNGQR",
    "AVDHSVSARAGQR",
    "RVHRSSSKRGPQK",
    "AVDRSYSPRGGQR",
    "TVRSSRSQIHPQK",
    "YVQRSHSSHTGSR",
    "FTGPSYSGHGGAR",
    "HVRSSRSQIGIQK"
]

# ---- simple family grouping (basic rich vs mixed) ----
def basic_fraction(seq):
    return sum(a in "KRH" for a in seq) / len(seq)

basic_family = [p for p in PEPTIDES if basic_fraction(p) > 0.3]
mixed_family = [p for p in PEPTIDES if basic_fraction(p) <= 0.3]

ordered_peptides = basic_family + mixed_family

labels = ["WT"] + [f"P{i+1}" for i in range(len(ordered_peptides))]
sequences = [WT] + ordered_peptides

fig, ax = plt.subplots(figsize=(13,8))

n_rows = len(sequences)
n_cols = len(WT)

ax.set_xlim(0, n_cols)
ax.set_ylim(0, n_rows)
ax.invert_yaxis()

ax.set_xticks(range(n_cols))
ax.set_xticklabels(range(1, n_cols+1), fontsize=12)
ax.set_yticks(range(n_rows))
ax.set_yticklabels(labels, fontsize=12)

# --- region shading ---
# N-terminal diversification
ax.axvspan(0, 4, color="lightgrey", alpha=0.3)

# SxS core
ax.axvspan(4, 7, color="lightgreen", alpha=0.3)

# C-terminal binding region
ax.axvspan(7, 13, color="blue", alpha=0.15)

# --- plot residues ---
for row, seq in enumerate(sequences):
    for col, aa in enumerate(seq):
        is_mut = (row > 0 and aa != WT[col])

        color = "black"
        if aa in "KRH":
            color = "#06849e"
        elif aa in "STNQ":
            color = "#2ca02c"
        elif aa in "DE":
            color = "#d62728"

        if is_mut:
            color = "red"

        weight = "bold" if col in [4,6] else "normal"

        ax.text(
            col + 0.5,
            row + 0.5,
            aa,
            ha="center",
            va="center",
            fontsize=16,
            color=color,
            fontweight=weight
        )

        # subtle mutation dot
        if is_mut:
            ax.plot(col + 0.5, row + 0.85, "o", color="red", markersize=2)

ax.set_title(
    "Evolutionary Patterns in ALDE-Designed SCOOP Peptides\n"
    "Grey = N-terminal diversification | Green = Conserved SxS core | Blue = Binding region\n"
    "Red letters = mutation vs WT",
    fontsize=16
)

ax.tick_params(length=0)
ax.grid(False)

plt.tight_layout()
plt.savefig("/Users/vishnu/v_files/GT/Research/spring_26/alde_peptide_alignment_aesthetic.png", dpi=300)
plt.close()

print("Saved → alde_peptide_alignment_aesthetic.png")




# ============================================================================
# Entropy
# ============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import log2

# ---------------- PATHS ----------------
TRAIN_PATH = "/Users/vishnu/v_files/GT/Research/scoop_mik2_af_pipeline/alde_athaliana/alde_training_round0.csv"
ROUND1_PATH = "/Users/vishnu/v_files/GT/Research/scoop_mik2_af_pipeline/alde_athaliana/analysis_round1/round1_af3_summary.csv"
# ---------------------------------------

AA20 = list("ACDEFGHIKLMNPQRSTVWY")

def compute_entropy(seqs):
    L = len(seqs[0])
    entropies = []
    for i in range(L):
        column = [s[i] for s in seqs]
        probs = []
        for aa in AA20:
            p = column.count(aa) / len(column)
            if p > 0:
                probs.append(p)
        H = -sum(p * log2(p) for p in probs)
        entropies.append(H)
    return np.array(entropies)

# Load
train_seqs = pd.read_csv(TRAIN_PATH)["sequence"].astype(str).tolist()
r1_df = pd.read_csv(ROUND1_PATH)

# Only keep high binders (optional but smarter biologically)
r1_df = r1_df[r1_df["iptm"] >= 0.8]
r1_seqs = r1_df["sequence"].astype(str).tolist()

# Sanity
L = len(train_seqs[0])
assert all(len(s)==L for s in train_seqs+r1_seqs)

# Entropy
H_train = compute_entropy(train_seqs)
H_r1 = compute_entropy(r1_seqs)
dH = H_r1 - H_train

positions = np.arange(1, L+1)

plt.figure(figsize=(12,5))
plt.plot(positions, H_train, label="Training", linewidth=2)
plt.plot(positions, H_r1, label="Round1 (binders)", linewidth=2)
plt.xlabel("Position")
plt.ylabel("Shannon Entropy (bits)")
plt.legend()
plt.title("Entropy Landscape")
plt.tight_layout()
plt.show()

plt.figure(figsize=(12,5))
plt.bar(positions, dH)
plt.axhline(0, linestyle="--")
plt.xlabel("Position")
plt.ylabel("ΔEntropy (Round1 - Training)")
plt.title("Entropy Shift After ALDE Selection")
plt.tight_layout()
plt.show()




# ============================================================================
# Clustering
# ============================================================================

#!/usr/bin/env python3

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from collections import Counter

# --------------------------
# INPUT: your peptides here
# --------------------------
peptides = [
    "AVDHSSSPRHGQR",
    "PVRGSHSSHLPVR",
    "RVGRSKSLRGPQK",
    "GVDHSYSQRAGQR",
    "AVDHSYSTRNGQR",
    "AVDHSVSARAGQR",
    "RVHRSSSKRGPQK",
    "AVDRSYSPRGGQR",
    "TVRSSRSQIHPQK",
    "YVQRSHSSHTGSR",
    "FTGPSYSGHGGAR",
    "HVRSSRSQIGIQK"
]

L = len(peptides[0])
N = len(peptides)

# --------------------------
# 1. Pairwise Hamming Matrix
# --------------------------
def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))

dist_matrix = np.zeros((N, N))

for i in range(N):
    for j in range(N):
        dist_matrix[i, j] = hamming(peptides[i], peptides[j])

dist_df = pd.DataFrame(dist_matrix, index=peptides, columns=peptides)

plt.figure(figsize=(10,8))
sns.heatmap(dist_df, cmap="viridis")
plt.title("Pairwise Hamming Distance")
plt.tight_layout()
plt.show()

# --------------------------
# 2. Hierarchical Clustering
# --------------------------
# Convert full distance matrix to condensed form properly
from scipy.spatial.distance import squareform

# dist_matrix already computed above
condensed = squareform(dist_matrix)

Z = linkage(condensed, method="average")

plt.figure(figsize=(8,5))
dendrogram(Z, labels=peptides, leaf_rotation=90)
plt.title("Hierarchical Clustering (Hamming)")
plt.tight_layout()
plt.show()

# choose 2–3 clusters manually
clusters = fcluster(Z, t=3, criterion="maxclust")

cluster_df = pd.DataFrame({
    "sequence": peptides,
    "cluster": clusters
})

print("\nCluster assignments:")
print(cluster_df.sort_values("cluster"))


