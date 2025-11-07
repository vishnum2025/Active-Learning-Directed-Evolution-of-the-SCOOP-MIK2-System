#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
	1.	Compute the patristic distance matrix from the phylogenetic tree.
	2.	Perform hierarchical clustering directly on those distances.
	3.	Automatically choose a cutoff height using a “gap” or “tree-cut” method (so we don’t have to guess k).
	4.	Assign cluster IDs purely based on tree structure.
	5.	Select one representative per cluster — the sequence closest to that cluster’s centroid.
	6.	Visualize the clusters in 2D (MDS)
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from sklearn.manifold import MDS

# ------------------------------------------------------
# --- Paths ---
BASE = "/Users/vishnu/v_files/GT/sem1/Research_sem1/phylogenetic_analysis/ML"
TREE_FILE = f"{BASE}/MIK2_rooted.tree"
OUT_FIG = f"{BASE}/a.test_mik2_ml_tree_structure_clusters.png"
OUT_TABLE = f"{BASE}/a.test_mik2_ml_tree_structure_reps.tsv"
# ------------------------------------------------------

# --- 1. Load the phylogenetic tree ---
tree = Phylo.read(TREE_FILE, "newick")
terminals = [t.name for t in tree.get_terminals()]
n = len(terminals)

# --- 2. Compute pairwise patristic distances ---
D = np.zeros((n, n))
for i, t1 in enumerate(terminals):
    for j, t2 in enumerate(terminals):
        if i < j:
            dist = tree.distance(t1, t2)
            D[i, j] = D[j, i] = dist

# --- 3. Hierarchical clustering ---
Z = linkage(D, method="average")

# --- 4. Choose adaptive cutoff (TreeCut-like) ---
# Compute the first large gap in linkage distances
diffs = np.diff(Z[:, 2])
cut_height = np.median(Z[:, 2]) + np.std(Z[:, 2])
# cutoff is placed one sd above the median merge height 
labels = fcluster(Z, t=cut_height, criterion="distance")

print(f"Adaptive cut height ≈ {cut_height:.3f}, yielding {len(set(labels))} clusters.")

# --- 5. Create DataFrame ---
df = pd.DataFrame({"name": terminals, "cluster": labels})

# --- 6. Compute MDS projection for visualization ---
mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
emb = mds.fit_transform(D)
df["x"], df["y"] = emb[:, 0], emb[:, 1]

# # --- 7. Compute centroid representative for each cluster ---
# centroid_reps = []
# for cid, group in df.groupby("cluster"):
#     cx, cy = group["x"].mean(), group["y"].mean()
#     group["dist_to_centroid"] = np.sqrt((group["x"] - cx)**2 + (group["y"] - cy)**2)
#     rep = group.loc[group["dist_to_centroid"].idxmin(), "name"]
#     members = group["name"].tolist()
#     centroid_reps.append((cid, rep, len(group)))

#     print(f"\n Cluster {cid}: {len(group)} members → representative: {rep}")
#     print("   Members:")
#     for m in members:
#         print(f"     - {m}")


# # --- 7b. Compute how far A_thaliana_AT4G08850 is from its cluster centroid ---
# target = "A_thaliana_AT4G08850"
# if target in df["name"].values:
#     row = df[df["name"] == target].iloc[0]
#     cluster_id = row["cluster"]
#     cluster_group = df[df["cluster"] == cluster_id]
#     cx, cy = cluster_group["x"].mean(), cluster_group["y"].mean()
#     dist_target = np.sqrt((row["x"] - cx)**2 + (row["y"] - cy)**2)
#     print(f"\n {target} belongs to Cluster {cluster_id} and is {dist_target:.4f} units from its centroid.")
# else:
#     print(f"\n {target} not found in terminal names.")

# # --- 7c. Compute how far E_salsugineum_Thhalv10028381m is from its cluster centroid ---
# target2 = "E_salsugineum_Thhalv10028381m"
# if target2 in df["name"].values:
#     row2 = df[df["name"] == target2].iloc[0]
#     cluster_id2 = row2["cluster"]
#     cluster_group2 = df[df["cluster"] == cluster_id2]
#     cx2, cy2 = cluster_group2["x"].mean(), cluster_group2["y"].mean()
#     dist_target2 = np.sqrt((row2["x"] - cx2)**2 + (row2["y"] - cy2)**2)
#     print(f"{target2} belongs to Cluster {cluster_id2} and is {dist_target2:.4f} units from its centroid.")
# else:
#     print(f"{target2} not found in terminal names.")



# --- 7. Compute centroid representative for each cluster ---
centroid_reps = []
for cid, group in df.groupby("cluster"):
    cx, cy = group["x"].mean(), group["y"].mean()
    group["dist_to_centroid"] = np.sqrt((group["x"] - cx)**2 + (group["y"] - cy)**2)
    
    # Force A_thaliana_AT4G08850 as cluster center for cluster 1
    if cid == 1 and "A_thaliana_AT4G08850" in group["name"].values:
        rep = "A_thaliana_AT4G08850"
    else:
        rep = group.loc[group["dist_to_centroid"].idxmin(), "name"]
    
    members = group["name"].tolist()
    centroid_reps.append((cid, rep, len(group)))

    print(f"\nCluster {cid}: {len(group)} members → representative: {rep}")
    print("   Members:")
    for m in members:
        print(f"     - {m}")

# --- 8. Save summary table ---
summary = pd.DataFrame(centroid_reps, columns=["cluster", "representative", "size"])
summary.to_csv(OUT_TABLE, sep="\t", index=False)



# --- 9. Visualization: publication-quality with cluster ellipses (fixed small clusters) ---

from matplotlib.patches import Ellipse

sns.set_theme(style="whitegrid", context="talk", font_scale=1.3)
plt.figure(figsize=(10, 8), dpi=400)
ax = plt.gca()

# Define cluster colors and markers
unique_clusters = sorted(df["cluster"].unique())
palette = sns.color_palette("Set2", n_colors=len(unique_clusters))
markers = ["o", "s", "^", "D", "P", "v", "X"]

# Scatter each cluster with marker and color
for i, cid in enumerate(unique_clusters):
    subset = df[df["cluster"] == cid]
    ax.scatter(
        subset["x"], subset["y"],
        s=90,
        alpha=0.9,
        marker=markers[i % len(markers)],
        color=palette[i],
        edgecolor="black",
        linewidth=0.5,
        label=f"Cluster {cid}"
    )

    # Draw ellipse around each cluster (works even for small clusters)
    # Draw ellipse around each cluster (robust for small clusters)
    if len(subset) >= 3:
        cov = np.cov(subset[["x", "y"]].T)
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        vals, vecs = vals[order], vecs[:, order]
        theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
        width, height = 2 * np.sqrt(vals * 5.991)  # 95% CI
        ellipse = Ellipse(
            xy=(subset["x"].mean(), subset["y"].mean()),
            width=width,
            height=height,
            angle=theta,
            color=palette[i],
            alpha=0.25,
            lw=1.5
        )
        ax.add_patch(ellipse)

    elif len(subset) == 2:
        # Compute midpoint and distance between two points
        p1, p2 = subset[["x", "y"]].values
        midx, midy = (p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2
        dist = np.linalg.norm(p1 - p2)
        ellipse = Ellipse(
            xy=(midx, midy),
            width=dist * 1.3,    # modest elongation
            height=dist * 0.5,   # prevent thin line
            angle=np.degrees(np.arctan2(p2[1] - p1[1], p2[0] - p1[0])),
            color=palette[i],
            alpha=0.25,
            lw=1.5
        )
        ax.add_patch(ellipse)

    else:
        # Singleton cluster: small round bubble
        circ = Ellipse(
            xy=(subset["x"].iloc[0], subset["y"].iloc[0]),
            width=0.05,
            height=0.05,
            angle=0,
            color=palette[i],
            alpha=0.3,
            lw=1.0
        )
        ax.add_patch(circ)

# Highlight cluster representatives
for cid, rep, _ in centroid_reps:
    row = df[df["name"] == rep]
    ax.scatter(
        row["x"], row["y"],
        s=400,
        facecolors='none',
        edgecolors='black',
        linewidths=2.4,
        zorder=5
    )
    ax.text(
        row["x"] + 0.015, row["y"], rep,
        fontsize=13,  # enlarged for clarity
        fontweight="bold",
        ha="left",
        va="center",
        color="black",
        bbox=dict(facecolor="white", alpha=0.65, edgecolor="none", pad=1)
    )

# Style and labels
ax.set_facecolor("#f2f2f2")
sns.despine(left=False, bottom=False)
ax.grid(True, color="white", lw=1.0)
ax.set_title("MIK2 Phylogenetic Clusters (MDS Projection)", fontsize=18, weight="bold", pad=20)
ax.set_xlabel("MDS Dimension 1", fontsize=14)
ax.set_ylabel("MDS Dimension 2", fontsize=14)

# Legend formatting
leg = ax.legend(
    title="Cluster",
    title_fontsize=12,
    fontsize=11,
    frameon=True,
    loc="upper right",
    bbox_to_anchor=(1.18, 1.02),
    borderpad=0.8
)
leg.get_frame().set_edgecolor("gray")
leg.get_frame().set_linewidth(0.8)
leg.get_frame().set_alpha(0.9)

plt.tight_layout()

plt.savefig(OUT_FIG, dpi=600, bbox_inches="tight", transparent=False)
plt.show()