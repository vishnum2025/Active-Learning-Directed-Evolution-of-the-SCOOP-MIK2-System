#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Select best 13-mer SCOOP peptides (validated by true length==13)
for the MIK2 cluster-representative species
(E. salsugineum, D. strictus, C. violacea)
grouped by 'Arabidopsis SCOOP' (e.g., 'SCOOP8-cluster'),
picking the lowest E-value (HMMER significance) per (species, Arabidopsis SCOOP).
"""

import pandas as pd
import numpy as np
import os
import re

# ==========================================================
# --- Paths (hardcoded) ---
BASE = "/Users/vishnu/v_files/GT/sem1/Research_sem1"
SD01_PATH = f"{BASE}/sd01.csv"  # your Dataset S1 exported as CSV
OUT_DIR   = f"{BASE}/scoop_mik2_af_pipeline/scoop_cluster_representatives"
os.makedirs(OUT_DIR, exist_ok=True)
# ==========================================================

# -------- Helpers --------
def normalize_cols(cols):
    return [re.sub(r'_+', '_',
                   c.strip()
                    .replace('\u00a0',' ')     # no-break space -> space
                    .replace('(', '')
                    .replace(')', '')
                    .replace('.', '_')
                    .replace(' ', '_'))
            for c in cols]

def find_col(cols, must_include=None, any_of=None):
    """
    Return the first column name whose lowercase matches:
      - contains *all* tokens in must_include
      - contains *any* of tokens in any_of (if provided)
    """
    lc = {c: c.lower() for c in cols}
    must_include = [t.lower() for t in (must_include or [])]
    any_of = [t.lower() for t in (any_of or [])]
    for c, cl in lc.items():
        if all(t in cl for t in must_include) and (not any_of or any(t in cl for t in any_of)):
            return c
    return None

def coerce_float(s):
    # E-values may be strings like "4E-43"
    try:
        return float(s)
    except Exception:
        return np.nan

# -------- Load --------
df = pd.read_csv(SD01_PATH, dtype=str, keep_default_na=False)
df.columns = normalize_cols(df.columns)

# -------- Detect columns robustly --------
col_species = find_col(df.columns, must_include=["species"])
# 13-mer sequence column (often "13mer_alignment_SCOOP" after normalization)
col_13mer   = find_col(df.columns, must_include=["13", "scoop"], any_of=["align", "alignment", "13mer"])
# length column near 13-mer (optional; we will compute length anyway)
col_len13   = find_col(df.columns, must_include=["length", "13"])
# HMMER significance column ("score" in your sheet = E-value)
col_score   = find_col(df.columns, must_include=["score"])
# Arabidopsis SCOOP cluster label (e.g., "SCOOP8-cluster")
col_arab    = find_col(df.columns, must_include=["arabidopsis", "scoop"])

print("Detected columns:")
print(f"  species: {col_species}")
print(f"  13mer  : {col_13mer}")
print(f"  len13  : {col_len13}  (optional)")
print(f"  score  : {col_score}")
print(f"  arab   : {col_arab}")

# Fail-fast if the critical columns are missing
critical_missing = [n for n,v in {
    "species": col_species, "13mer": col_13mer, "score": col_score, "arab": col_arab
}.items() if v is None]
if critical_missing:
    raise SystemExit(f"Missing required columns: {', '.join(critical_missing)}")

# -------- Clean / filter --------
df[col_13mer] = df[col_13mer].astype(str).str.strip()
df[col_arab]  = df[col_arab].astype(str).str.strip()

# Remove rows with NA/empty 13-mer or Arabidopsis SCOOP label
df = df[(df[col_13mer] != "") & (df[col_arab] != "") & (df[col_arab].str.upper() != "NA")]

# Enforce peptide alphabet; compute length directly
df = df[df[col_13mer].str.fullmatch(r"[A-Za-z]+")]
df[col_13mer] = df[col_13mer].str.upper()
df["len_from_seq"] = df[col_13mer].str.len()

# If a dedicated length column exists, trust it only if it equals 13 and matches computed length
if col_len13 is not None:
    # coerce to numeric
    df[col_len13] = pd.to_numeric(df[col_len13], errors="coerce")
    # keep those that are explicitly 13 OR whose computed length is 13
    df = df[(df[col_len13] == 13) | (df["len_from_seq"] == 13)]
else:
    # no explicit length column; rely on computed length
    df = df[df["len_from_seq"] == 13]

# Coerce E-value (“score”) to float
df["E_value"] = df[col_score].apply(coerce_float)
df = df[np.isfinite(df["E_value"])]

# -------- Focus species (cluster representatives’ species) --------
# Normalize species names (strip punctuation, case-insensitive)
df["species_norm"] = (
    df[col_species]
      .astype(str)
      .str.replace(r"[\s._-]", "", regex=True)
      .str.lower()
)

# Targets (normalized keys -> pretty label for filenames)
targets = {
    "esalsugineum": "E_salsugineum",
    "dstrictus":    "D_strictus",
    "cviolacea":    "C_violacea",
}
mask = df["species_norm"].isin(targets.keys())
df = df[mask].copy()
df["species_label"] = df["species_norm"].map(targets)

print(f"Rows after species and 13-mer filtering: {len(df)}")

# -------- Pick best (lowest E-value) per (species, Arabidopsis SCOOP cluster) --------
# Note: many rows share identical 13-mer; we do not deduplicate by sequence here
# because we want one representative per Arabidopsis SCOOP cluster label.
df = df.sort_values(["species_label", col_arab, "E_value"], ascending=[True, True, True])
reps = df.drop_duplicates(subset=["species_label", col_arab], keep="first")

# Tidy selection
reps = reps[["species_label", col_arab, col_13mer, "E_value"]].rename(columns={
    "species_label": "species",
    col_arab: "arabidopsis_scoop_cluster",
    col_13mer: "scoop_13mer",
})

# -------- Write per-species and combined outputs --------
for sp, sub in reps.groupby("species"):
    out_path = os.path.join(OUT_DIR, f"{sp}_13mer_reps.csv")
    sub.sort_values("E_value").to_csv(out_path, index=False)

combined_path = os.path.join(OUT_DIR, "all_cluster_representative_scoops.csv")
reps.sort_values(["species","E_value"]).to_csv(combined_path, index=False)

print("\nSaved:")
for sp in sorted(reps["species"].unique()):
    print(f"  {os.path.join(OUT_DIR, f'{sp}_13mer_reps.csv')}")
print(f"  {combined_path}")

print("\nCounts (unique Arabidopsis SCOOP clusters selected per species):")
print(reps.groupby("species").size().reset_index(name="n_representatives").to_string(index=False))