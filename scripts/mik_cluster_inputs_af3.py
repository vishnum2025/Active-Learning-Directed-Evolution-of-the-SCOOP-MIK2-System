#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import pandas as pd
from Bio import SeqIO

# -------------------- PATHS (hard-coded) --------------------
BASE = "/Users/vishnu/v_files/GT/sem1/Research_sem1"
MIK2_FASTA = f"{BASE}/sd03.txt"
SCOOP_CSV  = f"{BASE}/scoop_mik2_af_pipeline/scoop_cluster_representatives/all_cluster_representative_scoops.csv"
OUT_DIR    = f"{BASE}/scoop_mik2_af_pipeline/scoop_cluster_representatives/af3_inputs"
os.makedirs(OUT_DIR, exist_ok=True)

# Receptor IDs (cluster representatives) exactly as present in sd03.txt
RECEPTORS = {
    "E_salsugineum": "E_salsugineum_Thhalv10028381m",
    "D_strictus":   "D_strictus_Distr0009s35100___LIKE",
    "C_violacea":   "C_violacea_Clevi0007s1839___MIK2",
}

# Target SCOOP species “keys” for filtering the master CSV
SCOOP_SPECIES_KEYS = {
    "E_salsugineum": ["Esalsugineum"],
    "D_strictus":    ["Dstrictus"],
    "C_violacea":    ["Cviolacea"],
}

# -------------------- ECD TRIMMER (your logic) --------------------
def find_first_long_hydrophobic_run(seq: str, window=21, min_frac=0.75, start_region_frac=0.6):
    hydro = set("VILMFWYAGC")
    start_scan_i = int(len(seq) * start_region_frac)
    for i in range(start_scan_i, max(len(seq) - window + 1, start_scan_i + 1)):
        frag = seq[i:i+window]
        if sum(aa in hydro for aa in frag) >= int(min_frac * window):
            return i
    return None

def trim_ecd(seq: str, start=None, end=None):
    if start is None:
        start = 24  # skip SP region crudely
    if end is None:
        tm = find_first_long_hydrophobic_run(seq, window=21, min_frac=0.75, start_region_frac=0.6)
        end = tm if tm is not None else min(len(seq), 710)
    end = max(start, min(end, len(seq)))
    return seq[start:end], (start, end)

# -------------------- LOAD INPUTS --------------------
# MIK2 homologs
mik2_records = {rec.id.split()[0]: str(rec.seq) for rec in SeqIO.parse(MIK2_FASTA, "fasta")}

# SCOOP representatives
df = pd.read_csv(SCOOP_CSV)
df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

# Expect columns: species, arabidopsis_scoop_cluster, scoop_13mer, e_value
need = {"species","arabidopsis_scoop_cluster","scoop_13mer"}
missing = need - set(df.columns)
if missing:
    raise SystemExit(f"Missing columns in SCOOP CSV: {missing}")

# Clean and enforce 13-mer peptides
df["scoop_13mer"] = df["scoop_13mer"].astype(str).str.strip().str.upper()
df = df[df["scoop_13mer"].str.fullmatch(r"[A-Z]{13}")]

# Normalize species token to match keys
def norm(s): return re.sub(r"[\s._-]", "", str(s)).lower()
df["species_norm"] = df["species"].apply(norm)

# Build per-species SCOOP tables
scoop_by_species = {}
for label, keys in SCOOP_SPECIES_KEYS.items():
    mask = df["species_norm"].apply(lambda s: any(k.lower() in s for k in keys))
    sub = df[mask].drop_duplicates(subset=["arabidopsis_scoop_cluster","scoop_13mer"])
    scoop_by_species[label] = sub

# -------------------- GENERATE AF3 INPUT FASTAS --------------------
rows = []
n_written = 0

for rec_label, rec_id in RECEPTORS.items():
    if rec_id not in mik2_records:
        print(f"[WARN] receptor {rec_id} not found in sd03; skipping")
        continue

    full = mik2_records[rec_id]
    ecd, (s_i, e_i) = trim_ecd(full, start=None, end=None)

    # Pair this receptor with *all* SCOOPs from all three species
    for spp_label, scoops in scoop_by_species.items():
        for _, r in scoops.iterrows():
            pep = r["scoop_13mer"]
            cl  = r["arabidopsis_scoop_cluster"]

            # FASTA single-line “peptide:receptor” format (your AF3 convention)
            # Header: >{recLabel}__{sppLabel}_{cluster}
            out_name = f"{rec_label}__{spp_label}_{cl}.fasta"
            out_path = os.path.join(OUT_DIR, out_name)

            with open(out_path, "w") as f:
                f.write(f">{rec_label}__{spp_label}_{cl}\n")
                f.write(f"{pep}:{ecd}\n")

            rows.append({
                "receptor_label": rec_label,
                "receptor_id": rec_id,
                "ecd_slice": f"{s_i+1}-{e_i}",  # 1-based
                "scoop_species": spp_label,
                "scoop_cluster": cl,
                "peptide": pep,
                "fasta": out_name
            })
            n_written += 1

print(f"Wrote {n_written} AF3 FASTAs → {OUT_DIR}")

# Summary table
out_table = os.path.join(OUT_DIR, "af3_input_manifest.csv")
pd.DataFrame(rows).to_csv(out_table, index=False)
print(f"Manifest: {out_table}")