#!/usr/bin/env python3
# scripts/cluster_diff_and_design.py
#
# - Compares WT pocket(S5/S7) vs CxC/TxT S5/S7 contact clusters.
# - Produces differential rankings and a small, principled mutational panel:
#     * Off-target-only sites  -> propose {A,E,K} to disrupt
#     * Pocket (WT) sites     -> propose {N,Q,Y,H} to reinforce S/T interactions
# - Adds receptor residue amino acid identity in all outputs.
#
# All paths are hard-coded for your machine.

import os, math
import pandas as pd
import numpy as np

# ================== PATHS (EDIT IF NEEDED) ==================
BASE = "/Users/vishnu/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline"
CONTACT_DIR = os.path.join(BASE, "analysis_contacts")

POCKET_WT = os.path.join(CONTACT_DIR, "pocket_cluster_WT_S5S7.csv")
OFF_CXC   = os.path.join(CONTACT_DIR, "offtarget_cluster_CxC_S5S7.csv")
OFF_TXT   = os.path.join(CONTACT_DIR, "offtarget_cluster_TxT_S5S7.csv")

# Baseline complex used to emit AF3 variants: expects "peptide:receptor" on one line after header.
BASE_FASTA = "/Users/vishnu/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/CxC_TxT/AF3_Cviol1838-SCOOP12.fasta"
OUT_DIR    = os.path.join(BASE, "designs")
os.makedirs(OUT_DIR, exist_ok=True)

# Known pocket residues in 1838 model indexing (paper = model + 68):
POCKET_MODEL = [178, 200, 224, 226, 248]   # Do not change
# Conservative subset to start:
POCKET_MODEL_PRIMARY = [224, 248]

# Panels
OFFTARGET_DISRUPT_PANEL = ["A","E","K"]
POCKET_REINFORCE_PANEL  = ["N","Q","Y","H"]

# How many top off-target-only sites to disrupt per variant:
N_OFFTARGET_PER_VARIANT = 8


# ================== UTILS ==================
def read_csv_safe(p):
    if not os.path.exists(p):
        return pd.DataFrame()
    try:
        return pd.read_csv(p)
    except Exception as e:
        print(f"[WARN] could not read {p}: {e}")
        return pd.DataFrame()

def prep(df, tag):
    if df.empty:
        return df
    need = {"resnum_model","n_contacts","min_dmin","pep_positions","pep_aas","aa"}
    missing = need - set(df.columns)
    if missing:
        print(f"[WARN] {tag} missing columns: {missing}")
    df = df.copy()
    df["source"] = tag
    df["resnum_model"] = df["resnum_model"].astype(int)
    df["n_contacts"]   = pd.to_numeric(df["n_contacts"], errors="coerce").fillna(0.0)
    df["min_dmin"]     = pd.to_numeric(df["min_dmin"], errors="coerce")
    return df


def rank_offtarget_only(df_wt, df_var):
    """Return df of residues present in variant but absent/weak in WT, ranked."""
    w = df_wt.set_index("resnum_model")
    v = df_var.set_index("resnum_model")
    common = set(w.index) & set(v.index)
    v_only = set(v.index) - set(w.index)

    rows=[]
    for rn in sorted(set(v.index)):
        n_v = float(v.loc[rn, "n_contacts"]) if rn in v.index else 0.0
        d_v = float(v.loc[rn, "min_dmin"])   if rn in v.index else np.nan
        n_w = float(w.loc[rn, "n_contacts"]) if rn in w.index else 0.0
        d_w = float(w.loc[rn, "min_dmin"])   if rn in w.index else np.nan
        if "aa" in v.columns and rn in v.index:
            aa_v = str(v.loc[rn, "aa"])
        else:
            aa_v = "X"

        if "aa" in w.columns and rn in w.index:
            aa_w = str(w.loc[rn, "aa"])
        else:
            aa_w = "X"

        near_bonus = 1.0/(d_v+1e-6) if not math.isnan(d_v) else 0.0
        wt_pen     = n_w * (1.0/(d_w+1e-6) if not math.isnan(d_w) else 0.0)
        score = n_v*near_bonus - 0.5*wt_pen

        rows.append({
            "resnum_model": rn,
            "aa_variant": aa_v,
            "aa_WT": aa_w,
            "n_contacts_variant": n_v,
            "min_dmin_variant": d_v,
            "n_contacts_WT": n_w,
            "min_dmin_WT": d_w,
            "variant_only": rn in v_only,
            "score": score
        })
    out = pd.DataFrame(rows).sort_values("score", ascending=False)
    return out


def rank_pocket_reinforce(df_wt):
    """Return WT pocket residues ranked by how consistently close/frequent they are (S5/S7)."""
    if df_wt.empty:
        return pd.DataFrame()
    df = df_wt.copy()
    df["reinforce_score"] = df["n_contacts"] * (1.0/(df["min_dmin"]+1e-6))
    df = df.sort_values(["reinforce_score","n_contacts"], ascending=[False, False])
    return df


def read_baseline_complex(fasta_path):
    name, seq = None, []
    with open(fasta_path) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith(">"):
                name = line[1:].strip()
            else:
                seq.append(line)
    if name is None:
        raise RuntimeError("No FASTA header found")
    full = "".join(seq).strip()
    if ":" not in full:
        raise RuntimeError("Expected 'peptide:receptor' in baseline FASTA")
    pep, rec = full.split(":",1)
    return name, pep, rec


def mutate(rec, pos_model, new_aa):
    i = pos_model - 1
    if not (0 <= i < len(rec)):
        raise ValueError(f"Model pos {pos_model} out of range for receptor length {len(rec)}")
    return rec[:i] + new_aa + rec[i+1:]


# ================== MAIN ==================
def main():
    wt = prep(read_csv_safe(POCKET_WT), "WT")
    cxc= prep(read_csv_safe(OFF_CXC), "CxC")
    txt= prep(read_csv_safe(OFF_TXT), "TxT")

    if wt.empty:
        print("[ERROR] WT pocket cluster file missing/empty.")
        return
    if cxc.empty and txt.empty:
        print("[ERROR] Variant cluster files missing/empty.")
        return

    if not cxc.empty:
        diff_cxc = rank_offtarget_only(wt, cxc)
        diff_cxc.to_csv(os.path.join(OUT_DIR, "diff_offtarget_only_CxC.csv"), index=False)
    if not txt.empty:
        diff_txt = rank_offtarget_only(wt, txt)
        diff_txt.to_csv(os.path.join(OUT_DIR, "diff_offtarget_only_TxT.csv"), index=False)

    pocket_rank = rank_pocket_reinforce(wt)
    pocket_rank.to_csv(os.path.join(OUT_DIR, "pocket_reinforce_rank.csv"), index=False)

    off_cxc_hot = []
    off_txt_hot = []
    if not cxc.empty:
        off_cxc_hot = diff_cxc.head(N_OFFTARGET_PER_VARIANT)["resnum_model"].astype(int).tolist()
    if not txt.empty:
        off_txt_hot = diff_txt.head(N_OFFTARGET_PER_VARIANT)["resnum_model"].astype(int).tolist()

    rows=[]
    pocket_primary = POCKET_MODEL_PRIMARY
    pocket_full    = POCKET_MODEL

    for rn in pocket_primary:
        aa = wt.loc[wt["resnum_model"]==rn, "aa"].values[0] if rn in wt["resnum_model"].values else "X"
        for mut in POCKET_REINFORCE_PANEL:
            rows.append({"class":"POCKET_PRIMARY","resnum_model":rn,"paper":rn+68,"aa_native":aa,"mut":mut})

    for rn in (set(pocket_full) - set(pocket_primary)):
        aa = wt.loc[wt["resnum_model"]==rn, "aa"].values[0] if rn in wt["resnum_model"].values else "X"
        for mut in POCKET_REINFORCE_PANEL:
            rows.append({"class":"POCKET_FULL","resnum_model":rn,"paper":rn+68,"aa_native":aa,"mut":mut})

    for rn in off_cxc_hot:
        aa = cxc.loc[cxc["resnum_model"]==rn, "aa"].values[0] if rn in cxc["resnum_model"].values else "X"
        for mut in OFFTARGET_DISRUPT_PANEL:
            rows.append({"class":"OFF_CxC","resnum_model":rn,"paper":rn+68,"aa_native":aa,"mut":mut})

    for rn in off_txt_hot:
        aa = txt.loc[txt["resnum_model"]==rn, "aa"].values[0] if rn in txt["resnum_model"].values else "X"
        for mut in OFFTARGET_DISRUPT_PANEL:
            rows.append({"class":"OFF_TxT","resnum_model":rn,"paper":rn+68,"aa_native":aa,"mut":mut})

    design_df = pd.DataFrame(rows).drop_duplicates(subset=["class","resnum_model","mut"])
    design_df = design_df.sort_values(["class","resnum_model","mut"])
    design_csv = os.path.join(OUT_DIR, "design_candidates.csv")
    design_df.to_csv(design_csv, index=False)
    print("Wrote design candidates:", design_csv)

    try:
        name, pep, rec = read_baseline_complex(BASE_FASTA)
        out_fa = os.path.join(OUT_DIR, "AF3_variants_singletons.fasta")
        with open(out_fa, "w") as f:
            f.write(f">{name}__WT\n{pep}:{rec}\n")
            for _, r in design_df.iterrows():
                rn = int(r["resnum_model"])
                mut = r["mut"]
                rec2 = mutate(rec, rn, mut)
                tag = f"{r['class']}_{rn}{mut}"
                f.write(f">{name}__{tag}\n{pep}:{rec2}\n")
        print("Wrote FASTA:", out_fa)
    except Exception as e:
        print(f"[WARN] FASTA generation skipped: {e}")

    summary_txt = os.path.join(OUT_DIR, "next_steps_summary.txt")
    with open(summary_txt, "w") as f:
        f.write("# Pocket reinforcement targets (model idx):\n")
        f.write("primary: " + ",".join(map(str, POCKET_MODEL_PRIMARY)) + "\n")
        f.write("full:    " + ",".join(map(str, POCKET_MODEL)) + "\n\n")
        if off_cxc_hot:
            f.write("# Off-target CxC hotspots (model idx):\n")
            f.write(",".join(map(str, off_cxc_hot)) + "\n\n")
        if off_txt_hot:
            f.write("# Off-target TxT hotspots (model idx):\n")
            f.write(",".join(map(str, off_txt_hot)) + "\n\n")
        f.write("Design table: design_candidates.csv\n")
        f.write("Singleton FASTA: AF3_variants_singletons.fasta\n")
    print("Wrote:", summary_txt)


if __name__ == "__main__":
    main()