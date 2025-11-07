#!/usr/bin/env python3
# scripts/make_contact_clusters.py
#
# Builds:
#   - pocket_cluster_WT_S5S7.csv (WT SCOOP12–MIK2 contacts at S5/S7)
#   - offtarget_cluster_CxC_S5S7.csv and offtarget_cluster_TxT_S5S7.csv
# Adds peptide position + amino acid information in all outputs.

import os, glob
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

# ========= PATHS =========
BASE = "/Users/vishnu/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline"

WT_JOB = os.path.join(BASE, "1838_SCOOPs", "fold_1838_scoop12")

VARIANT_JOB_DIRS = {
    "CxC": os.path.join(BASE, "CxC_TxT", "fold_cxc_1838"),
    "TxT": os.path.join(BASE, "CxC_TxT", "fold_txt_1838"),
}

OUTDIR = os.path.join(BASE, "analysis_contacts")
os.makedirs(OUTDIR, exist_ok=True)

# ========= SETTINGS =========
DIST_CUTOFF = 5.0
PAPER_OFFSET = 68
S_POSITIONS = [5, 7]
PEPTIDE_LEN_MIN = 5
PEPTIDE_LEN_MAX = 40


# ========= HELPERS =========
def heavy_coords(res):
    xyz = [a.get_coord() for a in res.get_atoms() if getattr(a, "element", "") != "H"]
    return np.array(xyz) if xyz else None


def residues(chain):
    return [r for r in chain.get_residues() if is_aa(r, standard=True)]


def pick_receptor_peptide(structure):
    chains = []
    for ch in structure.get_chains():
        n = sum(1 for r in ch.get_residues() if is_aa(r, standard=True))
        if n > 0:
            chains.append((ch, n))
    if len(chains) < 2:
        raise RuntimeError("Expected ≥2 protein chains")
    chains_sorted = sorted(chains, key=lambda t: t[1])
    pep_candidates = [ch for ch, n in chains_sorted if PEPTIDE_LEN_MIN <= n <= PEPTIDE_LEN_MAX]
    pep = pep_candidates[0] if pep_candidates else chains_sorted[0][0]
    rec = max(chains, key=lambda t: t[1])[0]
    return rec, pep


def contact_triplets(rec, pep, cutoff=DIST_CUTOFF, only_pep_positions=None):
    rec_res = residues(rec)
    pep_res = residues(pep)
    out = []
    for ip, p in enumerate(pep_res, start=1):
        if only_pep_positions and ip not in only_pep_positions:
            continue
        Pc = heavy_coords(p)
        if Pc is None:
            continue
        try:
            pep_aa = seq1(p.get_resname(), custom_map={"SEC": "C"})
        except Exception:
            pep_aa = "X"
        for r in rec_res:
            Rc = heavy_coords(r)
            if Rc is None:
                continue
            dmin = float(np.min(np.linalg.norm(Rc[:, None, :] - Pc[None, :, :], axis=2)))
            if dmin <= cutoff:
                rn = int(r.get_id()[1])
                try:
                    rec_aa = seq1(r.get_resname(), custom_map={"SEC": "C"})
                except Exception:
                    rec_aa = "X"
                out.append((rn, rec_aa, ip, pep_aa, dmin))
    return out


def scan_job(job_dir, only_pep_positions=None):
    parser = MMCIFParser(QUIET=True)
    cifs = sorted(glob.glob(os.path.join(job_dir, "*_model_*.cif")))
    all_contacts = []
    for cif in cifs:
        try:
            st = parser.get_structure("af3", cif)
            rec, pep = pick_receptor_peptide(st)
            trips = contact_triplets(rec, pep, DIST_CUTOFF, only_pep_positions)
            all_contacts.extend(trips)
        except Exception as e:
            print(f"[WARN] {os.path.basename(job_dir)}: {os.path.basename(cif)} -> {e}")
    return all_contacts


def tally_contacts(contacts):
    by_res = defaultdict(list)
    for rn, rec_aa, pep_pos, pep_aa, d in contacts:
        by_res[(rn, rec_aa)].append((pep_pos, pep_aa, d))
    rows = []
    for (rn, rec_aa), lst in by_res.items():
        n = len(lst)
        pep_pos_counts = Counter(pp for pp, _, _ in lst)
        pep_pos_str = ",".join(str(p) for p in sorted(pep_pos_counts.keys()))
        pep_aa_str = ",".join(sorted(set(a for _, a, _ in lst)))
        min_d = float(np.min([d for _, _, d in lst]))
        rows.append({
            "resnum_model": rn,
            "resnum_paper": rn + PAPER_OFFSET,
            "aa": rec_aa,
            "n_contacts": n,
            "pep_positions": ";".join(f"{p}:{c}" for p, c in sorted(pep_pos_counts.items())),
            "pep_positions_str": pep_pos_str,
            "pep_aas": pep_aa_str,
            "min_dmin": min_d
        })
    df = pd.DataFrame(rows).sort_values(["n_contacts", "min_dmin"], ascending=[False, True])
    return df


def main():
    print("\n=== WT pocket cluster (S5 & S7) ===")
    if os.path.isdir(WT_JOB):
        wt_contacts = scan_job(WT_JOB, only_pep_positions=S_POSITIONS)
        wt_df = tally_contacts(wt_contacts)
        wt_out = os.path.join(OUTDIR, "pocket_cluster_WT_S5S7.csv")
        wt_df.to_csv(wt_out, index=False)
        print(wt_df.head(25).to_string(index=False))
        print("Saved:", wt_out)
    else:
        print(f"[WARN] WT job dir not found: {WT_JOB}")

    # Variant clusters restricted to S5/S7
    for tag, jd in VARIANT_JOB_DIRS.items():
        if not os.path.isdir(jd):
            print(f"[WARN] Variant dir missing: {jd}")
            continue
        print(f"\n=== {tag} off-target cluster (S5 & S7) ===")
        contacts = scan_job(jd, only_pep_positions=S_POSITIONS)
        if not contacts:
            print(f"[WARN] No contacts found for {tag}")
            continue
        df = tally_contacts(contacts)
        out_csv = os.path.join(OUTDIR, f"offtarget_cluster_{tag}_S5S7.csv")
        df.to_csv(out_csv, index=False)
        print(df.head(25).to_string(index=False))
        print("Saved:", out_csv)


if __name__ == "__main__":
    main()