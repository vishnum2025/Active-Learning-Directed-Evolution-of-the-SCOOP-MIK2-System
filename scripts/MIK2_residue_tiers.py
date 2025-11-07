#!/usr/bin/env python3
"""
Aggregate per-residue interface stats from AF3 server outputs and assign Tier A/B/C.

Fallback strategy:
- Prefer TRUE per-token PAE from *_full_data_k.json (keys: full_pae, token_chain_ids).
- Else fallback to chain-pair PAE from *_summary_confidences_k.json (key: chain_pair_pae_min off-diagonal).
- Else compute tiers without PAE (will push more residues out of Tier A/B).

Hard-wired paths for your project. Adjust if you moved folders.

"""

import os, re, glob, math, json
from collections import defaultdict
import numpy as np
import pandas as pd
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

# ---------------- CONFIG ----------------
ROOT = "/Users/vishnu/Desktop/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/CxC_TxT/cxc_results/mik2_mutations"
PATTERN = "fold_cxc_*"
CONSENSUS = "/Users/vishnu/Desktop/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/CxC_TxT/cxc_results/mik2_mutations/consensus_interface_receptor_af3.csv"
OUTDIR = "/Users/vishnu/Desktop/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/CxC_TxT/cxc_results/mik2_mutations/2_tier_results"
PAPER_TO_MODEL = {246: 202, 268: 224, 292: 248, 294: 250, 316: 272}
MODEL_TO_PAPER = {v: k for k, v in PAPER_TO_MODEL.items()}
POCKET_PAPER = [246, 268, 292, 294, 316]
pocket_model   = [PAPER_TO_MODEL[p] for p in POCKET_PAPER if p in PAPER_TO_MODEL]
DIST_CUTOFF = 5.0
# ----------------------------------------

def residues(chain):
    return [r for r in chain.get_residues() if is_aa(r, standard=True)]

def atom_coords(res):
    xyz = [a.get_coord() for a in res.get_atoms() if a.element != "H"]
    return np.array(xyz) if xyz else None

def res_plddt(res):
    vals=[a.get_bfactor() for a in res.get_atoms() if a.element != "H"]
    return float(np.mean(vals)) if vals else float("nan")

def pick_receptor_peptide(structure):
    """Longest chain = receptor; shortest = peptide."""
    chains=[]
    for model in structure:
        for ch in model:
            chains.append((ch, len(residues(ch))))
        break
    chains.sort(key=lambda t: t[1])
    P = chains[0][0]
    R = chains[-1][0]
    return R, P

def contact_pairs(R, P, cutoff=DIST_CUTOFF):
    Rres, Pres = residues(R), residues(P)
    pairs=[]
    for Rr in Rres:
        Cr = atom_coords(Rr)
        if Cr is None: continue
        for Pr in Pres:
            Cp = atom_coords(Pr)
            if Cp is None: continue
            dmin = float(np.min(np.linalg.norm(Cr[:,None,:]-Cp[None,:,:], axis=2)))
            if dmin <= cutoff:
                pairs.append((Rr.get_id()[1], Pr.get_id()[1], dmin))
    return pairs

def find_models(job_dir):
    return sorted(glob.glob(os.path.join(job_dir, "*_model_*.cif")))

def base_from_cif(cif_path):
    """strip trailing _model_k.cif -> base prefix used by AF3 server."""
    base = os.path.basename(cif_path)
    return re.sub(r"_model_\d+\.cif$", "", base)

def find_jsons_for_model(job_dir, cif_path):
    """Return (full_data_json_or_None, summary_json_or_None, k)."""
    m = re.search(r"_model_(\d+)\.cif$", os.path.basename(cif_path))
    k = m.group(1) if m else "0"
    base = base_from_cif(cif_path)
    fd = os.path.join(job_dir, f"{base}_full_data_{k}.json")
    sj = os.path.join(job_dir, f"{base}_summary_confidences_{k}.json")
    return (fd if os.path.exists(fd) else None,
            sj if os.path.exists(sj) else None,
            int(k))

def load_full_token_info(full_json):
    """Return (chain_to_token_indices: dict[label->np.array], full_pae: np.ndarray)."""
    J = json.load(open(full_json))
    if "token_chain_ids" not in J or "full_pae" not in J:
        raise KeyError("missing token_chain_ids or full_pae")
    tc = np.array(J["token_chain_ids"], dtype=object)  # allow 'A','B' or ints
    chain_to_idx = {}
    for lab in pd.unique(tc):
        chain_to_idx[lab] = np.where(tc == lab)[0]
    full_pae = np.array(J["full_pae"], dtype=float)
    return chain_to_idx, full_pae

def chain_labels_by_length(structure, chain_to_token):
    counts = {lab: len(idx) for lab, idx in chain_to_token.items()}
    pep_lab = min(counts, key=lambda k: counts[k])
    rec_lab = max(counts, key=lambda k: counts[k])
    return rec_lab, pep_lab

def token_index_maps(R_chain, P_chain, rec_token_idx, pep_token_idx):
    Rres = residues(R_chain); Pres = residues(P_chain)
    nR = min(len(Rres), len(rec_token_idx))
    nP = min(len(Pres), len(pep_token_idx))
    rec_tok2res = { int(rec_token_idx[i]): Rres[i].get_id()[1] for i in range(nR) }
    pep_tok2res = { int(pep_token_idx[i]): Pres[i].get_id()[1] for i in range(nP) }
    rec_res2tok = { v:k for k,v in rec_tok2res.items() }
    pep_res2tok = { v:k for k,v in pep_tok2res.items() }
    return rec_res2tok, pep_res2tok

def read_chain_pair_pae(summary_json):
    """Return a scalar chain-pair PAE (Å) if present; else None."""
    try:
        S = json.load(open(summary_json))
        # Expect a 2x2 matrix for binary complex; use off-diagonal min
        m = S.get("chain_pair_pae_min", None)
        if m is None:
            return None
        M = np.array(m, dtype=float)
        if M.shape[0] >= 2 and M.shape[1] >= 2:
            # Take the smaller off-diagonal (i,j) or (j,i)
            v = float(min(M[0,1], M[1,0]))
            return v
        return None
    except Exception:
        return None

def main():
    os.makedirs(OUTDIR, exist_ok=True)

    # load consensus frequencies (from your earlier AF3 summary)
    cons = pd.read_csv(CONSENSUS)
    freq_map = dict(zip(cons["resnum"].astype(int), cons["frequency"].astype(float)))

    pocket_model = [PAPER_TO_MODEL[p] for p in POCKET_PAPER if p in PAPER_TO_MODEL]

    parser = MMCIFParser(QUIET=True)
    jobs = sorted([d for d in glob.glob(os.path.join(ROOT, PATTERN)) if os.path.isdir(d)])
    if not jobs:
        raise SystemExit(f"No job dirs found under {ROOT} matching '{PATTERN}'")

    jobs_with_contact = defaultdict(set)
    iface_pae_vals    = defaultdict(list)   # TRUE or fallback chain-pair PAE per receptor residue
    local_plddt_vals  = defaultdict(list)
    min_dist_vals     = defaultdict(list)

    n_models = 0
    n_with_full_pae = 0
    n_with_chainpair = 0

    for job in jobs:
        design = os.path.basename(job.rstrip("/"))
        cifs = find_models(job)
        if not cifs:
            print(f"[WARN] {design}: no *_model_*.cif files")
            continue

        per_model_min_dist = defaultdict(list)
        rn_had_contact = set()

        for cif in cifs:
            n_models += 1
            structure = parser.get_structure(design, cif)
            try:
                R, P = pick_receptor_peptide(structure)
            except Exception as e:
                print(f"[WARN] {design}: cannot pick chains for {os.path.basename(cif)} ({e})")
                continue

            # pLDDT per receptor residue
            Rres = residues(R)
            for r in Rres:
                rn = r.get_id()[1]
                local_plddt_vals[rn].append(res_plddt(r))

            # distance to (paper) pocket
            ca = {r.get_id()[1]: r["CA"].get_coord() for r in Rres if "CA" in r}
            pocket_coords = [ca.get(pm) for pm in (pocket_model) if pm in ca]
            pocket_coords = [x for x in pocket_coords if x is not None]
            for rn, coord in ca.items():
                if pocket_coords:
                    d = min(np.linalg.norm(coord - pc) for pc in pocket_coords)
                    per_model_min_dist[rn].append(d)

            # contacts
            pairs = contact_pairs(R, P)
            for rn, _pn, _d in pairs:
                rn_had_contact.add(int(rn))

            # try full-data JSON
            full_json, summary_json, mid = find_jsons_for_model(job, cif)
            used_pae = False
            if full_json:
                try:
                    chain_to_token, full_pae = load_full_token_info(full_json)
                    rec_lab, pep_lab = chain_labels_by_length(structure, chain_to_token)
                    rec_idx = chain_to_token[rec_lab]
                    pep_idx = chain_to_token[pep_lab]
                    rec_res2tok, pep_res2tok = token_index_maps(R, P, rec_idx, pep_idx)

                    for rn, pn, _d in pairs:
                        tR = rec_res2tok.get(int(rn))
                        tP = pep_res2tok.get(int(pn))
                        if tR is None or tP is None:
                            continue
                        pae = float(full_pae[int(tR), int(tP)])
                        iface_pae_vals[int(rn)].append(pae)
                        used_pae = True
                except Exception as e:
                    # fall through to chain-pair fallback
                    pass

            if (not used_pae) and summary_json:
                # chain-pair PAE fallback: assign same value to all contacting residues for this model
                cp = read_chain_pair_pae(summary_json)
                if cp is not None and len(pairs) > 0:
                    for rn, _pn, _d in pairs:
                        iface_pae_vals[int(rn)].append(float(cp))
                    n_with_chainpair += 1

            if used_pae:
                n_with_full_pae += 1

        for rn in rn_had_contact:
            jobs_with_contact[int(rn)].add(design)
        for rn, dlist in per_model_min_dist.items():
            if dlist:
                min_dist_vals[rn].append(float(np.mean(dlist)))

    # Build table
    all_res = set(freq_map.keys()) | set(local_plddt_vals.keys()) | set(min_dist_vals.keys()) | set(jobs_with_contact.keys())
    rows=[]
    n_jobs = len(jobs)
    for rn in sorted(all_res):
        freq = float(freq_map.get(rn, 0.0))
        jobs_pct = 100.0 * len(jobs_with_contact.get(rn, set())) / n_jobs if n_jobs > 0 else 0.0
        mean_pae = float(np.mean(iface_pae_vals[rn])) if iface_pae_vals.get(rn) else math.nan
        mean_plddt = float(np.mean(local_plddt_vals[rn])) if local_plddt_vals.get(rn) else math.nan
        min_d = float(np.min(min_dist_vals[rn])) if min_dist_vals.get(rn) else math.nan
        rows.append({
            "resnum_model": rn,
            "resnum_paper": MODEL_TO_PAPER.get(rn, np.nan),
            "frequency": freq,
            "jobs_with_contact_pct": jobs_pct,
            "mean_interface_PAE": mean_pae,     # Å (true or chain-pair fallback)
            "mean_local_pLDDT": mean_plddt,
            "min_distance_to_pocket": min_d
        })

    df = pd.DataFrame(rows)

    # Tiers (using real PAE if available, else chain-pair PAE which tends to be ~5–15 Å for good interfaces)
    def assign_tier(r):
        f   = r["frequency"]
        job = r["jobs_with_contact_pct"]
        pae = r["mean_interface_PAE"]
        pl  = r["mean_local_pLDDT"]
        d   = r["min_distance_to_pocket"]
        # thresholds: tight if per-token PAE available; still reasonable with chain-pair PAE
        good_conf = (not np.isnan(pae) and pae <= 5.0) and (not np.isnan(pl) and pl >= 70.0)
        good_jobs = job >= 20.0
        near      = (not np.isnan(d)) and d <= 8.0
        if f >= 0.01 and good_jobs and good_conf and near: return "Tier A"
        if f >= 0.01 and good_jobs and good_conf:         return "Tier B"
        if f >= 0.01 and ((8.0 < pae <= 12.0) or (60.0 <= pl < 65.0)): return "Tier C"
        return "Other"

    df["Tier"] = df.apply(assign_tier, axis=1)

    os.makedirs(OUTDIR, exist_ok=True)
    out_csv = os.path.join(OUTDIR, "MIK2_residue_tiers.csv")
    df.to_csv(out_csv, index=False)

    print(f"Wrote {out_csv}")
    print(f"Jobs={n_jobs}, Models={n_models}, with per-token PAE={n_with_full_pae}, with chain-pair fallback={n_with_chainpair}")
    #print(f"Pocket residues (paper): {POCKET_PAPER}, ECD start={ECD_START}")

    # show paper residues for quick sanity
    paper = df[df["resnum_paper"].isin(POCKET_PAPER)][
        ["resnum_paper","frequency","jobs_with_contact_pct","mean_interface_PAE","mean_local_pLDDT","min_distance_to_pocket","Tier"]
    ].sort_values("resnum_paper")
    print("\nPaper-pocket residues summary:")
    if len(paper):
        print(paper.to_string(index=False))
    else:
        print("(no rows matched)")
    
if __name__=="__main__":
    main()

