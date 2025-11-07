#!/usr/bin/env python3


# to run:
#
# python /Users/vishnu/Desktop/GT/sem1/Research/AF3_tools/AF3_summary.py \
#   --root /Users/vishnu/Desktop/GT/sem1/Research/AF3_jobs/1838_SCOOPs \
#   --pattern "fold_1838_scoop*" \
#   --outdir /Users/vishnu/Desktop/GT/sem1/Research/AF3_jobs/1838_results

## OR

# ./scripts/af3_summary_results.sh --root ./alde_athaliana --pattern "fold_athaliana_*" --outdir ./alde_athaliana





import os, re, glob, json, math, argparse
from collections import defaultdict
import numpy as np
import pandas as pd
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

DIST_CUTOFF = 5.0  # Å, heavy-atom contact cutoff

def residues(chain):
    return [r for r in chain.get_residues() if is_aa(r, standard=True)]

def atom_coords(res):
    xyz=[]
    for a in res.get_atoms():
        if a.element != "H":
            xyz.append(a.get_coord())
    return np.array(xyz) if xyz else None

def res_plddt(res):
    # AF3 stores pLDDT in B-factor (0..100)
    vals=[a.get_bfactor() for a in res.get_atoms() if a.element != "H"]
    return float(np.mean(vals)) if vals else float("nan")

def pick_receptor_peptide(struct):
    """Return (R_chain, P_chain) using longest=R (MIK2 ECD), shortest=P (SCOOP)."""
    chains = list(struct.get_chains())
    if len(chains) < 2:
        raise RuntimeError("Expected ≥2 chains (peptide + receptor)")
    lens = [len(residues(c)) for c in chains]
    R = chains[int(np.argmax(lens))]
    P = chains[int(np.argmin(lens))]
    return R, P

def contacts_and_iface_plddt(R, P):
    Rres, Pres = residues(R), residues(P)
    contacts=[]
    for ir, Rr in enumerate(Rres):
        Cr = atom_coords(Rr)
        if Cr is None: 
            continue
        for ip, Pr in enumerate(Pres):
            Cp = atom_coords(Pr)
            if Cp is None: 
                continue
            dmin = np.min(np.linalg.norm(Cr[:,None,:]-Cp[None,:,:], axis=2))
            if dmin <= DIST_CUTOFF:
                w = 0.5*(res_plddt(Rr)+res_plddt(Pr))/100.0  # 0–1 weight
                contacts.append({
                    "resiR": Rr.get_id()[1], "resiP": Pr.get_id()[1],
                    "dmin": float(dmin), "w": float(w)
                })
    iface_plddt = float(np.mean([c["w"]*100.0 for c in contacts])) if contacts else float("nan")
    return contacts, iface_plddt, len(Rres), len(Pres), R.id, P.id

def load_summary_conf(summary_json):
    """AF3 server 'summary_confidences_k.json' keys can vary; return what exists."""
    d = {}
    if os.path.exists(summary_json):
        try:
            d = json.load(open(summary_json))
        except Exception:
            d = {}
    # Keep only a permissive subset of known/confidence-ish fields
    keep = {}
    for k in ["iptm","ipTM","ptm","pTM","rankingScore","mean_plddt",
              "interface_confidence","interface_plddt",
              "overall_confidence","overall_plddt"]:
        if k in d:
            keep[k] = d[k]
    return keep

def summarize_job(job_dir):
    """Process one AF3 job folder; return (rows, rec_counter, pep_counter)."""
    parser = MMCIFParser(QUIET=True)
    rows=[]
    rec_counter=defaultdict(float)
    pep_counter=defaultdict(float)

    cif_list = sorted(glob.glob(os.path.join(job_dir, "*_model_*.cif")))
    if not cif_list:
        return rows, rec_counter, pep_counter  # empty job

    # Infer design name from folder or file
    design = os.path.basename(job_dir.rstrip("/"))

    for cif in cif_list:
        m = re.search(r"_model_(\d+)\.cif$", cif)
        mid = int(m.group(1)) if m else 0
        struct = parser.get_structure(f"{design}_m{mid}", cif)
        R, P = pick_receptor_peptide(struct)
        cont, iface_plddt, nR, nP, RcID, PcID = contacts_and_iface_plddt(R, P)

        # Tally consensus across models (global, not per-design here)
        for c in cont:
            rec_counter[c["resiR"]] += c["w"]
            pep_counter[c["resiP"]] += c["w"]

        # Load companion summary JSON if present
        summ_json = cif.replace("_model_", "_summary_confidences_").replace(".cif",".json")
        summ = load_summary_conf(summ_json)

        row = {
            "design": design,
            "model": mid,
            "mmcif": os.path.abspath(cif),
            "n_contacts": len(cont),
            "interface_pLDDT_mean": iface_plddt,
            "receptor_chain_id": RcID, "peptide_chain_id": PcID,
            "n_receptor_len": nR, "n_peptide_len": nP,
        }
        row.update(summ)
        rows.append(row)

    return rows, rec_counter, pep_counter

def to_df(counter):
    tot = sum(counter.values()) or 1.0
    items = sorted(counter.items(), key=lambda kv:(-kv[1], kv[0]))
    return pd.DataFrame([{"resnum":k, "weighted_count":v, "frequency":v/tot} for k,v in items])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="AF3_outputs", help="Folder with AF3 job subfolders")
    ap.add_argument("--pattern", default="fold_af3_mik2_scoop*", help="Glob for job dirs inside root")
    ap.add_argument("--outdir", default="results_af3", help="Output directory for CSVs")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    job_dirs = sorted([d for d in glob.glob(os.path.join(args.root, args.pattern)) if os.path.isdir(d)])
    if not job_dirs:
        raise SystemExit(f"No job dirs found under {args.root} matching '{args.pattern}'")

    all_rows=[]
    rec_counter=defaultdict(float)
    pep_counter=defaultdict(float)

    for jd in job_dirs:
        rows, rc, pc = summarize_job(jd)
        all_rows.extend(rows)
        for k,v in rc.items(): rec_counter[k] += v
        for k,v in pc.items(): pep_counter[k] += v

    if not all_rows:
        raise SystemExit("No models found across jobs.")

    df = pd.DataFrame(all_rows)

    # Choose a sorting key for "best model": prefer explicit rankingScore/ipTM if present; else interface pLDDT, then contacts
    sort_cols = [c for c in ["rankingScore","ipTM","iptm","interface_pLDDT_mean","n_contacts"] if c in df.columns]
    sort_asc  = [False]*len(sort_cols)
    df_sorted = df.sort_values(sort_cols, ascending=sort_asc)

    # Per-design top pick
    top = df_sorted.groupby("design", as_index=False).first()

    # Write outputs
    all_path = os.path.join(args.outdir, "AF3_summary_all.csv")
    top_path = os.path.join(args.outdir, "AF3_summary_top.csv")
    rec_path = os.path.join(args.outdir, "consensus_interface_receptor_af3.csv")
    pep_path = os.path.join(args.outdir, "consensus_interface_peptide_af3.csv")

    df_sorted.to_csv(all_path, index=False)
    top.to_csv(top_path, index=False)
    to_df(rec_counter).to_csv(rec_path, index=False)
    to_df(pep_counter).to_csv(pep_path, index=False)

    print(f"Wrote:\n  {all_path}\n  {top_path}\n  {rec_path}\n  {pep_path}")

if __name__ == "__main__":
    main()