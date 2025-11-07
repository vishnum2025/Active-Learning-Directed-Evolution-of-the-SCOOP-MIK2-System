#!/usr/bin/env python3
import argparse, os, sys, json, glob, subprocess, pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def run_colabfold(fasta, outdir, num_models=5, num_recycles=3, model_type="alphafold2_multimer_v3"):
    os.makedirs(outdir, exist_ok=True)
    cmd = [
        "colabfold_batch",
        "--model-type", model_type,
        "--num-models", str(num_models),
        "--num-recycle", str(num_recycles),
        # Mirror paper: paired+unpaired default, no templates
        fasta, outdir
    ]
    print(">> Running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)

def find_result_jsons(run_dir):
    # ColabFold writes result_model_*.json in the run folder
    return glob.glob(os.path.join(run_dir, "result_model_*.json"))

def load_metrics(result_json_path):
    with open(result_json_path, "r") as f:
        data = json.load(f)
    # iptm/pTM keys may vary slightly; use get with fallback
    iptm = data.get("iptm", data.get("ipTM", None))
    ptm  = data.get("ptm", data.get("pTM", None))
    pae  = data.get("pae", None)  # a 2D list if present
    plddt = data.get("plddt", None)
    ranking_score = None
    if iptm is not None and ptm is not None:
        ranking_score = 0.8*iptm + 0.2*ptm
    # Also store path to corresponding PDB if present
    base = os.path.basename(result_json_path).replace("result_", "").replace(".json","")
    pdb_guess = os.path.join(os.path.dirname(result_json_path), f"{base}.pdb")
    if not os.path.exists(pdb_guess):
        # fall back to ranked_*.pdb
        cand = glob.glob(os.path.join(os.path.dirname(result_json_path), "ranked_*.pdb"))
        pdb_guess = cand[0] if cand else None
    return {
        "json": result_json_path,
        "pdb": pdb_guess,
        "iptm": iptm,
        "ptm": ptm,
        "ranking_score": ranking_score,
        "pae": pae,
        "plddt": plddt
    }

def pae_heatmap(pae, png_path, title=None):
    if pae is None:
        return
    arr = np.array(pae)
    plt.figure(figsize=(5,4))
    im = plt.imshow(arr, origin="lower")
    plt.colorbar(im, fraction=0.046, pad=0.04, label="PAE (Å)")
    if title:
        plt.title(title)
    plt.tight_layout()
    plt.savefig(png_path, dpi=200)
    plt.close()

def process_single_fasta(fasta_path, outdir, num_models, num_recycles, model_type, skip_run=False):
    base = os.path.splitext(os.path.basename(fasta_path))[0]
    run_dir = os.path.join(outdir, base)
    if not skip_run:
        run_colabfold(fasta_path, run_dir, num_models, num_recycles, model_type)
    # Parse metrics
    rows = []
    best = None
    for js in find_result_jsons(run_dir):
        met = load_metrics(js)
        rows.append(met)
        if met["ranking_score"] is not None:
            if best is None or met["ranking_score"] > best["ranking_score"]:
                best = met
    # Save per-run summary
    df = pd.DataFrame(rows)
    df["fasta"] = fasta_path
    df["run_dir"] = run_dir
    df.to_csv(os.path.join(run_dir, "afm_results_summary.csv"), index=False)
    # Save overall top info + PAE plot
    top_csv = os.path.join(run_dir, "top_model.csv")
    if best is not None:
        pd.DataFrame([best]).to_csv(top_csv, index=False)
        if best["pae"] is not None:
            png = os.path.join(run_dir, "top_model_PAE.png")
            pae_heatmap(best["pae"], png, title=f"{base} | ipTM={best['iptm']:.3f}, pTM={best['ptm']:.3f}")
    return rows

def main():
    p = argparse.ArgumentParser(description="Run AF-Multimer via ColabFold and aggregate metrics (ipTM/pTM/PAE).")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--fasta", type=str, help="Path to a single FASTA (chains separated by ':')")
    g.add_argument("--fastadir", type=str, help="Directory of FASTA files to batch")
    p.add_argument("--outdir", type=str, required=True, help="Output directory root")
    p.add_argument("--num-models", type=int, default=5)
    p.add_argument("--num-recycles", type=int, default=3)
    p.add_argument("--model-type", type=str, default="alphafold2_multimer_v3")
    p.add_argument("--skip-run", action="store_true", help="Skip calling colabfold_batch; just parse existing outputs")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    all_rows = []
    if args.fasta:
        all_rows.extend(process_single_fasta(args.fasta, args.outdir, args.num_models, args.num_recycles, args.model_type, args.skip_run))
    else:
        fastas = [f for f in glob.glob(os.path.join(args.fastadir, "*")) if f.lower().endswith((".fa",".fasta",".faa"))]
        for fp in sorted(fastas):
            try:
                all_rows.extend(process_single_fasta(fp, args.outdir, args.num_models, args.num_recycles, args.model_type, args.skip_run))
            except subprocess.CalledProcessError as e:
                print(f"[WARN] colabfold_batch failed for {fp}: {e}", file=sys.stderr)

    # Global summary
    if all_rows:
        df = pd.DataFrame(all_rows)
        df.to_csv(os.path.join(args.outdir, "summary.csv"), index=False)
        print(f"Done. Wrote global summary to {os.path.join(args.outdir, 'summary.csv')}")
    else:
        print("No results found or produced.")

if __name__ == "__main__":
    main()
