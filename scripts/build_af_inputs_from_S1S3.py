#!/usr/bin/env python3
"""
Build AF-Multimer input FASTAs for MIK2–SCOOP (±BAK1) directly from the authors' datasets:

- Dataset S1: pnas.2400862121.sd01.xlsx  (SCOOP peptides)
- Dataset S3: pnas.2400862121.sd03.txt   (MIK2 homologues; includes Arabidopsis AT4G08850)

Outputs FASTA files named like the authors' folders, e.g.:
  AF2_MIK2-SCOOP12.fasta                (binary: SCOOP:MIK2_ECD)
  AF3_MIK2-SCOOP12.fasta                (binary: SCOOP:MIK2_ECD)
  AF3_MIK2-SCOOP12_BAK1.fasta           (ternary: SCOOP:MIK2_ECD:BAK1_ECD)

Usage examples
--------------
1) Minimal (uses heuristic ECD trimming, no BAK1):
   python build_af_inputs_from_S1S3.py --s1 pnas.2400862121.sd01.xlsx --s3 pnas.2400862121.sd03.txt --outdir af_inputs

2) Specify which SCOOPs to generate:
   python build_af_inputs_from_S1S3.py --s1 sd01.xlsx --s3 sd03.txt --scoops 2,5,7,12 --outdir af_inputs

3) Provide exact MIK2 ECD boundaries (recommended if you know them, e.g., 24–709):
   python build_af_inputs_from_S1S3.py --s1 sd01.xlsx --s3 sd03.txt --mik2-ecd-start 24 --mik2-ecd-end 709 --outdir af_inputs

4) Include BAK1 and auto-trim its ECD from a FASTA you provide:
   python build_af_inputs_from_S1S3.py --s1 sd01.xlsx --s3 sd03.txt --bak1-fasta BAK1_Q94F62.fasta --outdir af_inputs

5) Zip the outputs:
   python build_af_inputs_from_S1S3.py --s1 sd01.xlsx --s3 sd03.txt --outdir af_inputs --zip
"""
import argparse, re, os
from pathlib import Path
import pandas as pd

def read_text(p: str) -> str:
    return Path(p).read_text()

def extract_fasta(text: str, header_exact: str) -> str:
    """Extract sequence for an exact FASTA header line (without trailing spaces)."""
    seq_lines, cap = [], False
    for line in text.splitlines():
        if line.startswith(">"):
            cap = (line.strip() == header_exact)
            continue
        if cap:
            if line.startswith(">"):  # shouldn't happen inside capture
                break
            seq_lines.append(line.strip())
    return "".join(seq_lines)

def find_first_long_hydrophobic_run(seq: str, window=21, min_frac=0.75, start_region_frac=0.6):
    hydro = set("VILMFWYAGC")
    for i in range(len(seq)-window):
        frag = seq[i:i+window]
        if sum(aa in hydro for aa in frag) >= int(min_frac*window):
            if i > len(seq)*start_region_frac:
                return i
    return None

def trim_ecd(seq: str, start=None, end=None) -> tuple[str, tuple[int,int]]:
    """Trim ectodomain. If start/end not given, use 24 and first long C-term hydrophobic run as TM start."""
    if start is None:
        start = 24
    if end is None:
        tm = find_first_long_hydrophobic_run(seq, window=21, min_frac=0.75, start_region_frac=0.6)
        end = tm if tm is not None else min(len(seq), 710)
    return seq[start:end], (start, end)

def load_scoop_sequence_from_S1(xlsx_path: str, scoop_label: str) -> str | None:
    xl = pd.ExcelFile(xlsx_path)
    df = xl.parse(xl.sheet_names[0])
    # rows containing the label in any column
    mask = df.apply(lambda col: col.astype(str).str.contains(scoop_label, case=False, na=False))
    rows = df[mask.any(axis=1)]
    # prefer A. thaliana entries if column exists
    if "species" in rows.columns:
        ath = rows[rows["species"].astype(str).str.contains("Athaliana", na=False)]
        if not ath.empty:
            rows = ath
    if rows.empty:
        return None
    # The sheet we saw had peptide sequences in 'Unnamed: 6'; fallback to AA-only regex if schema differs
    if "13mer alignment SCOOP " in rows.columns and rows["13mer alignment SCOOP "].dropna().shape[0] > 0:
        return str(rows.iloc[0]["13mer alignment SCOOP "]).strip()
    for c in rows.columns:
        vals = rows[c].dropna().astype(str).tolist()
        for v in vals:
            if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{8,30}", v):
                return v
    return None

def load_mik2_full_from_S3(txt_path: str) -> str:
    text = read_text(txt_path)
    return extract_fasta(text, ">A_thaliana_AT4G08850")

def load_bak1_from_fasta(fasta_path: str) -> str:
    text = read_text(fasta_path)
    # Accept either a specific header or just take the first sequence in file
    seqs, cur = [], []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
                cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    if not seqs:
        raise ValueError("No sequence found in BAK1 FASTA")
    return seqs[0]

def write_joint_fasta(peptide: str, mik2_ecd: str, out_path: str, bak1_ecd: str | None = None):
    stem = Path(out_path).stem
    if bak1_ecd:
        rec = f">design_{stem}\n{peptide}:{mik2_ecd}:{bak1_ecd}\n"
    else:
        rec = f">design_{stem}\n{peptide}:{mik2_ecd}\n"
    Path(out_path).write_text(rec)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--s1", required=True, help="Path to Dataset S1 Excel (pnas.2400862121.sd01.xlsx)")
    ap.add_argument("--s3", required=True, help="Path to Dataset S3 FASTA/TXT (pnas.2400862121.sd03.txt)")
    ap.add_argument("--outdir", required=True, help="Output folder for FASTAs")
    ap.add_argument("--scoops", default="", help="Comma-separated list of SCOOP numbers (e.g., 2,5,7,12). If empty, use authors' AF2/AF3 union.")
    ap.add_argument("--mik2-ecd-start", type=int, default=None, help="Override MIK2 ECD start (e.g., 24)")
    ap.add_argument("--mik2-ecd-end", type=int, default=None, help="Override MIK2 ECD end (e.g., 709)")
    ap.add_argument("--bak1-fasta", default=None, help="Optional FASTA file for BAK1; if provided, ternary AF3_*_BAK1 FASTAs will be created for selected SCOOPs")
    ap.add_argument("--zip", action="store_true", help="Zip the output folder at the end")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    # Load receptor sequences
    mik2_full = load_mik2_full_from_S3(args.s3)
    if not mik2_full:
        raise RuntimeError("Could not find Arabidopsis MIK2 in S3 (expected header '>A_thaliana_AT4G08850').")
    mik2_ecd, (s, e) = trim_ecd(mik2_full, start=args.mik2_ecd_start, end=args.mik2_ecd_end)
    print(f"MIK2 ECD length = {len(mik2_ecd)} (start={s}, end={e})")

    bak1_ecd = None
    if args.bak1_fasta:
        bak1_full = load_bak1_from_fasta(args.bak1_fasta)
        bak1_ecd, (bs, be) = trim_ecd(bak1_full, start=None, end=None)
        print(f"BAK1 ECD length = {len(bak1_ecd)} (auto-trim; override with --mik2-ecd-start/end for MIK2 only)")

    # Which SCOOPs to build
    if args.scoops.strip():
        scoops = sorted({int(x.strip()) for x in args.scoops.split(",") if x.strip()})
        af3_with_bak1 = set()  # if user custom list, we won't guess which have BAK1
    else:
        # Authors' union we saw in their folders
        af2 = [2,5,7,12,17,21,22,23,24,41,46,47,10,11,44,50]  # included a few seen under AF2_* too
        af3 = [2,5,7,10,11,12,17,21,22,23,24,27,31,34,41,44,46,47,50]
        scoops = sorted(set(af2+af3))
        af3_with_bak1 = {5,11,12,21,27,34,46,50}

    made = []
    for n in scoops:
        label = f"SCOOP{n}"
        pep = load_scoop_sequence_from_S1(args.s1, label)
        if pep is None:
            print(f"[WARN] Could not find sequence for {label} in S1; skipping.")
            continue
        # AF2-style
        out_fa = outdir / f"AF2_MIK2-SCOOP{n}.fasta"
        write_joint_fasta(pep, mik2_ecd, out_fa, bak1_ecd if (bak1_ecd and n in af3_with_bak1) else None)
        made.append(str(out_fa))
        # AF3-style
        out_fa3 = outdir / f"AF3_MIK2-SCOOP{n}.fasta"
        write_joint_fasta(pep, mik2_ecd, out_fa3, None)
        made.append(str(out_fa3))
        # AF3 ternary with BAK1 where applicable and available
        if bak1_ecd and n in af3_with_bak1:
            out_fa3b = outdir / f"AF3_MIK2-SCOOP{n}_BAK1.fasta"
            write_joint_fasta(pep, mik2_ecd, out_fa3b, bak1_ecd)
            made.append(str(out_fa3b))

    print(f"Wrote {len(made)} FASTAs to {outdir}")
    if args.zip:
        zpath = outdir.with_suffix(".zip")
        import zipfile
        with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
            for p in made:
                z.write(p, arcname=Path(p).name)
        print(f"Zipped to {zpath}")

if __name__ == "__main__":
    main()
