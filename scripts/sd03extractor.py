#!/usr/bin/env python3
"""
Build AF3 input FASTAs for three non-Arabidopsis receptors:
  - C. violacea Clevi0007s1839
  - C. violacea Clevi0007s1838
  - T. cacao Thecc08G107800 (outgroup)

Each output = SCOOP peptide + receptor ECD:
  AF3_<receptorShort>-SCOOP<n>.fasta
    >design_AF3_<receptorShort>-SCOOP<n>
    <PEPTIDE>:<RECEPTOR_ECD>

Usage
-----
python sd03extractor.py \
  --s1 sd01.xlsx \
  --s3 sd03.txt \
  --outdir af3_three \
  --scoops 2,5,7,10,11,12,17,21,22,23,24,27,31,34,41,44,46,47,50
"""
import argparse, re
from pathlib import Path
import pandas as pd

# --- robust ECD trim (same logic you used) ------------------------------------
def find_first_long_hydrophobic_run(seq: str, window=21, min_frac=0.75, start_region_frac=0.6):
    hydro = set("VILMFWYAGC")
    for i in range(len(seq) - window):
        frag = seq[i:i+window]
        if sum(aa in hydro for aa in frag) >= int(min_frac * window):
            if i > len(seq) * start_region_frac:
                return i
    return None

def trim_ecd(seq: str, start=24, end=None):
    if end is None:
        tm = find_first_long_hydrophobic_run(seq, window=21, min_frac=0.75, start_region_frac=0.6)
        end = tm if tm is not None else min(len(seq), 710)
    return seq[start:end], (start, end)

# --- S1 reader: get SCOOP peptide (13-mer preferred) --------------------------
def load_scoop_sequence_from_S1(xlsx_path: str, scoop_label: str) -> str | None:
    xl = pd.ExcelFile(xlsx_path)
    df = xl.parse(xl.sheet_names[0])

    # any row mentioning the label
    mask = df.apply(lambda col: col.astype(str).str.contains(scoop_label, case=False, na=False))
    rows = df[mask.any(axis=1)]
    if rows.empty:
        return None

    # prefer A. thaliana entries if 'species' exists
    if "species" in rows.columns:
        ath = rows[rows["species"].astype(str).str.contains("Athaliana", na=False)]
        if not ath.empty:
            rows = ath

    # prefer the 13-mer column if present
    col_13 = [c for c in rows.columns if c.strip().lower().startswith("13mer alignment")]
    if col_13:
        val = str(rows.iloc[0][col_13[0]]).strip()
        if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{8,30}", val):
            return val

    # otherwise scan all cells for an AA-only string
    for c in rows.columns:
        vals = rows[c].dropna().astype(str).tolist()
        for v in vals:
            if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{8,30}", v):
                return v
    return None

# --- S3 reader: extract receptor sequences by header tokens --------------------
def extract_by_header_tokens(s3_path: str, required_tokens: list[str]) -> str | None:
    """
    Find a header line starting with '>' containing all required_tokens (exact case),
    then return the sequence (concatenated lines) until the next header.
    """
    lines = Path(s3_path).read_text().splitlines()
    hit_idx = None
    for i, line in enumerate(lines):
        if line.startswith(">") and all(tok in line for tok in required_tokens):
            hit_idx = i
            break
    if hit_idx is None:
        return None
    seq = []
    for j in range(hit_idx+1, len(lines)):
        if lines[j].startswith(">"):
            break
        seq.append(lines[j].strip())
    return "".join(seq)

RECEPTORS = {
    # tokens to match in the header  -> short name used in output filenames
    tuple(["C_violacea", "Clevi0007s1839"]): "Cviol1839",
    tuple(["C_violacea", "Clevi0007s1838"]): "Cviol1838",
    tuple(["T_cacao",    "Thecc08G107800"]): "Tcacao",
}

def write_joint_fasta(peptide: str, receptor_ecd: str, out_path: Path):
    stem = out_path.stem
    rec = f">design_{stem}\n{peptide}:{receptor_ecd}\n"
    out_path.write_text(rec)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--s1", required=True, help="Dataset S1 Excel (SCOOPs)")
    ap.add_argument("--s3", required=True, help="Dataset S3 FASTA/TXT (homologues)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--scoops", default="2,5,7,10,11,12,17,21,22,23,24,41,44,46,47,50",
                    help="Comma-separated SCOOP numbers to build")
    ap.add_argument("--ecd-start", type=int, default=24, help="ECD start (default 24)")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    scoops = sorted({int(x) for x in args.scoops.split(",") if x.strip()})

    # Load receptors
    rec_seqs = {}
    for tokens, short in RECEPTORS.items():
        full = extract_by_header_tokens(args.s3, list(tokens))
        if not full:
            print(f"[WARN] Could not find receptor with tokens {tokens} in {args.s3}; skipping.")
            continue
        ecd, (s,e) = trim_ecd(full, start=args.ecd_start, end=None)
        print(f"{short}: ECD len={len(ecd)} (start={s}, end≈{e})")
        rec_seqs[short] = ecd

    made = []
    for n in scoops:
        label = f"SCOOP{n}"
        pep = load_scoop_sequence_from_S1(args.s1, label)
        if pep is None:
            print(f"[WARN] No sequence for {label} in S1; skip.")
            continue
        for short, ecd in rec_seqs.items():
            out_fa = outdir / f"AF3_{short}-SCOOP{n}.fasta"
            write_joint_fasta(pep, ecd, out_fa)
            made.append(out_fa.name)

    print(f"Wrote {len(made)} FASTAs to {outdir}")
    for nm in made[:6]:
        print("  ", nm)
    if len(made) > 6:
        print("  ...")

if __name__ == "__main__":
    main()