#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# from itertools import combinations, product

# # --------- HARD-CODED INPUT/OUTPUTS ---------
# IN_FASTA = "/Users/vishnu/Desktop/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/CxC_TxT/AF3_Cviol1838-SCOOP12.fasta"
# OUT_PRIMARY = "/Users/vishnu/Desktop/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/CxC_TxT/1838_electroneg_primary.fasta"
# # Model positions for 1838 MIK2 (mapped from paper): 246→202, 268→224, 292→248, 294→250, 316→272
# POS_PRIMARY = [224, 248]   # mutate only these
# ACIDS = ["D", "E"]         # electronegative choices
# # --------------------------------------------

# def read_single_record_fasta(path):
#     name, seq = None, []
#     with open(path) as f:
#         for line in f:
#             line = line.strip()
#             if not line: 
#                 continue
#             if line.startswith(">"):
#                 if name is not None:
#                     raise SystemExit("FASTA has multiple records; expected one 'PEPTIDE:RECEPTOR'")
#                 name = line[1:].strip()
#             else:
#                 seq.append(line)
#     if name is None:
#         raise SystemExit("No FASTA header found")
#     full = "".join(seq).upper()
#     if ":" not in full:
#         raise SystemExit("Expected 'PEPTIDE:RECEPTOR' in sequence")
#     pep, rec = full.split(":", 1)
#     return name, pep, rec

# def mutate_rec(rec: str, model_pos: int, to_aa: str) -> str:
#     i = model_pos - 1  # model numbering is 1-based
#     if not (0 <= i < len(rec)):
#         raise ValueError(f"Model pos {model_pos} out of range for receptor length {len(rec)}")
#     return rec[:i] + to_aa + rec[i+1:]

# def wrap80(s: str):
#     return "\n".join(s[i:i+80] for i in range(0, len(s), 80))

# def write_fasta(records, out_path):
#     with open(out_path, "w") as f:
#         for name, seq in records:
#             f.write(f">{name}\n{wrap80(seq)}\n")

# def build_primary(name, pep, rec):
#     out = []
#     out.append((f"{name}__WTrec", f"{pep}:{rec}"))  # WT passthrough
#     # singles
#     for p in POS_PRIMARY:
#         for aa in ACIDS:
#             r2 = mutate_rec(rec, p, aa)
#             out.append((f"{name}__{p}{aa}", f"{pep}:{r2}"))
#     # doubles
#     for (p1, p2) in combinations(POS_PRIMARY, 2):
#         for aa1, aa2 in product(ACIDS, ACIDS):
#             r2 = mutate_rec(rec, p1, aa1)
#             r3 = mutate_rec(r2,  p2, aa2)
#             out.append((f"{name}__{p1}{aa1}_{p2}{aa2}", f"{pep}:{r3}"))
#     return out

# def main():
#     name, pep, rec = read_single_record_fasta(IN_FASTA)
#     records = build_primary(name, pep, rec)
#     write_fasta(records, OUT_PRIMARY)
#     print(f"Wrote {len(records)} variants → {OUT_PRIMARY}")
#     print("Positions (model):", POS_PRIMARY, "Acids:", ACIDS)

# if __name__ == "__main__":
#     main()

####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Map AtMIK2 pocket residues (D246/N268/S292/H294/H316) onto the 1838 MIK2
# (Cleome violacea) by ECD alignment, and also report approximate LRR repeat numbers.

# - Input FASTA is hard-coded (sd03.txt).
# - AtMIK2 header and 1838 header patterns are hard-coded.
# - Uses Biopython pairwise2 + BLOSUM62 if present; else a simple global aligner.
# """

# import os, sys, re
# from collections import namedtuple

# INPUT_FASTA = "/Users/vishnu/Desktop/v_files/GT/sem1/Research_sem1/sd03.txt"  # put this script in the same folder as sd03.txt, or adjust the path
# ATMIK2_ID   = "A_thaliana_AT4G08850"           # AtMIK2
# MIK1838_ID_PATTERN = "C_violacea_Clevi0007s1838___MIK2 OR LIKE"              # the 1838 homolog to map onto (Cleome violacea)

# PAPER_SITES = [("D",246), ("N",268), ("S",292), ("H",294), ("H",316)]  # AtMIK2 full-length numbering

# # ----------------- FASTA I/O -----------------
# def read_fasta(path):
#     seqs = {}
#     hdr = None
#     chunks = []
#     with open(path, "r") as f:
#         for line in f:
#             if line.startswith(">"):
#                 if hdr:
#                     seqs[hdr] = "".join(chunks).replace(" ", "").replace("\n","")
#                 hdr = line[1:].strip()
#                 chunks = []
#             else:
#                 chunks.append(line.strip())
#         if hdr:
#             seqs[hdr] = "".join(chunks).replace(" ", "").replace("\n","")
#     return seqs

# def pick_headers(seqs, at_id, target_pat):
#     at_hdr = None
#     tg_hdr = None
#     for h in seqs:
#         if at_id in h:
#             at_hdr = h
#         if target_pat in h and tg_hdr is None:
#             tg_hdr = h
#     if not at_hdr:
#         raise SystemExit(f"Could not find AtMIK2 header containing '{at_id}' in {INPUT_FASTA}")
#     if not tg_hdr:
#         raise SystemExit(f"Could not find 1838 header containing '{target_pat}' in {INPUT_FASTA}")
#     return at_hdr, tg_hdr

# # ----------------- ECD finder (heuristic TM) -----------------
# HYDRO = set("AILMFWYV")  # hydrophobic
# def find_tm_start(seq, win=19, frac=0.7):
#     best_i, best_f = None, 0.0
#     for i in range(0, len(seq)-win+1):
#         w = seq[i:i+win]
#         f = sum(aa in HYDRO for aa in w) / float(win)
#         if f > best_f:
#             best_f, best_i = f, i
#     # Choose first window with high hydrophobicity as TM start if clearly hydrophobic
#     if best_f >= frac:
#         return best_i
#     # fallback: try motif "RRK" a bit upstream of kinase; if present, TM likely just before it; use the most hydrophobic 19-mer before RRK
#     m = re.search(r"R{2,}K", seq)
#     if m:
#         jmax = max(0, m.start()-win-5)
#         best_i2, best_f2 = None, 0.0
#         for i in range(0, jmax+1):
#             w = seq[i:i+win]
#             f = sum(aa in HYDRO for aa in w) / float(win)
#             if f > best_f2:
#                 best_f2, best_i2 = f, i
#         if best_i2 is not None:
#             return best_i2
#     # last resort: near the largest hydrophobic window
#     return best_i or len(seq)-1

# def slice_ecd(seq):
#     tm = find_tm_start(seq)
#     return seq[:tm], tm  # ECD and TM index (0-based, in full-length)

# # ----------------- Alignment -----------------
# have_biopy = False
# try:
#     from Bio import pairwise2
#     from Bio.SubsMat import MatrixInfo as mat
#     have_biopy = True
# except Exception:
#     have_biopy = False

# def simple_global_align(a, b, match=1, mismatch=-1, gap_open=-2, gap_extend=-1):
#     """
#     Minimal Needleman–Wunsch with affine-like (very simplified) gap penalty:
#     gap_open for starting a gap, gap_extend for continuing.
#     Returns aligned_a, aligned_b
#     """
#     # DP matrices
#     n, m = len(a), len(b)
#     # scores
#     M = [[0]*(m+1) for _ in range(n+1)]
#     X = [[-10**9]*(m+1) for _ in range(n+1)]  # gaps in a
#     Y = [[-10**9]*(m+1) for _ in range(n+1)]  # gaps in b
#     # traceback
#     TB = [[None]*(m+1) for _ in range(n+1)]

#     M[0][0] = 0
#     for i in range(1, n+1):
#         X[i][0] = gap_open + (i-1)*gap_extend
#         M[i][0] = X[i][0]
#         TB[i][0] = 'U'  # up (gap in b)
#     for j in range(1, m+1):
#         Y[0][j] = gap_open + (j-1)*gap_extend
#         M[0][j] = Y[0][j]
#         TB[0][j] = 'L'  # left (gap in a)

#     def sc(x,y): return match if x==y else mismatch

#     for i in range(1, n+1):
#         ai = a[i-1]
#         for j in range(1, m+1):
#             bj = b[j-1]
#             # extend/ open gaps
#             X[i][j] = max(M[i-1][j] + gap_open, X[i-1][j] + gap_extend)
#             Y[i][j] = max(M[i][j-1] + gap_open, Y[i][j-1] + gap_extend)
#             # match/mismatch
#             MM = M[i-1][j-1] + sc(ai,bj)
#             # choose best
#             best = MM
#             tb = 'D'
#             if X[i][j] > best:
#                 best, tb = X[i][j], 'U'
#             if Y[i][j] > best:
#                 best, tb = Y[i][j], 'L'
#             M[i][j] = best
#             TB[i][j] = tb

#     # traceback
#     i, j = n, m
#     al_a, al_b = [], []
#     while i>0 or j>0:
#         tb = TB[i][j]
#         if tb == 'D':
#             al_a.append(a[i-1]); al_b.append(b[j-1]); i-=1; j-=1
#         elif tb == 'U':
#             al_a.append(a[i-1]); al_b.append('-'); i-=1
#         else:
#             al_a.append('-'); al_b.append(b[j-1]); j-=1
#     return "".join(reversed(al_a)), "".join(reversed(al_b))

# def align_seqs(a, b):
#     if have_biopy:
#         # BLOSUM62 scoring
#         matrix = mat.blosum62
#         # pairwise2 returns multiple alignments; take best first
#         alns = pairwise2.align.globalds(a, b, matrix, -10, -0.5)  # gap open/extend
#         aa, bb, score, *_ = alns[0]
#         return aa, bb
#     else:
#         return simple_global_align(a, b)

# # ----------------- LRR indexing -----------------
# LRRmotif = re.compile(r"L..L.L..", re.I)
# def lrr_repeat_indices(seq):
#     """
#     Very rough LRR finder: returns a list of (start,end) for regions matching LxxLxLxx,
#     merges nearby hits (<=6 aa apart) into one repeat. Then assigns a "repeat number"
#     for each sequence position: pos_to_repeat[pos] = 1..N (or 0 if none).
#     """
#     hits = [m.span() for m in LRRmotif.finditer(seq)]
#     if not hits:
#         return [0]*len(seq), []
#     # merge close hits
#     merged = []
#     cur_s, cur_e = hits[0]
#     for s,e in hits[1:]:
#         if s - cur_e <= 6:
#             cur_e = max(cur_e, e)
#         else:
#             merged.append((cur_s, cur_e))
#             cur_s, cur_e = s, e
#     merged.append((cur_s, cur_e))
#     # assign numbers by order
#     pos_to_rep = [0]*len(seq)
#     for idx,(s,e) in enumerate(merged, start=1):
#         for k in range(s, min(e, len(seq))):
#             pos_to_rep[k] = idx
#     return pos_to_rep, merged

# # ----------------- Mapping helpers -----------------
# def build_gapped_to_ungapped_map(gapped):
#     """Return list mapping gapped index -> ungapped index (1-based), or 0 if gap."""
#     m = []
#     c = 0
#     for ch in gapped:
#         if ch == '-':
#             m.append(0)
#         else:
#             c += 1
#             m.append(c)
#     return m

# def ungapped_to_gapped_index(gapped, ungapped_pos1):
#     """
#     Given gapped seq string and a 1-based position in its ungapped form,
#     return the 0-based index in the gapped string where that residue sits.
#     """
#     c = 0
#     for i,ch in enumerate(gapped):
#         if ch != '-':
#             c += 1
#             if c == ungapped_pos1:
#                 return i
#     return None

# # ----------------- Main -----------------
# def main():
#     print(f"Reading: {INPUT_FASTA}")
#     seqs = read_fasta(INPUT_FASTA)
#     at_hdr, tg_hdr = pick_headers(seqs, ATMIK2_ID, MIK1838_ID_PATTERN)
#     at_full = seqs[at_hdr]
#     tg_full = seqs[tg_hdr]
#     print(f"Found AtMIK2: {at_hdr} (len={len(at_full)})")
#     print(f"Found 1838  : {tg_hdr} (len={len(tg_full)})")

#     at_ecd, at_tm_idx = slice_ecd(at_full)
#     tg_ecd, tg_tm_idx = slice_ecd(tg_full)
#     print(f"AtMIK2 ECD length ≈ {len(at_ecd)} (TM starts ~{at_tm_idx})")
#     print(f"1838   ECD length ≈ {len(tg_ecd)} (TM starts ~{tg_tm_idx})")

#     # Align ECDs
#     aA, aB = align_seqs(at_ecd, tg_ecd)
#     assert len(aA) == len(aB)
#     print(f"Aligned ECD length: {len(aA)}")

#     # LRR repeat numbering within ECDs
#     at_lrr_map, at_lrr_blocks = lrr_repeat_indices(at_ecd)
#     tg_lrr_map, tg_lrr_blocks = lrr_repeat_indices(tg_ecd)

#     # For each paper site (full-length At index), find its position in ECD and then map across alignment
#     Rec = namedtuple("Rec", "label aa_at at_idx_full at_ecd_1based align_col aa_1838 tg_idx_full tg_lrr at_lrr")
#     rows = []
#     for aa_lab, at_idx in PAPER_SITES:
#         if at_idx < 1 or at_idx > len(at_full):
#             print(f"[WARN] {aa_lab}{at_idx} outside AtMIK2 length")
#             continue
#         # position in ECD (1-based)
#         at_ecd_pos1 = at_idx if at_idx <= len(at_ecd) else None
#         if at_ecd_pos1 is None:
#             print(f"[WARN] {aa_lab}{at_idx} is past ECD (falls inside TM/after); skipping")
#             continue
#         # safety check amino acid label
#         aa_actual = at_full[at_idx-1]
#         if aa_actual.upper() != aa_lab:
#             print(f"[NOTE] AtMIK2 position {at_idx} is '{aa_actual}', not '{aa_lab}'. Using actual '{aa_actual}'.")
#             aa_lab = aa_actual.upper()

#         # gapped column index corresponding to this AtMIK2 ECD residue
#         col = ungapped_to_gapped_index(aA, at_ecd_pos1)
#         if col is None:
#             print(f"[WARN] could not map At {aa_lab}{at_idx} into alignment")
#             continue

#         # get aligned residue in 1838
#         aa_1838 = aB[col]
#         # map back to 1838 full-length index
#         # first: 1838 ECD 1-based position at this alignment column
#         tg_ecd_pos1 = ungapped_to_gapped_index(aB, 1)  # placeholder
#         # build a direct map once for speed
#         # (gapped_index -> ungapped_index)
#         # do it outside loop? small overhead anyway
#         g2u_at = build_gapped_to_ungapped_map(aA)
#         g2u_tg = build_gapped_to_ungapped_map(aB)
#         at_ecd_pos1_here = g2u_at[col]  # should equal at_ecd_pos1
#         tg_ecd_pos1_here = g2u_tg[col]  # 0 if gap

#         if tg_ecd_pos1_here == 0 or aa_1838 == '-':
#             # aligned to a gap in 1838
#             tg_full_idx = None
#         else:
#             tg_full_idx = tg_ecd_pos1_here  # ECD index 1-based
#             # convert to full-length 1-based index by adding 0 (ECD starts at 1)
#             # (because we sliced from N-terminus)
#             # If you ever offset ECD from mid-sequence, add that offset here.

#         # LRR repeat numbers (0 if none)
#         at_lrr_num = at_lrr_map[at_ecd_pos1-1] if 1 <= at_ecd_pos1 <= len(at_lrr_map) else 0
#         tg_lrr_num = tg_lrr_map[tg_ecd_pos1_here-1] if (tg_ecd_pos1_here and 1 <= tg_ecd_pos1_here <= len(tg_lrr_map)) else 0

#         rows.append(Rec(
#             label=f"{aa_lab}{at_idx}",
#             aa_at=aa_lab,
#             at_idx_full=at_idx,
#             at_ecd_1based=at_ecd_pos1,
#             align_col=col,
#             aa_1838=aa_1838 if aa_1838 != '-' else '(gap)',
#             tg_idx_full=(tg_ecd_pos1_here if tg_ecd_pos1_here else None),
#             tg_lrr=tg_lrr_num,
#             at_lrr=at_lrr_num
#         ))

#     print("\n=== Mapping (AtMIK2 → 1838 MIK2) ===")
#     print("Label  | At(full)  At(LRR) | 1838 AA  1838(full)  1838(LRR) | Note")
#     print("-------+--------------------+-------------------------------+------")
#     for r in rows:
#         note = ""
#         if r.aa_1838 == "(gap)":
#             note = "aligned-to-gap"
#         print(f"{r.label:<6} | {r.at_idx_full:>7}   R{r.at_lrr:>2} | {r.aa_1838:^7} {str(r.tg_idx_full or '-'):>11}   R{r.tg_lrr:>2} | {note}")

#     # Quick sanity: show how many LRR blocks we detected
#     print("\nApprox. LRR blocks detected:")
#     print(f"  AtMIK2: {len(at_lrr_blocks)} blocks  (crude motif-based)")
#     print(f"  1838  : {len(tg_lrr_blocks)} blocks  (crude motif-based)")

# if __name__ == "__main__":
#     main()


########################################################################################################################################################################



####################################################################################################################################################################################################################################################################################################################




# scripts/af3_pep_to_receptor_contacts.py
import numpy as np
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1  # <-- use this instead of three_to_one

def heavy_coords(res):
    xyz = [a.get_coord() for a in res.get_atoms() if getattr(a, "element", "") != "H"]
    return np.array(xyz) if xyz else None

def nearest_contacts(cif_path, cutoff=5.0):
    parser = MMCIFParser(QUIET=True)
    st = parser.get_structure("af3", cif_path)

    # receptor = longest protein chain; peptide = shortest
    chains = []
    for ch in st.get_chains():
        n_aa = sum(1 for r in ch.get_residues() if is_aa(r, standard=True))
        if n_aa > 0:
            chains.append((ch, n_aa))
    if len(chains) < 2:
        raise RuntimeError("Expected ≥2 protein chains in the CIF (receptor + peptide).")

    pep = min(chains, key=lambda t: t[1])[0]
    rec = max(chains, key=lambda t: t[1])[0]

    pep_res = [r for r in pep.get_residues() if is_aa(r, standard=True)]
    rec_res = [r for r in rec.get_residues() if is_aa(r, standard=True)]

    # Build receptor lookup (number, 1-letter, residue object)
    rec_table = []
    for r in rec_res:
        rn = r.get_id()[1]  # model numbering in the CIF
        try:
            aa = seq1(r.get_resname(), custom_map={"SEC": "C"})  # tolerant 3->1
        except Exception:
            aa = "X"
        rec_table.append((rn, aa, r))

    # For each peptide position (1..N), gather receptor hits within cutoff
    results = {}
    for ip, p in enumerate(pep_res, start=1):
        Pc = heavy_coords(p)
        if Pc is None:
            continue
        hits = []
        for rn, aa, r in rec_table:
            Rc = heavy_coords(r)
            if Rc is None:
                continue
            dmin = float(np.min(np.linalg.norm(Rc[:, None, :] - Pc[None, :, :], axis=2)))
            if dmin <= cutoff:
                hits.append((dmin, rn, aa))
        hits.sort(key=lambda x: x[0])  # by distance
        results[ip] = hits
    return results

if __name__ == "__main__":
    cif = "/Users/vishnu/Desktop/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/1838_SCOOPs/fold_1838_scoop12/fold_1838_scoop12_model_0.cif"
    res = nearest_contacts(cif, cutoff=5.0)
    print("1838 MIK2-SCOOP12 ( model 1)")
    for pep_pos, hits in res.items():
        if not hits:
            continue
        top = ", ".join([f"{aa}{rn}@{d:.1f}Å" for d, rn, aa in hits[:5]])
        print(f"Pep {pep_pos}: {top}")


##### +68 offset required for 1838 MIK2 to paper numbering 