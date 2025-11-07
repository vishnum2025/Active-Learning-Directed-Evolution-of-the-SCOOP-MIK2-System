# #!/usr/bin/env python3
# # check_paper_pockets_one_model.py
# # Usage:
# #   python af3_pockets.py --jobdir fold_af3_mik2_scoop12 --s3 sd03.txt --model 0
# # (omit --model to auto-pick best by ipTM/rankingScore)

# import os, re, glob, json, argparse, math
# import numpy as np
# from Bio.PDB.MMCIFParser import MMCIFParser
# from Bio.PDB.Polypeptide import is_aa, PPBuilder
# from Bio.SeqUtils import seq1
# from Bio import pairwise2

# PAPER_POS = [246,268,292,294,316]
# HBOND_DIST_MAX = 3.5  # Å (relaxed)

# def read_full_mik2_from_S3(s3_path):
#     seq, want = [], False
#     for line in open(s3_path):
#         line=line.rstrip("\n")
#         if line.startswith(">"):
#             want = (line.strip()==">A_thaliana_AT4G08850")
#             continue
#         if want: seq.append(line.strip())
#     if not seq: raise SystemExit("Arabidopsis MIK2 not found in S3.")
#     return "".join(seq)

# def residues(chain): return [r for r in chain.get_residues() if is_aa(r, standard=True)]

# def pick_best_model(jobdir, force_model):
#     cifs = sorted(glob.glob(os.path.join(jobdir,"*_model_*.cif")))
#     if not cifs: raise SystemExit(f"No *_model_*.cif in {jobdir}")
#     if force_model is not None:
#         for c in cifs:
#             m=re.search(r"_model_(\d+)\.cif$",c)
#             if m and int(m.group(1))==force_model: return c, force_model
#         raise SystemExit("Requested model not found.")
#     scored=[]
#     for c in cifs:
#         m=re.search(r"_model_(\d+)\.cif$",c); mid=int(m.group(1)) if m else 0
#         s=c.replace("_model_","_summary_confidences_").replace(".cif",".json")
#         iptm=rank=-1.0
#         if os.path.exists(s):
#             try:
#                 J=json.load(open(s)); iptm=J.get("ipTM",J.get("iptm",-1.0)); rank=J.get("rankingScore",-1.0)
#             except: pass
#         scored.append((iptm,rank,-mid,c))
#     scored.sort(reverse=True); _,_,negmid,cif=scored[0]
#     return cif,-negmid

# def longest_chain(struct):
#     chains=[]
#     for model in struct:
#         for ch in model: chains.append((ch,len(residues(ch))))
#         break
#     chains.sort(key=lambda t:t[1]); return chains[-1][0]

# def chain_seq(chain):
#     ppb=PPBuilder()
#     return "".join(str(pp.get_sequence()) for pp in ppb.build_peptides(chain))

# def map_full_to_model(full_seq, model_seq, full_positions):
#     aln = pairwise2.align.globalxx(full_seq, model_seq, one_alignment_only=True)
#     if not aln: return {p:None for p in full_positions}
#     aF,aM,*_=aln[0]; f_i=m_i=0; m={}
#     targ=set(full_positions)
#     for cF,cM in zip(aF,aM):
#         if cF!="-": f_i+=1
#         if cM!="-": m_i+=1
#         if cF!="- " and cM!="-":  # (space in quotes prevents accidental emdash transform)
#             pass
#         if cF!="-" and cM!="-":
#             if f_i in targ and f_i not in m: m[f_i]=m_i
#             if len(m)==len(full_positions): break
#     return {p:m.get(p,None) for p in full_positions}

# def atom_by_name(res, names):
#     if res is None: return None
#     for a in res.get_atoms():
#         if a.element!="H" and a.get_name() in names: return a
#     return None

# def main():
#     ap=argparse.ArgumentParser()
#     ap.add_argument("--jobdir", required=True)
#     ap.add_argument("--s3", required=True)
#     ap.add_argument("--model", type=int, default=None)
#     args=ap.parse_args()

#     full_mik2 = read_full_mik2_from_S3(args.s3)
#     parser = MMCIFParser(QUIET=True)
#     cif, mid = pick_best_model(args.jobdir, args.model)
#     struct = parser.get_structure("x", cif)

#     # detect receptor (longest) and peptide (shortest)
#     allchains = [ch for model in struct for ch in model]
#     allchains.sort(key=lambda c: len(residues(c)))
#     P = allchains[0]             # peptide (shortest)
#     R = allchains[-1]            # receptor (longest)

#     R_res = residues(R); P_res = residues(P)
#     R_seq = chain_seq(R); P_seq = chain_seq(P)

#     # peptide positions 5 and 7 must be Ser/Thr to test paper pockets
#     s5 = P_res[4] if len(P_res)>=5 and P_seq[4] in "ST" else None
#     s7 = P_res[6] if len(P_res)>=7 and P_seq[6] in "ST" else None

#     pos_map = map_full_to_model(full_mik2, R_seq, PAPER_POS)
#     targets = {}
#     for p in PAPER_POS:
#         mi = pos_map[p]
#         targets[p] = R_res[mi-1] if mi and 1<=mi<=len(R_res) else None

#     # print the identity at mapped positions
#     print(f"# {os.path.basename(args.jobdir)} model {mid}  R={R.id} P={P.id}")
#     print("# PaperPos  ModelIdx  ChainResnum  AA")
#     for p in PAPER_POS:
#         res = targets[p]
#         if res is None:
#             print(f"{p:9d} {str(None):>9} {str('-'):>12}  ?")
#         else:
#             try: aa = seq1(res.get_resname())
#             except: aa="X"
#             print(f"{p:9d} {pos_map[p]:9d} {res.get_id()[1]:12d}  {aa}")

#     # quick distance checks for “could-be” H-bonds (no angle), S/T OG/OG1 to D/N/His atoms
#     ser_like = ["OG","OG1"]; aspO=["OD1","OD2"]; asnOD1=["OD1"]; hisN=["ND1","NE2"]

#     def min_dist(resA, namesA, resB, namesB):
#         if resA is None or resB is None: return math.inf
#         A = [atom_by_name(resA,[n]) for n in namesA]; A=[a for a in A if a]
#         B = [atom_by_name(resB,[n]) for n in namesB]; B=[b for b in B if b]
#         if not A or not B: return math.inf
#         return float(min(np.linalg.norm(a.get_coord()-b.get_coord()) for a in A for b in B))

#     print("\n# Distances (Å) — smaller is better; ≲3.5 can be H-bond-like")
#     def show(label, d):
#         s = f"{d:4.2f}" if d<math.inf else "NA"
#         print(f"{label:12s}: {s}")

#     d_s5_d246 = min_dist(s5,ser_like, targets[246],aspO)
#     d_s5_n268 = min_dist(s5,ser_like, targets[268],asnOD1)
#     d_s7_s292 = min_dist(s7,ser_like, targets[292],ser_like)
#     d_s7_h294 = min_dist(s7,ser_like, targets[294],hisN)
#     d_s7_h316 = min_dist(s7,ser_like, targets[316],hisN)

#     show("S5–D246", d_s5_d246)
#     show("S5–N268", d_s5_n268)
#     show("S7–S292", d_s7_s292)
#     show("S7–H294", d_s7_h294)
#     show("S7–H316", d_s7_h316)

# if __name__ == "__main__":
#     main()



# #!/usr/bin/env python3
# """
# Quick check: does the AF3 receptor chain sequence align with full-length MIK2
# starting at paper residue 24? Shows offset and amino acid comparisons.

# Hardwired for your setup:
#   - job_dir: fold_af3_mik2_scoop12
#   - full_s3 : pnas.2400862121.sd03.txt
# """

# import os, re
# import numpy as np

# # ---------------- PATHS (EDIT IF MOVED) ----------------
# JOB_DIR = "/Users/vishnu/Desktop/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/AF3_MIK2_outputs/fold_af3_mik2_scoop12"
# FULL_S3 = "/Users/vishnu/Desktop/v_files/GT/sem1/Research_sem1/sd03.txt"
# PAPER_POSITIONS = [246, 268, 292, 294, 316]
# # -------------------------------------------------------

# AA3to1 = {
#     "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
#     "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
#     "THR":"T","TRP":"W","TYR":"Y","VAL":"V","SEC":"U","PYL":"O"
# }

# def read_full_mik2_from_S3(txt_path):
#     seq = []
#     with open(txt_path) as f:
#         cap = False
#         for line in f:
#             line = line.rstrip("\n")
#             if line.startswith(">"):
#                 cap = (line.strip() == ">A_thaliana_AT4G08850")
#                 continue
#             if cap:
#                 if line.startswith(">"):
#                     break
#                 seq.append(line.strip())
#     return "".join(seq).replace(" ", "")

# def pick_one_cif(job_dir):
#     cifs = [p for p in os.listdir(job_dir) if p.endswith(".cif")]
#     if not cifs:
#         raise RuntimeError(f"No .cif files in {job_dir}")
#     cifs.sort()
#     pref = [p for p in cifs if re.search(r"_model_0\.cif$", p)]
#     return os.path.join(job_dir, (pref[0] if pref else cifs[0]))

# def extract_receptor_chain_sequence_from_cif(cif_path):
#     cols = {}
#     rows = []
#     with open(cif_path) as f:
#         lines = f.readlines()

#     in_loop = False
#     cur_head = []
#     for i, line in enumerate(lines):
#         s = line.strip()
#         if s.startswith("loop_"):
#             in_loop = True
#             cur_head = []
#             continue
#         if in_loop and s.startswith("_"):
#             cur_head.append(s)
#             continue
#         if in_loop and cur_head and not s.startswith("_"):
#             if any(h.startswith("_atom_site.") for h in cur_head):
#                 for h in cur_head:
#                     if h in ["_atom_site.label_asym_id","_atom_site.label_seq_id","_atom_site.label_comp_id","_atom_site.type_symbol"]:
#                         cols[h] = cur_head.index(h)
#                 j = i
#                 while j < len(lines):
#                     t = lines[j].rstrip("\n")
#                     if t.strip().startswith("_") or t.strip().startswith("loop_") or t.strip() == "":
#                         break
#                     parts = t.split()
#                     try:
#                         asym = parts[cols["_atom_site.label_asym_id"]]
#                         seqid = parts[cols["_atom_site.label_seq_id"]]
#                         comp  = parts[cols["_atom_site.label_comp_id"]]
#                         elem  = parts[cols["_atom_site.type_symbol"]]
#                         rows.append((asym, seqid, comp, elem))
#                     except Exception:
#                         pass
#                     j += 1
#                 break

#     chain_res = {}
#     for asym, seqid, comp, elem in rows:
#         if elem == "H": continue
#         try:
#             rid = int(seqid)
#         except ValueError:
#             continue
#         chain_res.setdefault(asym, {})
#         chain_res[asym][rid] = comp

#     chain_seq = {}
#     for asym, resmap in chain_res.items():
#         order = sorted(resmap)
#         seq = [AA3to1.get(resmap[rid], "X") for rid in order]
#         chain_seq[asym] = ("".join(seq), order)

#     rec = max(chain_seq.items(), key=lambda kv: len(kv[1][0]))
#     rec_chain, (rec_seq, rec_order) = rec
#     return rec_chain, rec_seq, rec_order

# def align_local(a, b):
#     # crude local alignment using substring search
#     best = (-1, -1)
#     k = 12
#     for i in range(len(a)-k+1):
#         seed = a[i:i+k]
#         j = b.find(seed)
#         if j != -1:
#             return j, j+len(a), len(seed)
#     return -1, -1, 0

# def main():
#     cif = pick_one_cif(JOB_DIR)
#     rec_chain, rec_seq, rec_order = extract_receptor_chain_sequence_from_cif(cif)
#     full_seq = read_full_mik2_from_S3(FULL_S3)

#     start_b, end_b, score = align_local(rec_seq, full_seq)
#     offset = (start_b + 1)

#     print(f"Model: {os.path.basename(cif)}  receptor chain: {rec_chain}")
#     print(f"Receptor chain length: {len(rec_seq)}  Full MIK2 length: {len(full_seq)}")
#     print(f"Inferred offset = {offset}\n")

#     print("First 10 residues (model_idx → paper_pos : modelAA / paperAA):")
#     for i in range(10):
#         model_idx = rec_order[i]
#         paper_pos = offset + i
#         print(f"{model_idx:>4} → {paper_pos:>4} : {rec_seq[i]} / {full_seq[paper_pos-1]}")

#     print("\nPaper positions check:")
#     for pp in PAPER_POSITIONS:
#         delta = pp - offset
#         if 0 <= delta < len(rec_seq):
#             model_idx = rec_order[delta]
#             model_aa  = rec_seq[delta]
#             paper_aa  = full_seq[pp-1]
#             print(f"{pp:>4} → {model_idx:>4} : {paper_aa} / {model_aa}")
#         else:
#             print(f"{pp:>4} → (out of range)")

# if __name__ == "__main__":
#     main()

#!/usr/bin/env python3
"""
Check amino acids at key positions in the AF3 receptor chain
without applying any offset. Compare directly with paper residues.
"""

from Bio.PDB import MMCIFParser
from Bio.Data.IUPACData import protein_letters_3to1

# --- CONFIG ---
CIF_PATH = "/Users/vishnu/Desktop/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/AF3_MIK2_outputs/fold_af3_mik2_scoop12/fold_af3_mik2_scoop12_model_0.cif"
CHAIN_ID = "B"   # receptor chain

# Paper residues to check (full-length numbering + reported amino acid)
paper_residues = {
    246: "D",
    268: "N",
    292: "S",
    294: "H",
    316: "H"
}

# Convert 3-letter to 1-letter AA safely
def three_to_one_safe(resname):
    resname = resname.capitalize()
    if resname in protein_letters_3to1:
        return protein_letters_3to1[resname]
    return "X"  # unknown / non-standard residue

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("mik2", CIF_PATH)
chain = structure[0][CHAIN_ID]

print(f"Receptor chain {CHAIN_ID}, length = {len([r for r in chain])}\n")

print("Model index → amino acid (no offset) vs Paper residue (number, aa)")
for paper_pos, paper_aa in paper_residues.items():
    try:
        res = chain[(" ", paper_pos, " ")]
        aa_model = three_to_one_safe(res.get_resname())
        print(f"Model {paper_pos}: {aa_model}   Paper {paper_pos}: {paper_aa}")
    except KeyError:
        print(f"Model {paper_pos}: (not found in model chain)   Paper {paper_pos}: {paper_aa}")