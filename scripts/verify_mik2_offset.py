
import os, re
import numpy as np

# ---------------- PATHS (EDIT IF MOVED) ----------------
JOB_DIR = "/Users/vishnu/Desktop/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/AF3_MIK2_outputs/fold_af3_mik2_scoop12"
FULL_S3 = "/Users/vishnu/Desktop/v_files/GT/sem1/Research_sem1/sd03.txt"
PAPER_POSITIONS = [246, 268, 292, 294, 316]
# -------------------------------------------------------

AA3to1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V","SEC":"U","PYL":"O"
}

def read_full_mik2_from_S3(txt_path):
    seq = []
    with open(txt_path) as f:
        cap = False
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                cap = (line.strip() == ">A_thaliana_AT4G08850")
                continue
            if cap:
                if line.startswith(">"):
                    break
                seq.append(line.strip())
    return "".join(seq).replace(" ", "")

def pick_one_cif(job_dir):
    cifs = [p for p in os.listdir(job_dir) if p.endswith(".cif")]
    if not cifs:
        raise RuntimeError(f"No .cif files in {job_dir}")
    cifs.sort()
    pref = [p for p in cifs if re.search(r"_model_0\.cif$", p)]
    return os.path.join(job_dir, (pref[0] if pref else cifs[0]))

def extract_receptor_chain_sequence_from_cif(cif_path):
    cols = {}
    rows = []
    with open(cif_path) as f:
        lines = f.readlines()

    in_loop = False
    cur_head = []
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith("loop_"):
            in_loop = True
            cur_head = []
            continue
        if in_loop and s.startswith("_"):
            cur_head.append(s)
            continue
        if in_loop and cur_head and not s.startswith("_"):
            if any(h.startswith("_atom_site.") for h in cur_head):
                for h in cur_head:
                    if h in ["_atom_site.label_asym_id","_atom_site.label_seq_id","_atom_site.label_comp_id","_atom_site.type_symbol"]:
                        cols[h] = cur_head.index(h)
                j = i
                while j < len(lines):
                    t = lines[j].rstrip("\n")
                    if t.strip().startswith("_") or t.strip().startswith("loop_") or t.strip() == "":
                        break
                    parts = t.split()
                    try:
                        asym = parts[cols["_atom_site.label_asym_id"]]
                        seqid = parts[cols["_atom_site.label_seq_id"]]
                        comp  = parts[cols["_atom_site.label_comp_id"]]
                        elem  = parts[cols["_atom_site.type_symbol"]]
                        rows.append((asym, seqid, comp, elem))
                    except Exception:
                        pass
                    j += 1
                break

    chain_res = {}
    for asym, seqid, comp, elem in rows:
        if elem == "H": continue
        try:
            rid = int(seqid)
        except ValueError:
            continue
        chain_res.setdefault(asym, {})
        chain_res[asym][rid] = comp

    chain_seq = {}
    for asym, resmap in chain_res.items():
        order = sorted(resmap)
        seq = [AA3to1.get(resmap[rid], "X") for rid in order]
        chain_seq[asym] = ("".join(seq), order)

    rec = max(chain_seq.items(), key=lambda kv: len(kv[1][0]))
    rec_chain, (rec_seq, rec_order) = rec
    return rec_chain, rec_seq, rec_order

def align_local(a, b):
    # crude local alignment using substring search
    best = (-1, -1)
    k = 12
    for i in range(len(a)-k+1):
        seed = a[i:i+k]
        j = b.find(seed)
        if j != -1:
            return j, j+len(a), len(seed)
    return -1, -1, 0

def main():
    cif = pick_one_cif(JOB_DIR)
    rec_chain, rec_seq, rec_order = extract_receptor_chain_sequence_from_cif(cif)
    full_seq = read_full_mik2_from_S3(FULL_S3)

    start_b, end_b, score = align_local(rec_seq, full_seq)
    offset = (start_b + 1)

    print(f"Model: {os.path.basename(cif)}  receptor chain: {rec_chain}")
    print(f"Receptor chain length: {len(rec_seq)}  Full MIK2 length: {len(full_seq)}")
    print(f"Inferred offset = {offset}\n")

    print("First 10 residues (model_idx → paper_pos : modelAA / paperAA):")
    for i in range(10):
        model_idx = rec_order[i]
        paper_pos = offset + i
        print(f"{model_idx:>4} → {paper_pos:>4} : {rec_seq[i]} / {full_seq[paper_pos-1]}")

    print("\nPaper positions check:")
    for pp in PAPER_POSITIONS:
        delta = pp - offset
        if 0 <= delta < len(rec_seq):
            model_idx = rec_order[delta]
            model_aa  = rec_seq[delta]
            paper_aa  = full_seq[pp-1]
            print(f"{pp:>4} → {model_idx:>4} : {paper_aa} / {model_aa}")
        else:
            print(f"{pp:>4} → (out of range)")

if __name__ == "__main__":
    main()

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