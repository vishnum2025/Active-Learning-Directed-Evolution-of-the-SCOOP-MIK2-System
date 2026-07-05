# to get aligned FASTA, run this in terminal:

# mafft --auto mik2_ECDs_only.fasta > mik2_alignment.fasta

from Bio import AlignIO
from collections import Counter
import numpy as np

# --- INPUTS ---
alignment_file = "/Users/vishnu/EvolvePro/MIK2_ECD/mik2_alignment.fasta"
reference_id = "C_rubella_Carub0001s3337___LIKE"

interface_positions = [60, 61, 83, 85, 86, 107, 109, 110, 131, 133, 155, 157, 179, 181, 
                       201, 203, 205, 206, 225, 227, 229, 251, 253, 273, 275, 277, 278, 
                       297, 299, 302, 321, 323, 325, 326, 345, 347, 349, 369, 371, 373, 374, 393, 395]

# --- LOAD ALIGNMENT ---
alignment = AlignIO.read(alignment_file, "fasta")
print(f"Sequences in alignment: {len(alignment)}")
print(f"Alignment length: {alignment.get_alignment_length()}")

# --- FIND REFERENCE SEQUENCE ---
ref_index = None
for i, record in enumerate(alignment):
    if reference_id in record.id:
        ref_index = i
        break

if ref_index is None:
    raise ValueError("Reference sequence not found in alignment")

ref_seq = alignment[ref_index].seq

# --- MAP STRUCTURAL RESIDUE NUMBER → ALIGNMENT COLUMN ---
res_to_col = {}
res_counter = 0

for col_idx, aa in enumerate(ref_seq):
    if aa != "-":
        res_counter += 1
        res_to_col[res_counter] = col_idx

# --- ANALYZE VARIABILITY ---
print("\nPosition | Ref_AA | Most Common | Conservation | Unique AAs")

candidate_positions = []

for pos in interface_positions:

    if pos not in res_to_col:
        print(f"{pos:8} | --- missing in reference alignment ---")
        continue

    col = res_to_col[pos]
    column = alignment[:, col]

    residues = [aa for aa in column if aa != "-"]
    counts = Counter(residues)

    if len(residues) == 0:
        continue

    most_common_res, freq = counts.most_common(1)[0]
    conservation = freq / len(residues)
    unique_residues = len(counts)
    ref_aa = ref_seq[col]

    print(f"{pos:8} | {ref_aa:6} | {most_common_res:11} | {conservation:.2f}        | {unique_residues}")

    if conservation < 0.9:
        candidate_positions.append(pos)

print("\nVariable interface positions:")
print(candidate_positions)