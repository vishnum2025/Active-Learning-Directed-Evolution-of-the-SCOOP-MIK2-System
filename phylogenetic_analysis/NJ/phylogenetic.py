# # save as build_mik2_phylo_inputs.py and run with: python build_mik2_phylo_inputs.py sd03.txt
# import sys
# from pathlib import Path

# def parse_fasta(txt):
#     entries, header, seq = [], None, []
#     for line in txt.splitlines():
#         if line.startswith(">"):
#             if header:
#                 entries.append((header, "".join(seq)))
#             header = line.strip()[1:]
#             seq = []
#         else:
#             seq.append(line.strip())
#     if header:
#         entries.append((header, "".join(seq)))
#     return entries

# if __name__ == "__main__":
#     s3_path = sys.argv[1] if len(sys.argv) > 1 else "/Users/vishnu/Desktop/v_files/GT/sem1/Research_sem1/sd03.txt"
#     text = Path(s3_path).read_text()
#     entries = parse_fasta(text)

#     # Deduplicate exact sequence duplicates (some datasets include isoforms)
#     seen = set()
#     kept = []
#     for h, s in entries:
#         if s not in seen:
#             kept.append((h, s))
#             seen.add(s)

#     out = "/Users/vishnu/Desktop/v_files/GT/sem1/Research_sem1/MIK2_homologs_full.fasta"
#     with open(out, "w") as f:
#         for h, s in kept:
#             f.write(f">{h}\n")
#             for i in range(0, len(s), 60):
#                 f.write(s[i:i+60] + "\n")
#     print(f"Wrote {len(kept)} sequences to {out}")

from Bio import Phylo
Phylo.draw(Phylo.read("/Users/vishnu/Desktop/v_files/GT/sem1/Research_sem1/MIK2_clustalw.trimmed.ph", "newick"))
# For publication: upload the .dnd/.nwk to iTOL and style it