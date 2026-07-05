from Bio.PDB import MMCIFParser, DSSP
import numpy as np
import os

cif_path = "/Users/vishnu/Downloads/mik2x4_outputs_0/C_rubella_Carub0001s3337___LIKE-ALDE1/C_rubella_Carub0001s3337___LIKE-ALDE1_model.cif"

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("complex", cif_path)

model = structure[0]

# CHANGE THESE if needed after inspection
receptor_chain_id = "B"
peptide_chain_id = "A"

receptor = model[receptor_chain_id]
peptide = model[peptide_chain_id]

cutoff = 5.0
interface_residues = set()

for res_rec in receptor:
    if res_rec.id[0] != " ":
        continue
    for atom_rec in res_rec:
        for res_pep in peptide:
            if res_pep.id[0] != " ":
                continue
            for atom_pep in res_pep:
                dist = atom_rec - atom_pep
                if dist < cutoff:
                    interface_residues.add(res_rec.id[1])
                    break

interface_residues = sorted(list(interface_residues))

print("Number of interface residues:", len(interface_residues))
print("Residue indices:", interface_residues)



