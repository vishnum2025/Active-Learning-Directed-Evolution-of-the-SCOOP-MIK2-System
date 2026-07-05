from evolvepro.src.process import generate_single_aa_mutants
from Bio import SeqIO

wt_fasta_path = "/Users/vishnu/EvolvePro/c_rubella/C_rubella_Carub0001s3337.fasta"
output_path = "/Users/vishnu/EvolvePro/c_rubella/round1_single_mutants.fasta"

interface_positions = [61, 110, 181, 201, 273, 278, 302, 321, 326]

# Load WT sequence
record = SeqIO.read(wt_fasta_path, "fasta")
wt_seq = str(record.seq)

print("Verifying WT residues at selected positions:\n")

for pos in interface_positions:
    aa = wt_seq[pos - 1]   
    print(f"Position {pos}: {aa}")

print("\nGenerating single mutants...\n")

generate_single_aa_mutants(
    wt_fasta=wt_fasta_path,
    output_file=output_path,
    positions=interface_positions
)

# Confirm number of sequences generated
num_records = len(list(SeqIO.parse(output_path, "fasta")))
print(f"\nTotal sequences written (including WT): {num_records}")