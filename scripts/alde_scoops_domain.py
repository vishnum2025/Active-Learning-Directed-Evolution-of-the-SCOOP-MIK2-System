#!/usr/bin/env python3
"""
Generate peptide domain for ALDE optimization:
 - Enforces the conserved SxS motif
 - Includes known SCOOP and Fusarium peptides
 - Adds random and mutational diversity
"""

import random
import pandas as pd

AA20 = list("ACDEFGHIKLMNPQRSTVWY")

# ---------------- CONFIG ----------------
L = 13                    # SCOOP length
N_RANDOM = 8000           # number of random SxS peptides
N_MUTANTS_PER_BASE = 100  # number of mutants per seed peptide
BASE_PATH = "/Users/vishnu/v_files/GT/sem1/Research_sem1/scoop_mik2_af_pipeline/alde_athaliana"
TRAIN_FILE = f"{BASE_PATH}/alde_training_round0.csv"
OUT_FILE = f"{BASE_PATH}/alde_domain.csv"
# ----------------------------------------

# --- Load your known SCOOP and Fusarium peptides ---
df = pd.read_csv(TRAIN_FILE)
seed_peptides = [s.strip().upper() for s in df["sequence"]]
print(f"Loaded {len(seed_peptides)} seed peptides.")

domain = set(seed_peptides)

# --- Mutational variants (centered around real peptides) ---
for base in seed_peptides:
    base = list(base)
    # force motif fix before mutating
    if len(base) != L:
        continue
    base[4], base[6] = "S", "S"  # enforce canonical motif
    for _ in range(N_MUTANTS_PER_BASE):
        new = base.copy()
        for pos in random.sample([i for i in range(L) if i not in (4,6)], k=random.choice([1,2])):
            new[pos] = random.choice(AA20)
        domain.add("".join(new))

# --- Random peptides (de novo, but SxS fixed) ---
for _ in range(N_RANDOM):
    seq = [random.choice(AA20) for _ in range(L)]
    seq[4], seq[6] = "S", "S"
    domain.add("".join(seq))

# --- Finalize and save ---
domain = pd.DataFrame({"sequence": sorted(domain)})
domain.to_csv(OUT_FILE, index=False)
print(f"Saved domain: {len(domain)} peptides → {OUT_FILE}")
