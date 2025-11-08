# Active-Learning-Directed-Evolution-of-the-SCOOP-MIK2-System

<img width="565" height="433" alt="MIK2-SCOOP" src="https://github.com/user-attachments/assets/52f41885-3b9f-4764-9997-b1303156b46f" />

## Overview

This project focuses on engineering **synthetic immune-trigger peptides** that activate the plant immune receptor **MIK2** more effectively than natural SCOOP peptides. By integrating **phylogenetics**, **machine learning-based active learning**, and **AlphaFold 3 structural modeling**, this work demonstrates a framework for computational directed evolution of both:

- **Ligands (SCOOP-like peptides)**
- **Receptors (MIK2 variants)**

This project is a part of the **Plant Systems and Synthetic Biology Lab at Georgia Tech** under the guidance of **Dr. Lily Cheung**.

---

## Biological Background

Plants detect molecular signals of damage or infection using receptor-like kinases.  
**MIK2** is one such receptor that recognizes **SCOOP peptides**, small 13-amino-acid danger signals containing a conserved **S–x–S motif**. Pathogens like *Fusarium oxysporum* secrete **SCOOP-like mimic peptides** to inhibit this signaling. This results in an evolutionary arms race between plant immune receptors and microbial peptides.

---

## Computational Workflow

### **1. Phylogenetic Mapping of MIK2**
- Collected MIK2 sequences across Brassicaceae.
- Built a rooted **maximum-likelihood phylogenetic tree**.
- Clustered homologs → selected *Arabidopsis thaliana*, which was the cluster representative of the MIK2 clade for modeling.

### **2. Peptide Library Generation**
- Started with **39 natural SCOOP/SCOOP-like peptides** from *Arabidopsis thaliana* and *Fusarium oxysporum*.
- Generated:
  - 3900 single/double mutants with point mutations
  - 8000 de novo peptides (randomized, but with the functional SxS motif preserved)
- Total search space: **~11,500 peptides**

### **3. Active Learning Directed Evolution (ALDE)**
- Encoded peptides as **260 dimensional one-hot vectors**.
- Trained an ensemble of deep neural networks on seed peptides labeled by AlphaFold 3 predicted binding (iPTM scores) as the target variable "fitness".
- Used **Bayesian Optimization + Thompson Sampling + Hamming distance** to select diverse high-potential peptides from the sequence space.
- Output: **96 high-fitness peptide candidates**

### **4. Structural Screening with AlphaFold 3**
- Modeled all 96 peptides binding to MIK2 using **AlphaFold 3 on GT Phoenix HPC**.
- ~60% of peptides achieved **iPTM ≥ 0.8**, comparable or superior to natural SCOOPs.
- Novel peptides preserved the SxS core but explored new tolerated residues.


---
