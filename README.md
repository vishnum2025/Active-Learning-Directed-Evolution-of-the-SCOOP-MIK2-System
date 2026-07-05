# Active-Learning-Directed-Evolution-of-the-SCOOP-MIK2-System

<img width="565" height="433" alt="MIK2-SCOOP" src="https://github.com/user-attachments/assets/52f41885-3b9f-4764-9997-b1303156b46f" />

## Overview

This project focuses on engineering **synthetic immune-trigger peptides** that activate the plant immune receptor **MIK2** more effectively than natural SCOOP peptides. By integrating **phylogenetics**, **machine learning-based active learning**, and **AlphaFold 3 structural modeling**, this work demonstrates a framework for computational directed evolution of both:

- **Ligands (SCOOP-like peptides)**
- **Receptors (MIK2 variants)**

This project is a part of the **Plant Systems and Synthetic Biology Lab @ Georgia Institute of Technology** under the guidance of **Dr. Lily Cheung**.

---

## Biological Background

Plants detect molecular signals of damage or infection using receptor-like kinases.  
**MIK2** is one such receptor that recognizes **SCOOP peptides**, small 13-amino-acid danger signals containing a conserved **S–x–S motif**. Pathogens like *Fusarium oxysporum* secrete **SCOOP-like peptides** to inhibit this signaling. This results in an evolutionary arms race between plant immune receptors and microbial peptides.

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

### **5. Receptor Selection for Orthogonality**

- Screened all MIK2 clade orthologs and paralogs against top ALDE peptides (ALDE1, ALDE2) and native SCOOPs (SCOOP10B, SCOOP12) using AlphaFold 3.
- Selected via delta fitness: Δfitness = max(iPTM_ALDE) − max(iPTM_SCOOP).
- C. rubella MIK2 (Carub0001s3337) had the highest delta fitness → selected as engineering chassis (CrMIK2).

### **6. Interface Mapping**

- Identified 43 receptor residues within 5 Å of CrMIK2–ALDE1 complex (AlphaFold 3 structure).
- Filtered against MSA of all MIK2 orthologs for evolutionary variability.
- Intersection yielded 9 candidate mutagenesis sites: 61, 110, 181, 201, 273, 278, 302, 321, 326.

### **7. Single Mutant Screening and Additive Model**

- Modeled all single amino acid substitutions at the 9 positions against ALDE1 and SCOOP12 via AlphaFold 3.
- Built additive null model to test independence hypothesis of interface residues.
- Quantified epistasis via logit-transformed metric to correct bounded-scale artifacts.

### **8. EVOLVEpro for Double Mutant Generation**

- Trained random forest regression on ESM-2 embeddings of 171 single mutant sequences.
- Ranked all 12,996 possible double mutants.
- Top 20 EVOLVEpro nominations validated via AlphaFold 3: 18/20 met selectivity threshold (ALDE1 iPTM ≥ 0.80, SCOOP12 iPTM < 0.80).
- Tested EVOLVEpro nominations and additive model top predictions against all 14 ALDE peptides and 27 native SCOOPs.



---


_All the cited work has been mentioned in the References.md file_
