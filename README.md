
# Building a profile HMM for the Kunitz-Type Protease Inhibitor Domain (PF00014)

## Overview

This project provides a computational pipeline to build a Profile Hidden Markov Model (HMM) based on structural data, focusing on the Kunitz-type protease inhibitor domain (Pfam ID: PF00014), a small, evolutionarily conserved protein motif that serves as the active region of protease inhibitors. This project was developed as part of an assignment for the Laboratory of Bioinformatics 1 course, within the Master's Degree in Bioinformatics at the University of Bologna. 
The pipeline integrates:
- Data Collection and Filtering
- Structural Alignment
- HMM construction
- Model testing and performance evaluation

For details on data cleaning, manual curation, and the full processing workflow, see `docs/workflow_full.md`

---

## Pipeline Summary

### 1. **Data Collection**

- Go to [UniProt](https://www.uniprot.org/) :
	- Retrieve reviewed protein sequences containing the Kunitz domain (Pfam: PF00014): uniprotkb_ft_domain_Kunitz_AND_xref_pfa_2025_05_05.fasta
	- Download all proteins sequences from the SwissProt database (https://www.uniprot.org/help/downloads) in FASTA format: uniprot_sprot.fasta 
    
- Collect **3D structures** from [RCSB PDB](https://www.rcsb.org/) containing PF00014, filtered by:
    
    - Resolution ≤ 3.5 Å
        
    - Sequence length 45–80 residues
    Output: rcsb_pdb_custom_report_20250505053111.csv

### 2-3. **Dataset Preparation and Structural Alignment**

- Convert the rcsb_pdb_custom_report_20250505053111.csv file containing the structural data into a FASTA format: pdb_kunitz.fasta
    
- Remove redundant sequences using **CD-HIT**: pdb_kunitz_nr.fasta
    
- Visualization of the alignment using MUSCLE (https://www.ebi.ac.uk/jdispatcher/msa/muscle?stype=protein) followed by manual curation (if necessary) of the sequences: pdb_kunitz_nr_23.fasta
    
- Align final sequences with **PDBeFold**: pdb_uppercase.ali

### 4. **Build the HMM**

- Use [**HMMER**](http://hmmer.org/) to build the profile HMM: pdb_kunitz_nr.hmm

### 5. **Positive and Negative Dataset Generation**

- Create a database containing all sequences of the uniprotkb_ft_domain_Kunitz_AND_xref_pfa_2025_05_05.fasta file by using BLAST.
    
- Use BLASTP to **detect and flag redundant sequences** (>95% identity, ≥50 aligned residues). Obtain the file without the redundant sequences using the Python script [`get_ids.py`](get_ids.py) : ok_kunitz.fasta

- The negative dataset contains all the proteins of the uniprot_sprot.fasta except the proteins containing the Kunitz domain: negs.fasta

### 6. **Data splitting and Model Evaluation**

- Split both sets into two subsets each: pos_1.fasta, pos_2.fasta, neg_1.fasta, neg_2.fasta

- The trained HMM is used to scan both positive and negative datasets using HMMER (`hmmsearch`), extracting e-values and bit scores as features. Sequences with no match are reintegrated with default scores to preserve dataset completeness. Extract the **sequence ID**, **e-value**, and **bit score** from the `hmmsearch` output (`--tblout`) by filtering out comment lines (starting with `#`) and selecting the relevant columns. These values are saved in `.class` files: pos_1.class, pos_2.class, neg_1.class, neg_2.class

### 7. **Performance Analysis**

- two evaluation sets are created by concatenating the respective positive and negative `.class` files:

	- `set_1.class` = `pos_1.class` + `neg_1.class`
    
	- `set_2.class` = `pos_2.class` + `neg_2.class`

- The script [`performance.py`](performance.py) is used to assess model performance across a range of e-value thresholds (`1e-1` to `1e-12`). It reports standard classification metrics such as accuracy (Q2), Matthews Correlation Coefficient (MCC), True Positive Rate (TPR), False Positive Rate (FPR), and Positive Predictive Value (PPV).

## Requirements

- [HMMER](http://hmmer.org/)
    
- [CD-HIT](https://github.com/weizhongli/cdhit)
    
- [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
    
- Python 3 (for `get_ids.py`, `performance.py`)

## External Resources

- [UniProt](https://www.uniprot.org/)
    
- [RCSB PDB](https://www.rcsb.org/)
    
- PDBeFold
    
- MUSCLE
    
- WebLogo, [Skylign](https://skylign.org/)
