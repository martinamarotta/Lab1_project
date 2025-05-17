# Full Workflow

## 1. Data Collection and Filtering

- **UniProt (Swiss-Prot) Search**  
    Retrieved reviewed sequences annotated with:
    
    - _Feature_: `Family and Domains / Domains [FT]: Kunitz`
        
    - _Cross-reference_: `Pfam: PF00014`  
        **Output**: `all_uniprot_kunitz.fasta`
        
- **Full Swiss-Prot Download**  
    All reviewed protein sequences in FASTA format.  
    **Output**: `uniprot_sprot.fasta`
    
- **RCSB PDB Search**  
    Filters applied:
    
    - Pfam domain: `PF00014`
        
    - Resolution ≤ 3.5 Å
        
    - Sequence length: 45–80 residues  
        Exported custom CSV with: PDB ID, sequence, chain ID, and polymer entity ID.  
        **Output**: `pdb_kunitz.csv`
---

## 2. Dataset Preparation

**Convert the PDB report into FASTA format:**

```bash
cat pdb_kunitz.csv | tr -d '"' | \
awk -F ',' '{if (length($2)>0) {name = $2}; print name,$4,$5,$6}' | \
grep PF00014 | awk '{print ">"$1"_"$3; print $2}' > pdb_kunitz.fasta
```

**Redundancy reduction with CD-HIT:**

```bash
cd-hit -i pdb_kunitz.fasta -o pdb_kunitz.clst
mv pdb_kunitz.clst pdb_kunitz_nr.fasta
```

**Manual curation after visual alignment (via MUSCLE):**

* Removed sequence too long (`2ODY_E`)

```bash
awk '/^>/{f=($0 ~ /^>2ODY_E/)?1:0} !f' pdb_kunitz_nr.fasta > pdb_kunitz_senza_2ODY_E.fasta
```

* Removed sequence too short (`5JBT_Y`)

```bash
awk '/^>/ {keep = ($0 != ">5JBT_Y")} keep' pdb_kunitz_senza_2ODY_E.fasta > pdb_kunitz_nr_23.fasta
```

**Final output:**

```
pdb_kunitz_nr_23.fasta
```

---

## 3. Structural Alignment

**Create empty alignment file:**

```bash
touch pdb_kunitz_msa.ali
```

**Multiple structural alignment via PDBeFold** (manual step)

**Convert to uppercase:**

```bash
awk '/^>/ { print; next } { print toupper($0) }' pdb_correct_msa.ali > pdb_uppercase.ali
```

---

## 4. HMM Construction

**Build HMM profile:**

```bash
hmmbuild pdb_kunitz_nr.hmm pdb_uppercase.ali
```

---

## 5. Positive and Negative Dataset Generation

**BLASTP is used to identify highly similar Swiss-Prot Kunitz entries (≥95% identity, ≥50 residues) against all Kunitz:**

```bash
makeblastdb -in all_uniprot_kunitz.fasta -input_type fasta -dbtype prot -out all_uniprot_kunitz.fasta

blastp -query pdb_kunitz_nr_23.fasta -db all_uniprot_kunitz.fasta -out pdb_kunitz_nr_23.blast -outfmt 7
```

**Filter hits with ≥95% identity and ≥50 residues:**

```bash
grep -v "^#" pdb_kunitz_nr_23.blast | awk '$3 >= 95 && $4 >= 50 {print $2}' | sort -u > redundant_ids.txt
```

**Remove redundant sequences:**

```bash
grep "^>" all_kunitz.fasta | cut -d' ' -f1 | cut -c2- > all_kunitz_fasta.ids

comm -23 <(sort all_kunitz_fasta.ids) <(sort redundant.ids) > to_keep.ids

python3 get_ids.py to_keep.ids all_kunitz.fasta > ok_kunitz.fasta
```

**Create negative dataset:**

```bash
grep ">" uniprot_sprot.fasta | cut -d "|" -f2 > sp.ids

comm -23 <(sort sp.ids) <(sort all_kunitz_fasta.ids) > negs.ids

python3 get_ids.py negs.ids uniprot_sprot.fasta > negs.fasta
```

---

## 6. Data Splitting and Model Evaluation

**Split positives (ok_kunitz) and negatives (negs) into two sets each:**

```bash
sort -R negs.ids > random_sp_negs.ids
head -n 286417 random_sp_negs.ids > neg_1.ids
tail -n 286417 random_sp_negs.ids > neg_2.ids

sort -R to_keep.ids > random_ok_kunitz.ids
head -n 183 random_ok_kunitz.ids > pos_1.ids
tail -n 183 random_ok_kunitz.ids > pos_2.ids
```

**Extract sequences using get_ids.py:**

```bash
python3 get_ids.py pos_1.ids uniprot_sprot.fasta > pos_1.fasta
python3 get_ids.py pos_2.ids uniprot_sprot.fasta > pos_2.fasta
python3 get_ids.py neg_1.ids uniprot_sprot.fasta > neg_1.fasta
python3 get_ids.py neg_2.ids uniprot_sprot.fasta > neg_2.fasta
```

**Run HMMER search:**

```bash
hmmsearch -Z 1000 --max --tblout pos_1.out pdb_kunitz_nr.hmm pos_1.fasta
hmmsearch -Z 1000 --max --tblout pos_2.out pdb_kunitz_nr.hmm pos_2.fasta
hmmsearch -Z 1000 --max --tblout neg_1.out pdb_kunitz_nr.hmm neg_1.fasta
hmmsearch -Z 1000 --max --tblout neg_2.out pdb_kunitz_nr.hmm neg_2.fasta
```

> `--max` disables acceleration filters for exhaustive search
> `-Z 1000` sets a fixed database size for e-value normalization

**Extract features for classification:**

```bash
grep -v "^#" pos_1.out | awk '{split($1,a,"|"); print a[2],1,$5,$8}' | tr " " "\t" > pos_1.class
grep -v "^#" pos_2.out | awk '{split($1,a,"|"); print a[2],1,$5,$8}' | tr " " "\t" > pos_2.class
grep -v "^#" neg_1.out | awk '{split($1,a,"|"); print a[2],1,$5,$8}' | tr " " "\t" > neg_1.class
grep -v "^#" neg_2.out | awk '{split($1,a,"|"); print a[2],1,$5,$8}' | tr " " "\t" > neg_2.class
```

**Re-add unmatched negatives with e-value = 10.0 (to balance sets):**

```bash
comm -23 <(sort neg_1.ids) <(cut -f 1 neg_1.class | sort) | awk '{print $1"\t0\t10.0\t10.0"}' >> neg_1.class
comm -23 <(sort neg_2.ids) <(cut -f 1 neg_2.class | sort) | awk '{print $1"\t0\t10.0\t10.0"}' >> neg_2.class
```

---

## 7. Performance Analysis

**Merge datasets:**

```bash
cat pos_1.class neg_1.class > set_1.class
cat pos_2.class neg_2.class > set_2.class
```

**Evaluate model at increasing sensitivity thresholds:**

```bash
for i in $(seq 1 12); do python3 performance.py set_1.class 1e-$i; done
for i in $(seq 1 12); do python3 performance.py set_2.class 1e-$i; done
```


