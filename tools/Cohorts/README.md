# README for CountAllele.py and SummaryAlleles.py

This README provides an overview of how to use `CountAllele.py`, and `SummaryAlleles.py`,as well as describes its output format in detail.

## 1. Requirements
- **Python** version 3.7 or higher.
- **Numpy**.
- **Pandas**.

## 2. Function Description
- **CountAllele.py**: convert genotyping results into copy numbers of allele-types
- **SummaryAlleles.py**: summary cohort outputs of CountAllele.py and include annotation

## 3. Running CountAllele.py

To run the script, use the following command:

```bash
python CountAllele.py -i <Inputfile> -t <Allele_type_Annotation> -o <Outputfile>
```
or
```bash
python CountAllele.py -f <results_folder> -t <Allele_type_Annotation> -n <numthreads>
```

## 3. Running SummaryAlleles.py

To run the script, use the following command:

```bash
python SummaryAlleles.py -f <results_folder> -t <Allele_type_Annotation> -o <Outputfile>
```

## 4. Output description

It output the cohort copy number table that summarizes the genotyped copy number of each sample for each allele-type.

The row represents each allele-type.  

The first 10 columns are for annotations. The rest columns are for individuals.  


- **Name**
- **Allele Type**
- **Relationship to Reference Genes**:
  - **Ref**: same allele type as reference.
  - **Alt**: different allele type.
  - **Dup**: allele types of similar paralogs.
  - **Novel**: allele types of diverged paralogs.
- **Overlap Genes**
- **Mapped transcripts counts**
- **Relationship to Reference Genes**:
  - **Ref**: same allele type as reference.
  - **Alt**: different allele type.
  - **Dup**: allele types of similar paralogs.
  - **Novel**: allele types of diverged paralogs.
- **Most likely Location in Pangenome Assemblies**:
  - Format: `ref_allelename_overlap:contig:start-end[strand]`
- **Additional annotation**:
  - Format: `Other smaller allele-types that this type captulate | other annotion information like gene consertion`
- **SV annotation**:
  - Insertion or deletion found, as well as their number of unique bases.
- **members**:
  - List of all its pangenome-allele members

## 5. Rest columns are copy numbers of samples. 












