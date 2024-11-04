# README for SampleAnnotate.py

This README provides an overview of how to use `SampleAnnotate.py` and describes its output format in details.

## 1. Requirements
- **Python** version 3.7 or higher.

## 2. Running SampleAnnotate.py

To run the script, use the following command:

```bash
python SampleAnnotate.py -i <Inputfile> -a <Sample_Annotation> > result.txt
```

Replace `<Inputfile>` with the path to your input file and `<Sample_Annotation>` with the sample annotation file.

## 3. Output Description

The output will be a tab-separated values (TSV) file where each row contains information about a genotyped pangenome allele. The columns of the output table include:

- **Name**
- **Allele Type**
- **Overlap Genes**
- **Mapped Transcripts**
- **Relationship to Reference Genes**:
  - **Ref**: same allele type as reference.
  - **Alt**: different allele type.
  - **Dup**: allele types of similar paralogs.
  - **Novel**: allele types of diverged paralogs.
- **Location in Pangenome Assemblies**:
  - Format: `contig:start-end[strand]`
- **Liftover on the Reference**:
  - Format: `chromosome:start-end[strand]:CIGAR-extended`
  - **CIGAR-extended**: extended version of CIGAR strings (see section below).

## 4. CIGAR-Extended Description for Liftover Pairwise Alignment

The **CIGAR-extended** format is an enhanced version of the CIGAR string, allowing it to encode all variant information within the alignment string.

- For each variant (e.g., substitution, insertion, deletion), variant details are appended to the CIGAR segment.
  - For example, `6D` becomes `6DAAATTT`, meaning `AAATTT` was deleted.
- For mismatches, both the original and altered sequences are included sequentially.

You can convert CIGAR-Extended to normal CIGAR simply via python:

CIGAR = "".join(re.findall(r'\d+[=MIDXSH]',CIGAR-Extended))






