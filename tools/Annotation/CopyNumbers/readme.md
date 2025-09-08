# Aggregate Copy Number Calling Script

This script computes **gene-level aggregate copy number (CN)** from transcript-level allele similarity records, using a GENCODE GFF3 gene annotation.

---

## 🔧 Usage

```bash
python script.py -i <input.txt> -g <gencode.gff3> -o <output.tsv>
```

### Arguments

- `-i`, `--input` : Input text file with transcript:gene:similarity records (4th column).
- `-g`, `--gene`  : GENCODE GFF3 annotation file.
- `-o`, `--output`: Output TSV file with aggregate CN summary.

---

## 📥 Input Format

### Input table (`-i`)

A tab-delimited file from SampleAnnotate.py. 


### GENCODE annotation (`-g`)

Standard GENCODE GFF3 annotation file (e.g. `gencode.v33.annotation.gff3`).

Only `HAVANA`-sourced lines are parsed for genes, transcripts, and exons.

We recommend use gencode.v33, and version conflicts might cause problems in some cases. 

---

## 📤 Output

A tab-delimited file with one row per transcript, summarizing CN contribution.

**Note:** The aggregate copy numbers are based on **MANE transcripts**. So some genes may appear on **multiple lines** if multiple transcripts are present in the input. Also, some large genes might be distributed in multiple blocks and may have partial CNVs.

### Columns:

```
genename	id	type	total_size	total_exons	aggregate_copy_number	alleles
```

- `genename`: Transcript ID (ENST)
- `id`: Gene ID (ENSG)
- `type`: Gene biotype (e.g. protein_coding)
- `total_size`: Total length of the MANE transcript
- `total_exons`: Number of exons
- `aggregate_copy_number`: Sum of overlap-normalized similarity
- `alleles`: List of supporting alleles and exon ranges

---

## 🧪 Example

```bash
python script.py \
  -i input_matches.txt \
  -g gencode.v38.annotation.gff3 \
  -o aggregate_cn_summary.tsv
```

---



## 📚 Dependencies

- Python 3.x
- Standard libraries only (no external packages)
