# Ctyper Genotyping FASTA Converter

This script converts Ctyper genotyping annotation outputs into standard FASTA sequence files by extracting the aligned reference regions based on CIGAR strings and reference genome data.
---

## 🔧 Requirements

- Python 3.6+
- No external Python packages (uses only built-in libraries)

Note:
1. To make sure of inputting both GRCH38 and CHM13 references because some sequences can only be mapped to CHM13. References need to be separated with comma.
2. CHM13 reference needs to have chromosome names with NC_* instead of chr*, to distinguish from GRCH38.
3. Also consider including all alternative loci and input $CTYPER_PATH/data/Allalters.fa alongside with your references,
4. e.g:  -r HG38_main.fa,AllAlters.fa,CHM13.fa
---

## 🚀 Usage

```bash
python getFASTA.py -i <Inputfile> -r <References> -a <Sample_Annotation> -o <Outputfile>

Sample_Annotation need to be a tsv file or gzipped tsv.gz file
References are GRCH38_References (including main chromosomes and all alternative loci) ,can be single or multiple fasta files (comma separated).
e.g 
```bash
python getFASTA.py -i genotype.tsv -r HG38_main.fa,AllAlters.fa,CHM13.fa -a PangenomeAlleles_annotationfix.tsv.gz -o output.fa

To study the unmapped sequences:
awk -F'\t' '/^>/{keep = ($3 == "Unmap")} keep' <Outputfile> > unmap_sequences.fa
