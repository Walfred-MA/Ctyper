# Ctyper Genotyping FASTA Converter

This script converts Ctyper genotyping annotation outputs into standard FASTA sequence files by extracting the aligned reference regions based on CIGAR strings and reference genome data.
---

## 🔧 Requirements

- Python 3.6+
- No external Python packages (uses only built-in libraries)

Note:
To make sure of inputting all alternative loci, consider also input $CTYPER_PATH/data/Allalters.fa alongside with your reference(s), separated with comma. 
---

## 🚀 Usage

```bash
python getFASTA.py -i <Inputfile> -r <GRCH38_References> -a <Sample_Annotation> -o <Outputfile>

Sample_Annotation need to be a tsv file or gzipped tsv.gz file
References are GRCH38_References (including main chromosomes and all alternative loci) ,can be single or multiple fasta files (comma separated). 
