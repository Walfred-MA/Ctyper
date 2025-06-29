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
python3 ctyper_to_fasta.py -i input.txt -r reference1.fa,Allalters.fa -a annotation.tsv.gz -o output.fasta
