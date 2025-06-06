# README for GenotypetoNomenclature.py

This README provides an overview of how to use `GenotypetoNomenclature.py`.

## 1. Requirements
- **Python** version 3.7 or higher.

## 2. Function
It takes the genotype table file from SampleAnnotate.py as input and outputs the public nomenclatures for a list of genes based on annotation results. 

## 3. Running SampleAnnotate.py

To run the script, use the following command:

```bash
python GenotypetoNomenclature.py -i <Inputfile> -a <Nomenclature_Annotation> > result.txt
```

Replace `<Inputfile>` with the path to your input file and `<Sample_Annotation>` with the sample annotation file.
