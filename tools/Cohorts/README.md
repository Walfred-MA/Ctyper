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


