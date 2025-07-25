# README for GenotypetoVCF.py

This README provides an overview of how to use `GenotypetoVCF.py`.

## 1. Requirements
- **Python** version 3.7 or higher.

## 2. Function
It converts the genotype table file from SampleAnnotate.py to VCF format. Also a BED4 file to indicate the liftovered regions. 

## 3. Running SampleAnnotate.py

To run the script, use the following command:

```bash
python GenotypetoVCF.py -i <Inputfile> -a <Sample_Annotation> -o <Outputfile>
```
Sample_Annotation need to be a tsv file or gzipped tsv.gz file

## 4. Note:
1. **The genotyping regions of ctyper have overlaps between different targeting groups if they are close or share genes, which would cause duplicated variants on VCF, consider picking only one group at a time**
2. **It has to note this conversion does not fully perseve the original table information. Duplications and non-reference sequences may be removed**
3. **The coordinates provide are based on reference liftover. However, it may not be ideal. In the events of duplications and gene conversions, the liftover is always controversy and may subject to different locations regarding different purposes of liftover.**
4. **Genotyping itself has lower base accuracy than NGS alignment, so use by cautious**
5. **Some sequences are masked with lower case, those sequence are likely have high error rates**
6. **It is suggest to use the results to help as phasing templates, or to locate the known variants on its pangenome-alleles. We do not recommend to use in population studies**
     
