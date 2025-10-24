# README for matrixparition.py

This README provides an overview of how to use `matrixparition.py` and describes its output format in details.

## 1. Requirements
- **Python** version 3.7 or higher.

## 2. Running matrixparition.py

To run the script, use the following command:

```bash
python matrixparition.py -i <Inputfile> -p <number of partitions> 
```
If you only have 16G RAM, you need to partition into two parts;  
If you only have 8G RAM, you need to partition into three parts.  
Then you can concatenate the genotyping output from different parts together.  



Alternatively, if you only interested in certain gene matrix, you may use -g or -G option to select them out:
```bash
python matrixparition.py -i <Inputfile> -g <interested matrices, separated by comma> -G <the file of interested matrix list , separate by line>
```
You can find the matrix information at data/select_files.txt.

for example if you only interested in SMN genes:

```bash
cat Ctyper/data/select_files.txt | grep "SMN"
newGeneCopies/SMN/SMN_partitions/SMN_group1.fa
newGeneCopies/SMN/SMN_partitions/SMN_group1_SMN1OOOSMN2OOOCH17-64J14v7.fa       CH17-64J14.7,AC140134.2,RP11-974F13.3,CH17-133H2.1,RP11-846E15.5,SMN1,SMN2,
```
The output shows SMN can be found in matrix "SMN_group1"

```bash
python matrixparition.py -i HprcCpcHgsvc_final42_matrix.v1.0.txt -g "SMN_group1"
```

If you only interested in HLAs, distract all HLA matrices into a list text file
```bash
cat Ctyper/data/select_files.txt | grep -w "HLA" | cut -f1 | rev | cut -d/ -f1 | rev | cut -d. -f1 > HLA_matrix.list
```

```bash
python matrixparition.py -i HprcCpcHgsvc_final42_matrix.v1.0.txt -G HLA_matrix.list
```

It will output and index partions. 
Each partition can be run separately to reduce RAM usage
