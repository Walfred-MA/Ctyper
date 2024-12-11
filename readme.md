<a id="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->


<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/Walfred-MA/Ctyper">
    <img src="images/figure1.png" width="500" >
  </a>

  <h3 align="center">ctyper</h3>

  <p align="center">
    A pangenome allele-specific and copy number specific genotyping tool
    <br />
    <a href="https://www.biorxiv.org/content/10.1101/2024.08.11.607269v1"><strong>Cite the paper»</strong></a>
    <br />
    <br />
    <a href="https://github.com/Walfred-MA/Ctyper/issues/new?labels=bug&template=bug-report---.md">Report Bug</a>
    ·
    <a href="https://github.com/Walfred-MA/Ctyper/issues/new?labels=enhancement&template=give-feedback---.md">Give Feedback</a>
     ·
    <a href="https://github.com/Walfred-MA/Ctyper/issues/new?labels=enhancement&template=ask-for-help---.md">Ask For Help</a>
     ·
    <a href="https://github.com/Walfred-MA/Ctyper/issues/new?labels=enhancement&template=suggest-genes---.md">Suggest include new genes</a>
  </p>
</div>

<!-- TABLE OF CONTENTS -->
## Table of Contents
1. [Getting Started](#getting-started)
2. [About ctyper](#about-ctyper)
3. [Repository Overview](#repository-overview)
4. [Prerequisites](#prerequisites)
5. [Inputs](#inputs)
6. [Installation](#installation)
7. [Usage](#usage)
8. [Results Annotation](#results-annotation)
9. [Results Visualization](#Results-visualization)
10. [Cohort Analysis](#cohort-analysis)
11. [License](#license)
12. [Contact](#contact)


<!-- Getting Started -->
# Getting Started

Let us get start to genotype an example sample NA12718

first obtain its cram file 

```bash
wget "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram"

wget "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram.crai"
```

Download relying databases

```bash
wget "https://zenodo.org/records/14399353/files/HprcCpcHgsvc_final42_matrix.v1.0.txt.gz"
gunzip HprcCpcHgsvc_final42_matrix.v1.0.txt.gz

wget "https://zenodo.org/records/14399353/files/HprcCpcHgsvc_final42_matrix.v1.0.txt.index"

wget "https://zenodo.org/records/14399353/files/PangenomeAlleles_annotationfix.v1.0.tsv.gz"
gunzip PangenomeAlleles_annotationfix.v1.0.tsv.gz
```

Download and compile ctyper, making sure you have gcc >= 8 and eigen, zlib, and HTSlib installed

$EIGEN_ROOT is the root path of EIGEN, usally /usr/

```bash
git clone https://github.com/Walfred-MA/Ctyper 

export EIGEN_ROOT=$EIGEN_ROOT

cd Ctyper/src && make

mv ctyper ..
```

Now, let's us start to genotype it, $LIBRARY_PATH is your path of LD_LIBRARY, if you are using conda to install HTSlib or samtools, then it should be usally at /home/$user_name/miniconda3/lib/, otherwise have a try at /usr/local/lib/

```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH

Ctyper/ctyper -i NA12718.final.cram -m HprcCpcHgsvc_final42_matrix.v1.0.txt -c 1 -b Ctyper/data/backgrounds.list -o ctyper.out > log.txt &
```

After finishing genotyping, let us start to visualize and analyze, making sure you have python3 installed. 

first, change the raw output to a text table. The headers are: 
allele name, allele type, HG38 gene(s), transcriptname(exon_number):transcript_genename:similarity, type category(ref: in same type of HG38 gene, alt: different orthology type, dup: duplicated paralog, novel: diverged paralog), pangenome assembly location, reference liftover location and alignment CIGAR string. 
```bash
python Ctyper/tools/Annotation/SampleAnnotate.py -i ctyper.out -a PangenomeAlleles_annotationfix.v1.0.tsv -o genotype.txt
```

second, if you want to have genotypes of HLA, CYP2D6, and KIR. 
```bash
python Ctyper/tools/Annotation/Nomenclature/GenotypetoNomenclature.py -i  genotype.txt -a Ctyper/data/all_nomenclature.txt
```

example output: (the first column are names of pangenome allele, the second column is their public Nomenclature)
```
CYP2D_group1_HG01358_h2_162 *1
CYP2D_group1_GW00005_h1_12 *1
HLA-A_group2_HG01071_h2_1333 HLA-A*03:01:01:01
HLA-A_group2_HG00741_h1_1304 HLA-A*11:01:01:01
HLA-B_group10_HG02717_h2_2635 HLA-B*15:03:01:01
HLA-B_group10_HG02055_h2_2334 HLA-B*07:02:01:14
HLA-C_group8_chr6_2356 HLA-C*07:02:01:03
HLA-DM_group1_HG01361_h1_160 HLA-DMA*01:01:01:13,HLA-DMB*01:04
HLA-DM_group1_NA12878_h1_216 HLA-DMA*01:03:01:02,HLA-DMB*01:07:01:02
HLA-DO_group1_NA12878_h2_435 HLA-DOB*01:01:01:24
HLA-DO_group1_HG02572_h1_374 HLA-DOB*01:01:01:11
HLA-DO_group2_chr6_GL000251v2_alt_485 HLA-DOA*01:01:05
HLA-DO_group2_chr6_GL000253v2_alt_490 HLA-DOA*01:01:01:01
HLA-DP_group2_HG00733_h2_790 HLA-DPA1*02:01:01:01,HLA-DPB1*11:01:01:01
HLA-DP_group3_HG00733_h2_789 HLA-DPB2*03:01:01:03
HLA-DQ_group84_HG01071_h1_58107 HLA-DQA2*01:04:01:01,HLA-DQB2*01:01:01:04
HLA-DRA_group1_HG02818_h2_198 HLA-DRA*01:02:02:11
HLA-DRA_group1_chr6_GL000253v2_alt_244 HLA-DRA*01:01:01:02
HLA-E_group6_GW00048_h1_639 HLA-E*01:01:01:01
HLA-E_group6_chr6_GL000255v2_alt_2012 HLA-E*01:01:01:52
HLA-F_group1_HG01258_h2_666 HLA-F*01:03:01:01
HLA-F_group1_GW00057_h1_494 HLA-F*01:01:02:09
HLA-G_group1_HG01071_h2_771 HLA-G*01:01:01:05
HLA-G_group1_GW00038_h1_416 HLA-G*01:01:03:03
HLA-H_group3_HG02630_h2_1530 HLA-H*02:04:01
HLA-H_group3_NA12878_h2_2195 HLA-H*02:07:01:01
HLA-J_group1_GW00027_h1_325 HLA-J*01:01:01:04
HLA-J_group1_GW00043_h1_497 HLA-J*01:01:01:02
HLA-K_group2_HG01258_h2_650 HLA-K*01:01:01:01
HLA-L_group1_GW00022_h2_200 HLA-L*01:03
HLA-L_group1_GW00052_h2_462 HLA-L*01:03
HLA-N_group1_NA24385_h1_230 HLA-N*01:01:01:01
HLA-P_group1_HG01071_h2_425 HLA-P*02:01:01:02
HLA-S_group1_GW00018_h1_38 HLA-S*01:01:01:02
HLA-S_group1_chr6_GL000253v2_alt_241 HLA-S*01:01:01:01
HLA-T_group1_chr6_463 HLA-T*03:01:01:01
HLA-U_group1_HG00732_h1_347 HLA-U*01:05
HLA-U_group1_HG01175_h1_389 HLA-U*01:04
HLA-V_group1_HG01071_h2_641 HLA-V*01:01:01:01
HLA-W_group1_HG01071_h2_308 HLA-W*01:01:01:05
KIR2DL_group1_HG01952_h1_188 KIR2DL1*0040101,KIR2DL2*0010102,KIR2DL4*0080204,KIR2DL5B*0020109
KIR2DL_group1_chr19_KI270933v1_alt_348 KIR2DL1*0020103,KIR2DL3*0020103,KIR2DL4*0080204
KIR3DL1_group1_HG01123_h2_259 KIR3DL1*0040201
KIR3DL1_group1_HG00741_h2_246 KIR3DL1*0040101
```

Then let Plot SMN genotyping results as mutant map (make sure you have pandas and matplotlib installed)

(optional) download genecode gff3 annotation
```bash
wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gff3.gz"
zcat gencode.v47.annotation.gff3.gz | grep "gene_name=SMN" > SMN.gff3
```

Get names of genotyped SMN alleles
```bash
grep "^SMN" genotype.txt | cut -f1 | tr "$\n" "," | sed  's/.$/\n/'

SMN_group1_GW00024_h2_152,SMN_group1_chr5_718,SMN_group1_chr5_719,SMN_group1_HG01175_h1_475
```


Otain all SMN annotation from database
```bash
grep "^SMN" PangenomeAlleles_annotationfix.v1.0.tsv > SMN_annotation.txt  &
```

Plot 
```bash
python Ctyper/tools/Plot/typemutant.py -i SMN_allannotations.txt  (optional) -g SMN.gff3 -n SMN_group1_GW00024_h2_152,SMN_group1_chr5_718,SMN_group1_chr5_719,SMN_group1_HG01175_h1_475 -o SMN.png
```

example output of SMN mutant map: 

Each row is each allele, and each vertically location is the consensus in multiple sequence alignments.  
The gap is the deletion/missing sequences on each allele and each black dots is each variants on HG38.  
Each color represents each allele-type we defined.  
The location of genotyped alleles are slightly bolded in color (may need to zoom in sometimes) and highlighted with a red dot on the left.  
(optional) The genecode gene and exons are highlight at the top.

![Alt Text](images/exampleSMN.png)




<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ABOUT THE PROJECT -->
# About ctyper

Ctype is a command line tool to perform copy number sensitive versatile genotyping for NGS (Next-Generation Sequencing) data using pangenome database. It is designed to work on complex CNV (copy number variation) genes, but can also work as genotyping, local phasing or SV-calling tools for other genes. 

The results will be represented as pangenome alleles, which is defined as genic segments with locally phased variants that are combinatorially heritable mostly range in 15-50 kb, about the size of "Haplotype blocks" or "LD-blocks". The pangenome alleles are further classified as allele-types among populations, which are highly similar subgroups, and can be used to represent complex genetic variations such as structrual variations, gene conversion, duplication, translocation and etc. 

Ctyper is highly efficient, accurate and visualizable, thus allows high resolution large cohort association studies on complex CNV genes as well as complex genetic variations. 

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- GETTING STARTED -->
# Repository Overview

This repository includes the following components:

1. **Ctyper (src/)**  
   A C++ program for performing genotyping on NGS data.
2. **Auxiliary Analyze Tool (tools/Annotation/)**  
    Python tools to analyze genotyping results.
    - `tools/Annotation/SampleAnnotate.py`: Convert the raw output file to more interpretable table with annotations.
    - `tools/Annotation/Nomenclature/GenotypetoNomenclature.py`: Output the public nomenclatures from genotyping results. The nomenclatures annotation can be found at Data/all_nomenclature.txt.
    - `tools/Annotation/VCF/GenotypetoVCF.py`: Converting the genotyping results to VCF format. This requires Individual Sample Annotation Database (eg. PangenomeAlleles_annotationfix.v1.0.tsv). Note: this not recommended to be used in association studies, because of the limitation of VCF file on representing pangenome, but maybe used for locating known important variants.
      
4. **Visualization Tool (tools/Plot/)**  
   A Python tool for visualizing genotyping results using multiple sequence alignments.

5. **Cohort Analysis Tool (tools/Cohort/)**  
   A Python tool for summarizing results across all samples in a cohort study, including annotation information.

6. **Pangenome Allele Database**  
   A database required by ctyper for genotyping. Files: `HprcCpcHgsvc_final42_matrix.v1.0.txt` and `HprcCpcHgsvc_final42_matrix.v1.0.txt.index` available at [https://zenodo.org/records/13381931](https://zenodo.org/records/13381931).

7. **Population Analysis Annotation Database**  
   A database containing annotation information to support population-level analysis. File: `PangenomeAlleles_typefix.v1.0.tsv` available at [https://zenodo.org/records/13381931](https://zenodo.org/records/13381931).

8. **Individual Sample Annotation Database**  
   A comprehensive database with detailed annotations for individual pangenome alleles, supporting individual sample studies and visualization. File: `PangenomeAlleles_annotationfix.v1.0.tsv` available at [https://zenodo.org/records/13381931](https://zenodo.org/records/13381931).

9. **Additional Data Files (data/)**  
   - `backgrounds.list`: A list of k-mers used as backgrounds to determine NGS coverage if not predetermined.
   - `select_files.txt`: The catalog for all included genes and matrices, which can be used to locate the gene of interest.
   - `all_nomenclature.txt`: The public nomenclatures used by GenotypetoNomenclature.py, currently including HLAs, CYP2D6, and KIRs.  
     - **HLAs include**: HLA-A, HLA-B, HLA-C, HLA-DMA, HLA-DMB, HLA-DOA, HLA-DOB, HLA-DPA1, HLA-DPA2, HLA-DPB1, HLA-DPB2, HLA-DQA1, HLA-DQA2, HLA-DQB1, HLA-DQB2, HLA-DRA, HLA-E, HLA-F, HLA-G, HLA-H, HLA-J, HLA-K, HLA-L, HLA-N, HLA-P, HLA-S, HLA-T, HLA-U, HLA-V, HLA-W, HLA-Y.
     - **KIRs include**: KIR2DL1, KIR2DL2, KIR2DL3, KIR2DL4, KIR2DL5A, KIR2DL5B, KIR2DP1, KIR2DS1, KIR2DS2, KIR2DS3, KIR2DS4, KIR2DS5, KIR3DL1, KIR3DL2, KIR3DL3.

10. **Test Cases (tests/)**  
   Simple test cases for validating the tools and pipeline.




<!-- Prerequisites -->
# Prerequisites  

Ctyper is officially supported only in a **Linux** environment. Although it may run on UNIX-like system like MACOS as well, but we may not provide support for it.  

### System Required
You need RAM > 20G to run all genes at once. 
If you have a less RAM, you may split the whole database to smaller partitions. See in tools/Partition/.


### Required Software

1. **GCC 8**
2. **Eigen Library**  
   [https://gitlab.com/libeigen/eigen](https://gitlab.com/libeigen/eigen)
3. **HTSlib**  
   [https://github.com/samtools/htslib](https://github.com/samtools/htslib)  
   Alternatively, you can install the entire Samtools package:  
   [https://github.com/samtools/samtools](https://github.com/samtools/samtools)

### Python Tools (Optional)

For additional Python tools, you need:

1. **Python 3.7+**
2. **NumPy**
3. **Pandas**  
   [https://pandas.pydata.org/](https://pandas.pydata.org/)
4. **Matplotlib**  
   [https://matplotlib.org/stable/install/index.html](https://matplotlib.org/stable/install/index.html)

  

<!-- Inputs-->
# Inputs  
  
Ctyper takes five types of files as input:  
  
1. **CRAM files** (`*.cram`) — must be indexed (recommended due to lower I/O intensity). 
2. **BAM files** (`*.bam`) — must be indexed.
3. **SAM files** (`*.sam`).
4. **FASTQ files** (`*.fastq`).
5. **FASTA files** (`*.fa`, `*.fasta`).
  
<!-- Installation -->
# Installation  
  
1. **Install all prerequisites.**
2. **Navigate to the `src/` directory:**

   ```bash
   cd src/
   ```

3. **Compile ctyper:**

   ```bash
   make
   ```

4. **Move ctyper to your installation folder:**

   ```bash
   mv ctyper /path/to/your/install/folder
   ```

5. **Download the pangenome allele database file** (we'll refer to it as `$Database`) **and its index file.**
6. **(Optional) Download the allele-type annotation file.**
7. **(Optional) Download the full alleles annotation file.**
8. **You're ready to go!**


<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
# Usage

Ctyper requires:

1. **Input file(s)**
2. **Corresponding output file(s)**
3. **The pangenome allele database file** (must be indexed)
   - the database of pangenome allele for 3,351 CNV genes and 212 medically challenging genes. The information of those genes and their belonged matrices can be found at data/select_files.txt.
   - If you prefer only run selected genes or you want to run in smaller RAM, you can use a tool called matrixpartion.py we included at tools/Partition
   
4. **Sequencing coverage information**, either by:
   - Providing background k-mers for ctyper to determine sequencing coverage (recommended)
   - Directly providing sequencing coverage information (useful when you don't have WGS data)



## Running Ctyper

Ctyper can process:

### 1. A Single File

Using background k-mers:

```bash
ctyper -i $Inputfile -m $Database -o $Outputfile -b background.list -n 1
```

Or providing sequencing coverage of 31-mers: 
31-mers = read_coverage * (read_length - 30) / read_length, for example a NGS at read coverage = 30, its 31-mers coverage is 30 * (150-30)/150  = 24.

```bash
ctyper -i $Inputfile -m $Database -o $Outputfile -d $sequencing_coverage -n 1
```

If you only interested in genes at certain region, for example chr1:100-1000

```bash
samtools view -b input.bam chr1:100-1000 -f 4 -o subset.bam
samtools sort -o sorted_input.bam input.bam
samtools index sorted_input.bam

ctyper -i sorted_input.bam -m $Database -o $Outputfile -d $sequencing_coverage -n 1
```


### 2. A Cohort of Files

Using background k-mers:

```bash
ctyper -I $AllInputs -m $Database -o $AllOutputs -b background.list -n $threads
```

Or providing sequencing coverage:

```bash
ctyper -I $AllInputs -m $Database -o $AllOutputs -D $All_sequencing_coverages -n $threads
```

- **$AllInputs**: A text file where each line is the path to an input file.
- **$AllOutputs**: A text file where each line is the output file corresponding to the input file (same line number).
- **$All_sequencing_coverages**: A text file where each line is the sequencing coverage information for the corresponding input file (same line number).

---

## Parameters

Supported parameters:

- **Help:**
  - `-h`: Print help information.

- **Inputs:**
  - `-i <string>`: Path to an individual input file.
  - `-I <string>`: Path to a file listing multiple input files (one per line).

- **Database:**
  - `-m <string>`: Path to the database file used for genotyping.
  - `-b <string>`: Path to the background k-mer list.

- **Coverage Information** (cannot be used with `-b`):
  - `-d <float>`: Sequencing coverage of the input file.
  - `-D <string>`: Path to a file listing sequencing coverages (one per input file, corresponding by line number).

- **Outputs:**
  - `-o <string>`: Path to the individual output file.
  - `-O <string>`: Path to a file listing multiple output files (one per input file, corresponding by line number).

- **Multithreading:**
  - `-n <int>`: Number of threads to use (default is 1).

- **Bias Correction:**
  - `-c <bool>`: Perform bias correction for Illumina data.


---

# Results Annotation

Here shows the commands to make the genotyping results interpretable:

1. **Summary the genotyping results into a table with annotation information**

   ```bash
   python tools/Annotation/SampleAnnotate.py -i $ctyper_outputs.txt -a PangenomeAlleles_annotationfix.tsv > genotype.txt
   ```

2. **(optional) Obtain public nomenclatures for important genes,including HLA, CYP2D6, and KIR**

   ```bash
   python tools/Annotation/Nomenclature/GenotypetoNomenclature.py -i genotype.txt -a data/all_nomenclature.txt > nomenclature.txt
   ```
3. **(optional) Convert genotyping results into VCF file**

   ```bash
   python tools/Annotation/VCF/GenotypetoVCF.py -i genotype.txt -o genotype.vcf
   ```

---
---

# Results Visualization

Visualization is performed on a **gene-by-gene** basis (not genome-wide).

For example, to visualize the gene **AMY1A**:

1. **Identify the gene group for AMY1A:**

   ```bash
   cat data/select_files.txt | grep -w "AMY1A" | cut -f2,3
   ```

   Output:

   ```
   newGeneCopies/AMY/AMY_partitions/AMY_group1_AMY1BOOOAMY1COOOAMY1A.fa	AMY1B,AMY1C,RNPC3,AMY1A,AMYP1,ACTG1P4,RP5-1108M17.5,AMY2B,AMY2A,
   ```

   This shows that **AMY1A** is in **AMY_group1**, along with other amylase genes.

2. **Extract the annotation for AMY_group1 from the full annotation table:**

   ```bash
   grep "^AMY_group1_" PangenomeAlleles_annotationfix.tsv > AMY_group1_annotationfix.tsv
   ```

3. **Extract the genotyping results for AMY_group1:**

   ```bash
   grep "^result: AMY_group1_" genotype.txt
   ```

   Output:

   ```
   result: AMY_group1_GW00031_h1_556,AMY_group1_GW00051_h2_891,
   ```

   The genotyping result is: `AMY_group1_GW00031_h1_556,AMY_group1_GW00051_h2_891`.

4. **Visualize the results:**

   ```bash
   python typemutant.py -i AMY_group1_annotationfix.tsv -n "AMY_group1_GW00031_h1_556,AMY_group1_GW00051_h2_891" -o output.png
   ```

**Optional:** To visualize the GENCODE annotation on the MSA:

1. **Obtain the GENCODE annotation:**

   ```bash
   grep "gene_name=AMY" genecode.gff3 > AMY.gff3
   ```

2. **Run the visualization with GENCODE annotation:**

   ```bash
   python typemutant.py -i AMY_group1_annotationfix.tsv -g AMY.gff3 -n "AMY_group1_GW00031_h1_556,AMY_group1_GW00051_h2_891," -o output.png
   ```


## Cohort Analysis

There are two scripts in the `tools/Cohort` folder that work together for cohort analysis.

1. **Download the allele-type annotation table:**

   ```
   PangenomeAlleles_typefix.tsv
   ```

2. **Run `CountAllele.py` on each sample to get allele-type results:**

   ```bash
   # For each result in $results
   for result in $results; do
       python CountAllele.py -i $result -t PangenomeAlleles_typefix.tsv -o "${result}_alleletype.out"
   done
   ```

   **Or**, to run in parallel:

   ```bash
   python CountAllele.py -f $results_folder/ -t PangenomeAlleles_typefix.tsv -n $numthreads
   ```

3. **Summarize results into a single file and add annotations:**

   ```bash
   python SummaryAlleles.py -f $results_folder/ -t PangenomeAlleles_typefix.tsv -o cohort_results.tsv
   ```

---





<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Walfred MA - wangfeim@usc.edu
Mark Chaisson - mchaisso@usc.edu
<p align="right">(<a href="#readme-top">back to top</a>)</p>






