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
    <img src="images/ctyperlogo.png" width="500" height="500" >
  </a>

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
1. [About ctyper](#about-ctyper)
2. [Getting Started](#getting-started)
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


<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ABOUT THE PROJECT -->
# About ctyper

Ctyper is a versatile genotyping tool designed for copy number-sensitive analysis of NGS (Next-Generation Sequencing) data using a pangenome database. It excels at genotyping complex CNV (copy number variation) and highly polymorphic genes, but can also be used for standard genotyping, local phasing, and structural variant (SV) calling across other genes.

Ctyper enables rapid genotyping of target genes with state-of-the-art accuracy, often completing each target in seconds, and is fully scalable for biobank-level analyses.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- Getting Started -->
# Getting Started

ctyper now is available in conda:

```bash
conda create -n ctyper_env  -c conda-forge -c bioconda ctyper
conda activate ctyper_env
```

Download relying databases

```bash
wget "https://zenodo.org/records/16340156/files/HprcCpcHgsvc_cmr_matrix.txt.gz"
wget "https://zenodo.org/records/16340156/files/HprcCpcHgsvc_cmr_matrix.txt.index"
wget "https://zenodo.org/records/16340156/files/HprcCpcHgsvc_cmr_matrix.txt.bgd"
wget "https://zenodo.org/records/16340156/files/PangenomeAlleles_annotationfix.tsv.gz"
gunzip HprcCpcHgsvc_final42_matrix.v1.0.txt.gz
gunzip PangenomeAlleles_annotationfix.v1.0.tsv.gz
```

Post-analysis tools are included at $CONDA_PREFIX/share/ctyper/tools/ 

if you cannot find it, you may obtain it from GitHub folders here. 


Let us get start to genotype an example sample NA12718

```bash
wget "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram"
wget "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239480/NA12718.final.cram.crai"
```

$LIBRARY_PATH is your path of LD_LIBRARY, if you are using conda to install HTSlib or samtools, then it should be usually at /home/$user_name/miniconda3/lib/, otherwise have a try at /usr/local/lib/

```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH
Ctyper/ctyper -i NA12718.final.cram -m HprcCpcHgsvc_cmr_matrix.txt -o ctyper.out > log.txt &
```

Alternatively, for only a few target gene, for example HLAs
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH
Ctyper/ctyper -i NA12718.final.cram -m HprcCpcHgsvc_cmr_matrix.txt -g HLA* -B TargetRegions1KG.bed -o ctyper.out > log.txt &
```

After finishing genotyping, let us start to visualize and analyze, making sure you have python3 installed. 

first, change the raw output to a text table. The headers are: 
allele name, allele type, HG38 gene(s), transcriptname(exon_number):transcript_genename:similarity, type category(ref: in same type of HG38 gene, alt: different orthology type, dup: duplicated paralog, novel: diverged paralog), pangenome assembly location, reference liftover location and alignment CIGAR string. 
```bash
python Ctyper/tools/Annotation/SampleAnnotate.py -i ctyper.out -a PangenomeAlleles_annotationfix.v1.0.tsv -o genotype.txt
```

If you want to have fasta sequences. 
```bash
wget "https://zenodo.org/records/16340156/files/Allalters.fa"
python getFASTA.py -i genotype.txt -r HG38_main.fa,Allalters.fa,CHM13.fa -a PangenomeAlleles_annotationfix.tsv.gz -o output.fa
```

If you want to have nomenclature of HLA, CYP2D6, and KIR. 
```bash
python Ctyper/tools/Annotation/Nomenclature/GenotypetoNomenclature.py -i  genotype.txt -a Ctyper/data/all_nomenclature.txt > output__nomenclature.txt
```

If you want to have vcf file of HLA, CYP2D6, and KIR. 
```bash
python Ctyper/tools/Annotation/Nomenclature/GenotypetoVCF.py -i genotype.txt -a PangenomeAlleles_annotationfix.tsv.gz -o output.vcf
```

If you want to have mutant plot visualizations, for example SMN genes. 
```bash
python Ctyper/tools/Plot/typemutant.py -i genotype.txt (optional: -G SMN.gff3) -a PangenomeAlleles_annotationfix.tsv.gz -g SMN -o SMN.png
```

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
You need RAM > 8G to run all genes at once. 
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
6. **Jellyfish files** (`*.jy`, `*.txt`).
  
<!-- Installation -->
# Installation from conda

   ```bash
   conda create -n ctyper_env  -c conda-forge -c bioconda ctyper
   conda activate ctyper_env
   ```


# Installation from source

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

Required:
- **Database:**
  - `-m <file>`: Path to the matrix database (requires <file>.index). If not provided, runs in dry-run mode to estimate NGS read depth only.
    
- **Inputs:**
  - `-i <file>`: Path to an individual input file.
  - `-I <file>`: Path to a file listing multiple input files (one per line).

- **Outputs:**
  - `-o <file>`: Path to the individual output file. The output will be appended.
  - `-O <file>`: Path to a file listing multiple output files (one per input file, corresponding by line number).

Genotyping targets Options:
- **Gene Targets**:
  - `-g <string>`: Target gene name, prefix (ending with '\*', e.g. 'HLA\*'), or matrix (starting with '#', e.g. '#SMN_group1'). Can be specified multiple times.
  - `-G <file>`: Path to a file listing multiple genes (one per line).
- **Region specified**:
  - `-B <file>`: BED file to restrict region analysis. Make sure its 2nd name field matches your reference genome MD5 (use md5sum $reference). With -g/-G, only BED entries with names matching genes/matrices are used. One made from profiling run on EBI/GRCh38_full_analysis_set_plus_decoy_hla.fa is included in [github/$Ctyper/](https://github.com/Walfred-MA/Ctyper/data/). 
  - `-r <chr:start-end>`: Add a specific region for analysis. Can be specified multiple times. Regions will be merged if provided multiple times, can work with -B.
  - `-r gene`: Special value for -r; add regions from matrix database (must be combined with -g/-G). Not recommended if profile BED available or running global mode.
  - `-r Unmap`: Special value for -r; Force include all unmapped reads. No need to specify.
  - `-r HLA`: Special value for -r; Force include reads on all HLA decoys. No need to specify.

Running options:
- **Multithreading:**
  - `-n <int>`: Number of threads to run samples in parallel (default: 1). Parallel per-sample is more efficient for large cohorts, especially on slow disks.
  - `-N <int>`: Number of threads per sample. Due to file I/O bottleneck, suggest 1-4 for HDD. Parallelism within sample is memory-friendly (default: 1).

NGS file information:
- **Coverage Information**:
  - `-d <float>`: Fixed 31-mer depth value (incompatible with -D and -b). 31-mer depth = (1 - 30/read_length) × sequencing_depth. For 150 bp reads, this = 0.8 × sequencing_depth.
  - `-D <string>`: File of depth values (one per input, line-by-line, incompatible with -b and -d).
  - `-b <string>`: Background k-mer file for NGS coverage estimation (incompatible with -d/-D). In target runs, randomly generated 1M regions are used. Default: <matrix>.bgd
  - `-c <0/1>`: Enable NGS k-mer coverage bias correction (default: 1).
    
- **CRAM/BAM/SAM files:**
  - `-T <file>`: Reference FASTA file for reading CRAM files (default: use REF_CACHE and REF_PATH environment variables).

Supplementary profiling run to generate target regions:
- **Profiling Run options:**
  - `-p <file>`: Input aligned NGS file for profiling.
  - `-P <file>`: File listing multiple aligned NGS files for profiling (one per line). Can be used with both -O and -o; individual results go to -O paths, summary saved to -o.
---

# Results Annotation

Here shows the commands to make the genotyping results interpretable:

1. **Summary the genotyping results into a table with annotation information**

   ```bash
   python tools/Annotation/SampleAnnotate.py -i $ctyper_outputs.txt -a PangenomeAlleles_annotationfix.tsv > genotype.txt
   ```

To understand the output file, genotype.txt is a tsv table file. Its headers are:

| Column Name     | Description                           |
|-----------------|---------------------------------------|
| allele_name       | Unique identifier for the pangenone allele, format is $prefix_$groupindex_$sample_$haplotype_$index   |
| clade_name        | which clade (type of an allele, alias to HLA nameclature) this pangenone allele belongs to |
| transcript      | The exonic DNA it contains, format is $transcript_id:$gene_name:$similarity, if have multiple transcripts, then separated by semicolon, and if only contains a part of transcript, then it would be $transcript_id($start_exon_index-$end_exon_index):$gene_name:$similarity |
| classfication        | Ref: the alleles very similar to reference genes, Alt: alternative version of reference genes, Dup: duplicated paralogs, Novel: paralogs with novel sequences |
| pangenome_location | the location of the assembly, format: $sample#$haplotype#$contig:$start-$end$strand |
| liftover_location | The liftover location of the allele, format: $chromosome:$start-$end$strand:$alignment_CIGAR  |

Example row:
A4GALT_group1_GW00056_h1_227    A4GALT_group1_9 A4GALT  ENST00000401850:A4GALT:99.68;   Alt     GW00056#1#GWHBKGU00000010:748831-766661-        chr22:42689120-42706939-:2058M2I5D303M1X1463M4I3116M1X261M1X279M1X82M1X448M1X459M1X124M1X85M2D225M1X1240M1X650M9I794M1X144M18D20I929M1X399M1X793M1X83M1X709M1X182M1X296M1X136M1X650M1X96M1X59M3I291M1I95M1X276M1X245M1X263M17D16I36M1X184M1X299M


2. **(optional) Obtain public nomenclatures for important genes,including HLA, CYP2D6, and KIR**

   ```bash
   python tools/Annotation/Nomenclature/GenotypetoNomenclature.py -i genotype.txt -a data/all_nomenclature.txt > nomenclature.txt
   ```

nomenclature.txt is also a table file with two columns: first column the name the genotyped pangenome-allele, and the second is the public nameclature this pangenome-allele contains. 


3. **(optional) Convert genotyping results into VCF file**

   ```bash
   python tools/Annotation/VCF/GenotypetoVCF.py -i genotype.txt -o genotype.vcf
   ```

4. **(optional) Convert genotyping results into fasta file**

   ```bash
   python getFASTA.py -i genotype.txt -r HG38_main.fa,Allalters.fa,CHM13.fa -a PangenomeAlleles_annotationfix.tsv.gz -o output.fa
   ```


---

# Results Visualization

Visualization is performed on a **gene-by-gene** basis (not genome-wide).

For example, to visualize the gene **SMN**:
		
1. **Visualize the results:**

	```bash
	python typemutant.py -a PangenomeAlleles_annotationfix.tsv.gz  -g SMN -i genotype.txt -o output.png
	```

**Optional:** To visualize the GENCODE genes:
	
1. **Obtain the GENCODE annotation:**

	```bash
	grep "gene_name=SMN" genecode.gff3 > SMN.gff3
	```
	
2. **Run the visualization with GENCODE annotation:**

	```bash
	python typemutant.py -a PangenomeAlleles_annotationfix.tsv.gz -g SMN -G SMN.gff3 -n genotype.txt -o output.png
	```

**Optional:** Substract the visualization with certain genes:
	
	```bash
	zcat PangenomeAlleles_annotationfix.tsv.gz | grep "SMN1:" > Substract_annotationfix.tsv
	```
	
	```bash
	python typemutant.py -a Substract_annotationfix.tsv -g SMN -G SMN.gff3 -n genotype.txt -o output.png
	```
example output of SMN mutant map: 

Each row is each allele, and each vertically location is the consensus in multiple sequence alignments.  
The gap is the deletion/missing sequences on each allele and each black dots is each variants on HG38.  
Each color represents each allele-type we defined.  
The location of genotyped alleles are slightly bolded in color (may need to zoom in sometimes) and highlighted with a red dot on the left.  
(optional) The genecode gene and exons are highlight at the top.

![Alt Text](images/exampleSMN.png)

## Cohort Analysis

There are two scripts in the `tools/Cohort` folder that work together for cohort analysis.

1. **Download the allele-type nomenclature table:**

   ```
   wget "https://zenodo.org/records/16340156/files/PangenomeAlleles_typefix.tsv"
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
   python SummaryAlleles.py -f $results_folder/ -t PangenomeAlleles_typefix.tsv -o cohort_results.vcf
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






