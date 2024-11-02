<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
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
    <a href="https://github.com/Walfred-MA/Ctyper/issues/new?labels=enhancement&template=suggest-a-gene---.md">Suggest include new gene</a>
  </p>
</div>


<!-- ABOUT THE PROJECT -->
## About ctyper

Ctype is a command line tool to perform copy number sensitive versatile genotyping for NGS (Next-Generation Sequencing) data using pangenome database. It is designed to work on complex CNV (copy number variation) genes, but can also work as genotyping, local phasing or SV-calling tools for other genes. 

The results will be represented as pangenome alleles, which is defined as genic segments with locally phased variants that are combinatorially heritable mostly range in 15-50 kb, about the size of "Haplotype blocks" or "LD-blocks". The pangenome alleles are further classified as allele-types among populations, which are highly similar subgroups, and can be used to represent complex genetic variations such as structrual variations, gene conversion, duplication, translocation and etc. 

The output format is in text format, but may be converted into VCF format (however, not recommended). 

Ctyper is highly efficient, accurate and visualizable, thus allows high resolution large cohort association studies on complex CNV genes as well as complex genetic variations. 

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- GETTING STARTED -->
This repository includes the following components:
	
ctyper (src)  
A C++ program for performing genotyping on NGS data.  

Visualization Tool (tools/Plot)  
A Python tool for visualizing genotyping results using multiple sequence alignments.  

Cohort Analysis Tool (tools/Cohort)  
A Python tool for summarizing results across all samples in a cohort study, including annotation information.  

Pangenome Allele Database (link in data)  
A database required by ctyper for genotyping.  

Population Analysis Annotation Database (link in data)  
A database containing annotation information to support population-level analysis.  

Individual Sample Annotation Database (link in data)  
A comprehensive database with detailed annotations for individual pangenome alleles, supporting individual sample studies and visualization.  

Test Cases (tests)  
Simple test cases for validating the tools and pipeline.  




<!-- Prerequisites -->
## Prerequisites

Ctyper is officially supported only in Linux environment.

It requires:
1. gcc-8
2. eigen (https://gitlab.com/libeigen/eigen)
3. htslib (https://github.com/samtools/htslib), you can also install the whole samtools instead (https://github.com/samtools/).


For additional python tools, they require:

1. Python 3.7+
2. numpy
3. pandas (https://pandas.pydata.org/)
4. matplotlib (https://matplotlib.org/stable/install/index.html)


<!-- Inputs-->
## Inputs

Ctyper takes five types of files as input:

1. CRAM file (*.cram), needed to be indexed  (recommended due to less I/O intensive).
2. BAM file (*.bam), needed to be indexed.
3. SAM file (*.sam).
4. FASTQ file (*.fastq).
5. Fasta file (*.fa, *.fasta).

<!-- Installation -->
## Installation

1. Install all Prerequisites.
2. go to src/ directory.

$ cd src/

3. Compile ctyper. 

$ make 

4. Move the ctyper to install folder

$ mv ctyper $folder_install

6. Download the pangenome allele database file (we call it $Database) and its index file. 
   
7. (optional) Download the allele-type annotation file
   
8. (optional) Download the full alleles annotation file
   
9. ready to go

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

Ctyper needs:
1. inputfile(s),
2. corresponding outputfile(s),
3. the pangenome allele database file, needed to be indexed
4. either directly proving sequencing coverage information (when you do not have WGS data) or proving background kmers for ctyper to determine sequencing coverage (recommended).  



Ctyper takes two kinds of inputs:

1. genotyping single file.

ctyper -i $Inputfile -m $Database -o $Outputfile -b background.list -n 1

or 

ctyper -i $Inputfile -m $Database -o $Outputfile -d $sequencing_coverage -n 1


2. genotyping a cohort of files.

ctyper -I $AllInputs -m $Database -o $AllOutputs -b background.list -n $threads

or 

ctyper -I $AllInputs -m $Database -o $AllOutputs -D $All_sequencing_coverages -n $threads

AllInputs is a text where each line is the path of each input file.
AllOutputs is a text where each line is the output file of the correponding input file (the input file with the same row number).
All_sequencing_coverages is a text where each line is the sequencing coverage information for the correponding input file (the input file with the same row number).


## Parameters

Currently support parameters:

Help:
  -h print for help

Inputs:

  -i string, the path of individual input file
  
  -I string, the path of a mega input file, where each line corresponds to each path of single files

Database:

  -m string, the database file use for genotyping
  
  -b string, the path of the background kmer list

Coverage information (does not co-exit with -b):

  -d float, the sequencing coverage of the input file
  
  -D string, the path of all sequencing coverages, where each line corresponds to each input file with the same row index.
  
Outputs:

  -o string, the path of individual output file
  
  -O string, the path of a mega output file, where each line corresponds to the output path for each input file with the same row index.

multhreads:

  -n int, number of thread use, default is 1

Bias correction:

  -c bool, if performs bias correction for illumina data


## Visualization

This visualization is gene by gene and not genome-wide. 

For example, if you interested in visualizing gene AMY1A:

First, Look for AMY1A to find out which gene group it in.

$ cat data/select_files.txt | grep -w "AMY1A" | cut -f2,3

newGeneCopies/AMY/AMY_partitions/AMY_group1_AMY1BOOOAMY1COOOAMY1A.fa	AMY1B,AMY1C,RNPC3,AMY1A,AMYP1,ACTG1P4,RP5-1108M17.5,AMY2B,AMY2A,


this shows it is in AMY_group1, together with other amylase genes. 


Second, download full alleles annotation, and distract the annotation of AMY_group1 from the full annotation table.

$ cat PangenomeAlleles_annotationfix.tsv  | grep "^AMY_group1_" > AMY_group1_annotationfix.tsv


Third, distract the AMY_group1 genotyping results:

$ cat genotype.txt  |  grep "^result: AMY_group1_" 

result: AMY_group1_GW00031_h1_556,AMY_group1_GW00051_h2_891,


this shows genotyping result is: AMY_group1_GW00031_h1_556,AMY_group1_GW00051_h2_891,


Last, visualization:

$ python typemutant.py -i AMY_group1_annotationfix.tsv -n "AMY_group1_GW00031_h1_556,AMY_group1_GW00051_h2_891,"  -o output.png

(optional)
if you also want to visualize the genecode annotation on MSA, you can do following steps:

first obtain genecode annotation:
$ cat genecode.gff3| grep "gene_name=AMY" > AMY.gff3


Then use it as the input, run:
$ python typemutant.py -i AMY_group1_annotationfix.tsv -g AMY.gff3 -n "AMY_group1_GW00031_h1_556,AMY_group1_GW00051_h2_891,"  -o output.png


## Cohort analysis

There are two scripts in the tools/Cohorts folder, which work together for Cohort analysis. 

first, download the allele-type annotation table: PangenomeAlleles_typefix.tsv. 

second, running CountAllele.py to each sample to get the allele-type results

$ for ($result in $results); do python CountAllele.py -i $result -t PangenomeAlleles_typefix.tsv -o "$result"_alleletype.out ; done

or if you want to run in parallel

$ python CountAllele.py -f $results_folder/ -t  PangenomeAlleles_typefix.tsv  -n $numthreads

Third, summarize results to a mega file and adding annotation

$ python summaryalleles.py -f $results_folder/ -t  PangenomeAlleles_typefix.tsv  -o cohort_results.tsv





<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Walfred MA - wangfeim@usc.edu
Mark Chaisson - mchaisso@usc.edu
<p align="right">(<a href="#readme-top">back to top</a>)</p>






