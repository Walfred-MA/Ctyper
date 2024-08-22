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
    <a href="https://github.com/Walfred-MA/Ctyper/issues/new?labels=enhancement&template=feature-request---.md">Request Feature</a>
  </p>
</div>



<!-- ABOUT THE PROJECT -->
## About ctyper

Ctype is a command line tool to genotype complex CNVs (copy number variation) in NGS (Next-Generation Sequencing) data using pangenome database.

The results will be represented as pangenome alleles, which is defined as genic segments with locally phased variants that are combinatorially heritable mostly range in 15-50 kb, about the size of "Haplotype blocks" or "LD-blocks". 

The pangenome alleles are further classified as pangenome allele-types, which are highly similar subgroups, and can be used to represent complex genetic variations such as structrual variations, gene conversion, duplication, translocation and etc. 

Ctyper is highly efficient, accurate and visualizable, thus allows high resolution large cohort association studies on complex CNV genes as well as complex genetic variations. 

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<!-- GETTING STARTED -->
## Getting Started

This repository contains:
1. ctyper itself (src)

A c++ program use to perform genotyping on NGS data

2. a tool to visualize the ctyper's results (tools/Plot)

A python tool use to visualize genotyping results using multiple sequence alignments

3. a tool to help summary cohort study results (tools/Cohort)

A python tool use to summary results from all samples and include annotation information.

4. a pangenome allele database ctyper relies on (link in data)

The database used by ctyper for genotyping.

5. the annotation of pangenome allele-type to help population analysis (link in data)

The database with annotation information used mostly by cohort study tool.

6. A more detailed annotation of all individual pangenome alleles to help individual sample study (link in data)

The database with annotation information used in visualization.

7. simple test cases (tests)





<!-- USAGE Prerequisites -->
### Prerequisites

Ctyper only officially support Linux environment.

It requires:
1. gcc-8
2. eigen (https://gitlab.com/libeigen/eigen)
3. htslib (https://github.com/samtools/htslib), you can also install the whole samtools instead (https://github.com/samtools/).

<!-- USAGE Installation -->
### Inputs

Ctyper takes five types of files as input:

1. CRAM file (*.cram), needed to be indexed  (recommended due to less I/O intensive).
2. BAM file (*.bam), needed to be indexed.
3. SAM file (*.sam).
4. FASTQ file (*.fastq).
5. Fasta file (*.fa, *.fasta).

<!-- USAGE Installation -->
### Installation

1. Install all Prerequisites.
2. go to src/ directory.

cd src/

3. Compile ctyper. 

make 

4. Move the ctyper to install folder

mv ctyper $folder_install

6. Download the pangenome allele database file (we call it $Database)
   
7. Download the pangenome allele-type annotation file
   
8. (optional) Download the full pangenome alleles annotation file
   
9. ready to go

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

Ctyper needs:
1. inputfile(s),
2. corresponding outputfile(s),
3. the pangenome allele database file,
4. either directly proving sequencing coverage information or proving background kmers to determine sequencing coverage.  



Ctyper takes two kinds of input formats:

1. genotyping single file.

ctyper -i $Inputfile -m $Database -o $Outputfile -b background.list -c 1
or 
ctyper -i $Inputfile -m $Database -o $Outputfile -d $sequencing_coverage -c 1


2. genotyping a cohort of files.

ctyper -I $AllInputs -m $Database -o $AllOutputs -b background.list -c $threads
or 
ctyper -I $AllInputs -m $Database -o $AllOutputs -D $All_sequencing_coverages -c $threads

AllInputs is a text where each line is the path of each input file.
AllOutputs is a text where each line is the output file of the correponding input file (the input file with the same row number).
All_sequencing_coverages is a text where each line is the sequencing coverage information for the correponding input file (the input file with the same row number).


## Parameters



## Addtional tools


## Demo

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Walfred MA - wangfeim@usc.edu
Mark Chaisson - mchaisso@usc.edu
<p align="right">(<a href="#readme-top">back to top</a>)</p>






