# README for BatchRegionMerger.py

This README provides an overview of how to use `BatchRegionMerger.py` and describes its output format in detail.

## 1. Requirements
- Python version 3.7 or higher
- numpy
- numba

## 2. Running BatchRegionMerger.py

To run the script, use the following command:

    python BatchRegionMerger.py -i <InputFolder> -o <OutputFile> -n <Threads>

Replace `<InputFolder>` with the path to your input folder and `<OutputFile>` with the desired output file.
`<Threads>` is optional (default: 8).

Example:

    python BatchRegionMerger.py -i ./bed_files -o merged_regions.txt -n 8

## 3. Input Format

Each input file should be a tab-delimited file with four columns:
    chrom  start  end  group

Example:

    chr1    100     200     group1
    chr1    180     300     group1
    chr2    50      150     group2

## 4. Output Format

The output is a tab-delimited file with four columns:
    chrom   merged_start   merged_end   group

Example:

    chr1    100     300     group1
    chr2    50      150     group2

Overlapping or nearby regions (within 400 bp by default) are merged.

## 5. Description

BatchRegionMerger.py merges all regions for each (chrom, group) key across multiple files in parallel, using a specified number of threads.
