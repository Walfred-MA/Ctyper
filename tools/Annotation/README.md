This is the readme file for SampleAnnotate.py, and mainly decribed its output table. 

1.Requirement: python >3.7

2.Running SampleAnnotate.py:

bash
python SampleAnnotate.py -i $Inputfile -a $Sample_Annotation > result.txt

3.Output description

Output will be a text TSV table. Each row represents the information of each genotyped pangenome-alleles. 

The columns of this table are: 

Name 
Allele type 
Overlap genes
Mapped transcripts
relationship to reference genes: 
  -Ref: same allele-type; 
  -Alt: different alleles-type; 
  -Dup: alleles-types of similar paralog; 
  -Novel: alleles-types of diverged paralog
Location in the pangenome assemblies:
  -format: contig:start-end[strand]
Liftover on the reference:
  -format: chromsome:start-end[strand]:CIGAR-extended
  -CIGAR-extended: extended version of CIGAR strings

4. Description of the CIGAR-extended used to encode liftover pairwise alignment.

CIGAR-extended is an extended version of CIGAR. The extension allow it to include all variant information in its text. 

For each variant (X/I/D), the variant information will be appended to the end of their cigar segment. For example 6D will extended to 6DAAATTT, means AAATTT has been deleted. 
For mismatch, both orginal sequences and changed sequences will be appended in sequential. 









