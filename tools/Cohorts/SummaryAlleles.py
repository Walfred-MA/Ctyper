#!/usr/bin/env python3

import os
import argparse
import collections as cl
import numpy as np
import multiprocessing as mul




def main(args):
	
	allfiles = [args.input+"/"+x for x in os.listdir(args.input) if x[-len("_allelecount.out"):] == "_allelecount.out"]
	
	num_samples = len(allfiles)
	
	allele_text = dict()
	allele_counts = cl.defaultdict(lambda : np.full(num_samples, 0.0, dtype = 'float' ))
	allele_sortidx = {}
	headers = []
	
	with open(args.typefile, mode = 'r') as f:
		for line in f:
			if len(line )==0:
				continue
			lines= line.split()
			if lines[0] == "Unknown" or lines[0] == "NA":
				continue
			allele_text[lines[0]] = line.strip()
			
			
			
	for index,tablefile in enumerate(allfiles):
		
		header = tablefile.split("/")[-1].split("_")[0]
		headers.append(header)
		
		with open(tablefile, mode = 'r') as f:
			
			for line in f:
				
				if len(line )==0:
					continue
				line= line.split()
				
				if line[0] == "Unknown" or line[0] == "NA":
					continue
				
				allele_counts[line[0]][index] += sum([int(float(x)+0.5) for x in line[1].split(",") if len(x)])
				allele_sortidx[line[0]] = index
				
	with open(args.output, mode = 'w') as f:
		
		f.write("1KGsamples\t"+"\t".join(headers)+"\n")
		
		all_allele = sorted(list(allele_counts.keys()), key = lambda x: int(x.split("_")[-1]))
		
		for allele in all_allele:
			
			counts = allele_text[allele]+"\t"+"\t".join(list(map(str, allele_counts[allele])))
			
			f.write(counts+"\n")
			
			
			
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program summarize cohort genotyping results")
	parser.add_argument("-i", "--input", help="path to input data file",
						dest="input", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output",
						type=str,required=True)
	parser.add_argument("-t", "--type", help="path to annotation file", dest="typefile",type=str,required=True)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
