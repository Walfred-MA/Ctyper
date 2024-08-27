#!/usr/bin/env python3

import pandas as pd
import collections as cl
import os
import argparse
import math
import numpy as np
import multiprocessing as mul

def runjob(inputfile, output, name_to_allele, allele_groups, multi):
	
	groupname = ""  
	allele_counts = cl.defaultdict(int)
	with open(inputfile, mode  ='r') as f:
		
		for line in f:
			
			if len(line) == 0 or  ( "result:" not in line and line[0] != '#' ):
				continue
			
			if line[0] == '>':
				groupname = line.split()[0][2:]
				continue
			
			line = line.split("result:")[1] 
			genes = [x for x in line.strip().split(",")[:-1]]
			
			if len(line.strip()) == 0:
				alleles = allele_groups[groupname]
				
				for allele in alleles:
					allele_counts[allele] = -1
					
			for gene in genes:
				
				#genename, genecount = gene.split(":")
				
				if "(aux)" in gene:
					print(gene)
					
				genename = gene
				
				if genename in multi:
					alleles = multi[genename]
				else:
					alleles = [ name_to_allele.get(genename,"NA") ]
					
				for allele in alleles:
					if allele == "NA":
						continue
					if "(aux)" in gene:
						allele_counts[allele] = -1
					elif allele_counts[allele] != -1:
						allele_counts[allele] += float(1.0) 
	with open(output, mode = 'w') as f:
		
		for allele, counts in allele_counts.items():
			
			f.write("{}\t{}\n".format(allele, counts))
			
			
def loadtable(tablefile,fixchunk):
	
	multi = dict()
	
	allele_groups = cl.defaultdict(list)
	name_to_allele = cl.defaultdict(lambda:"NA")
	with open(tablefile, mode = 'r') as f:
		for line in f:
			line = line.split()
			names = line[-1].split(",")
			
			if fixchunk:
				for name in names:
					if ":" in name:
						multi[name.split(":")[0]] = name.split(":")[1].split(";")
						
			names = [x for x in names if "@" not in x and ":" not in x]
			alignto = line[0]
			
			alignto_group = alignto
			allele_groups[alignto_group].append(alignto)
			for name in names:
				name_to_allele[name] = alignto
				
				
				
	return name_to_allele, allele_groups, multi

def main(args):
	
	manager = mul.Manager()
	
	name_to_allele, allele_groups,multi = loadtable(args.table, args.fixchunk)
	
	
	
	if len(args.input):
		
		runjob(args.input, args.output, name_to_allele, allele_groups,multi)
		
	else:
		name_to_allele[''] = ''
		allele_groups[''] = ''  
		name_to_allele=  manager.dict(dict(name_to_allele))
		allele_groups = manager.dict(dict(allele_groups))
		
		allfiles = [args.folder+"/"+x for x in os.listdir(args.folder) if x.endswith(".txt")]
		p=mul.Pool(processes=args.threads)
		
		
		for thefile in allfiles:
			#result = p.apply_async(runjob, (thefile, thefile + "_allelecount.out", name_to_allele, allele_groups, multi))
			runjob(thefile, thefile + "_allelecount.out", name_to_allele, allele_groups,multi)
			
		p.close()
		p.join()
		
		
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program determines allele specific gene counts from ctyper result")
	parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, default = "")
	parser.add_argument("-f", "--folder", help="path to input data file", dest="folder", type=str, default = "")
	parser.add_argument("-t", "--table", help="path to annotation table file", dest="table", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output", type=str, default = "")
	parser.add_argument("-n", "--threads", help="number if threads", dest="threads", type=int, default = 4)
	parser.add_argument("-c", "--fixchunk", help="if fix inconsistent truncating", dest="fixchunk", type=int, default = 1)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
if __name__ == "__main__":
	run()
	