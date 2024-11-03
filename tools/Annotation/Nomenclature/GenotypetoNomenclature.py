#!/usr/bin/env python3

import argparse
import parser
import collections as cl


def nomenclature(inputfile, annofile ):
	
	nametonomen = dict()
	with open(annofile, mode = 'r') as f:
		
		for line in f:
			
			line = line.split()
			prefix = line[0].split('_')[0]
			
			nomenclatures = ",".join([x for x in line[-1].split(",") if x.startswith(prefix) or x.startswith("*") or x.startswith("None")])
			
			
			if len(nomenclatures):
				nametonomen[line[0]] = nomenclatures
	
	with open(inputfile, mode = 'r') as f:
		
		for line_ in f:
							
			allele = line_.strip().split()[0]
			if allele in nametonomen:
				print(allele, nametonomen[allele])
	


def main(args):
	
	nomenclature(args.input, args.anno)
	
def run():
	"""
			Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="the program to change ctyper genotype to public nomenclature, currently including HLA, CYP2D6 and KIR")
	parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
	parser.add_argument("-a", "--nomenclature", help="path to nomenclature annotation file",dest="anno", type=str, required=True)
	

	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
