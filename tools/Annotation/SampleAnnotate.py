#!/usr/bin/env python3

import argparse
import parser
import collections as cl

def summary(inputfile, annofile, outputfile, iffilter):
	
	results = set()
	with open(inputfile, mode = 'r') as f:
		
		for line in f:
			
			if line.startswith("result: "):
				
				results.update(line[8:].split(",")[:-1])
	
	with open(outputfile, mode = 'w') as w:	
		
		with open(annofile, mode = 'r') as f:
			
			for line in f:
				
				line = line.strip().split('\t')
				
				if iffilter and line[6] != "Exon":
					
					continue
				
				if line[0] in results:
					
					line = line[:3] + line[4:6] + [line[-1]]
					if ":" in line[2]:
						line[2] = line[2].split(":")[0]
						
					if ":" in line[-1]:
						line[-1] = ":".join(line[-1].split(":")[1:])
					
					w.write("\t".join(line)+"\n")
			



def main(args):
	
	summary(args.input,  args.anno, args.output, args.filter)
	
def run():
	"""
			Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program to summary ctyper genotyping results for a NGS data")
	parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
	parser.add_argument("-a", "--anno", help="path to input data file",dest="anno", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
	parser.add_argument("-f", "--filter", help="if filter decoys and introns", dest="filter",type=int, default = 1)

	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()