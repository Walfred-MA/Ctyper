#!/usr/bin/env python3

import os
import argparse
import collections as cl


def indexfile(inputfile):
	
	names = []
	locations = []
	kmer_nums = []
	eles_nums = []
	with open(inputfile, mode = 'r') as f:
		
		line = f.readline()
		kmer_num = 0
		eles_num = 0
		while line:
			if len(line):
				if line[0] == '#':
					locations.append(f.tell() - len(line))
					names.append(line.split()[0])
					kmer_nums.append(kmer_num)
					eles_nums.append(eles_num)
					kmer_num = 0 
					eles_num = 0
				elif line[0] == "&":
					kmer_num += 1
					eles_num += line.count(",")+6
					
					
			line = f.readline()
			
		locations.append(f.tell())
		kmer_nums.append(kmer_num)
		eles_nums.append(eles_num)
		
	with open(inputfile+".index", mode = 'w') as f:
		
		last_location = 0
		
		for name, location,kmer_num, eles_num in zip(names, locations[1:], kmer_nums[1:], eles_nums[1:]):
			
			f.write("{}\t{}\t{}\t{}\t{}\n".format(name, last_location, location, kmer_num, eles_num))
			last_location = location


def partition(args):
	
	
	filepaths = []
	outfiles = []
	for index in range(args.partnum):
		
		filepath = args.input+"_part"+str(index+1)
		
		filepaths.append(filepath)
		outfiles.append(open(filepath, mode = 'w'))
	
	genelist = set(args.genes.split(","))
	if len(args.genelist):
		
		with open(args.genelist, mode = 'r') as f:
			
			genelist.update([x.strip() for x in f.readlines()])
	
	skip = 0
	matrix = ""
	index = 0
	with open(args.input, mode = 'r') as f:
		
		line = f.readline()
		kmer_num = 0
		eles_num = 0
		while line:
			if len(line):
				if line[0] == '#':
					
					currentfile = outfiles[index%args.partnum]
					currentfile.write(matrix)
					
					if len(genelist) and line.split()[0][1:] not in genelist:
						del matrix
						matrix = ""
						skip = 1
							
					else:
						skip = 0
						index += 1
						del matrix
						matrix = line
					
				elif skip == 0:
					kmer_num += 1
					matrix += line
					
			line = f.readline()
		
		currentfile.write(matrix)
		
	for file in outfiles:
		file.close()
		
	for filepath in filepaths:
		indexfile(filepath)
		
	
def main(args):

	partition(args)

			
			
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program run partition")
	parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
	parser.add_argument("-p", "--part", help="number of partitions", dest="partnum", type=int, default = 1)
	parser.add_argument("-g", "--genes", help="distracted gene matrice(s), separated by comma", dest="genes", type=str, default = "")
	parser.add_argument("-G", "--genelist", help="file path of distracted gene matrix list, separated by line", dest="genelist", type=str, default = "")
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
