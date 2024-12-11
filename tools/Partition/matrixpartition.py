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


def partition(inputfile, partnum = 2):
	
	
	filepaths = []
	outfiles = []
	for index in range(partnum):
		
		filepath = inputfile+"_part"+str(index+1)
		
		filepaths.append(filepath)
		outfiles.append(open(filepath, mode = 'w'))
		
	
	matrix = ""
	index = 0
	with open(inputfile, mode = 'r') as f:
		
		line = f.readline()
		kmer_num = 0
		eles_num = 0
		while line:
			if len(line):
				if line[0] == '#':
					
					currentfile = outfiles[index%partnum]
					currentfile.write(matrix)
					index += 1
					matrix = line
					
				else:
					matrix += line
					
			line = f.readline()
			
		currentfile.write(matrix)
		
	for file in outfiles:
		file.close()
		
	for filepath in filepaths:
		indexfile(filepath)
		
	
def main(args):

	partition(args.input, partnum = args.partnum)

			
			
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program run partition")
	parser.add_argument("-i", "--input", help="path to input data file", dest="input", type=str, required=True)
	parser.add_argument("-p", "--part", help="number of partitions", dest="partnum", type=int, required=True)
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
