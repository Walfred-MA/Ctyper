#!/usr/bin/env python3

import argparse
import gzip
import re
import collections as cl

NOTE_PATTERN = re.compile(r"\(([^)]+)\)")

def interval_intersection(X, Y, frac_threshold=0.5):

	i = j = 0
	result = []
	x_significant = set()
	X = sorted(X)
	Y = sorted(Y)

	while i < len(X) and j < len(Y):
		a_start, a_end = X[i][:2]
		b_start, b_end = Y[j][:2]
	
		# Compute overlap
		start = max(a_start, b_start)
		end = min(a_end, b_end)
	
		if start < end:  # Non-zero intersection
			result.append([start, end])

			# Check fraction of X[i] that is covered
			x_len = a_end - a_start
			overlap_len = end - start
			if x_len > 0 and (overlap_len / x_len) >= frac_threshold:
					x_significant.add(i)
					
		# Advance the interval that ends first
		if a_end < b_end:
				i += 1
		else:
				j += 1
				
	return result, x_significant


def compress_ranges(numbers):

	if len(numbers) == 0:
		return ""

	
	compressed = []
	start = prev = numbers[0]
	
	for n in numbers[1:]:
		if n == prev + 1:
			prev = n
		else:
			if start == prev:
				compressed.append(f"{start}")
			else:
				compressed.append(f"{start}-{prev}")
			start = prev = n
			
	# Add final group
	if start == prev:
		compressed.append(f"{start}")
	else:
		compressed.append(f"{start}-{prev}")
		
	return compressed

def parse_attributes(attr_str):
	"""Parse GFF3 attributes column into a dict."""
	attrs = {}
	for entry in attr_str.strip().split(";"):
		if "=" in entry:
			key, value = entry.split("=", 1)
			attrs[key] = value
	return attrs

def readgff3(inf, records):
	
	lastchrom = ""
	lastexon = 0
	transsize = cl.defaultdict(int)
	genetoanno = cl.defaultdict(list)
	transtogenes = cl.defaultdict(str)
	transtoexons = cl.defaultdict(list)
	exontogenenum = cl.defaultdict(int)
	genecounter = cl.defaultdict(int)
	manelist = cl.defaultdict(lambda: [[],[]])
	with open(inf, "r", encoding="utf-8") as f:
		
		for line in f:
			if line.startswith("#"):
				continue
			cols = line.strip().split("\t")
			if len(cols) != 9 or cols[1] != "HAVANA":
				continue
			
			chrom, source, feature, start, end, score, strand, phase, attributes = cols
			start, end = int(start), int(end)
						
			attrs = parse_attributes(attributes)
			
			if feature == "gene" and "ID" in attrs and "gene_type" in attrs:
				gene_id = attrs["gene_id"].split(".")[0]
				anno = attrs["gene_type"]
				name = attrs["gene_name"] if "gene_name" in attrs else ""
				
				genetoanno[gene_id] = [name, anno]
				lastexon = [0,0]
				
			if feature == "transcript" :
				gene_id = attrs["gene_id"].split(".")[0]
				enst_id = attrs["transcript_id"].split(".")[0]
				
				if "MANE_Select" in attrs:
					manelist[gene_id][0].append(enst_id)
				else:
					manelist[gene_id][1].append(enst_id)
			
			if feature == "exon" and "ID" in attrs and "transcript_name" in attrs and "gene_id" in attrs:
				ense_id = attrs["exon_id"].split(".")[0]
				enst_id = attrs["transcript_id"].split(".")[0]
				ensg_id = attrs["gene_id"].split(".")[0]
																								
				transtogenes[enst_id] = ensg_id
				transtoexons[enst_id].append((start, end, ense_id))
				transsize[enst_id] += end - start
				
	for gene_id, thelist in manelist.items():
		
		mane = sorted(thelist[0], key = lambda x: transsize[x], reverse= 1)+sorted(thelist[1], key = lambda x: transsize[x], reverse= 1)
		
		manelist[gene_id] = [mane[0], transtoexons[mane[0]]]
	
	return genetoanno, transtoexons, transtogenes, manelist

def readtable(inf):
	
	records = cl.defaultdict(list)
	with open(inf, "r", encoding="utf-8") as fin:
		for ln, line in enumerate(fin, 1):
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			parts = line.split()
			if len(parts) < 4:
				# Not enough columns; skip gracefully
				continue
			
			c1 = parts[0]
			c4 = [x for x in parts[4].split(";") if len(x)]  # the 4th column (0-based index 3)
		
			# Split 4th column on ';' and parse each non-empty token
			for token in c4:
				
				field1, field2, field3 = token.split(":")
				note = NOTE_PATTERN.search(field1)
				note = f"({note.group(1)})" if note else ""
				field1 = field1.split("(")[0]
				records[field1].append((note, field2, float(field3), c1))

	return records

def summary(records, genetoanno, transtoexons, transtogenes, manelist):
	
	gene_records = cl.defaultdict(list)
	gene_cn = cl.defaultdict(float)
	for trans, data in records.items():
		
		if trans not in transtogenes:
			continue
		
		genename = transtogenes[trans]
		
		
		
		exons = transtoexons[trans]
		for (note, field2, field3, c1) in data:
			
			if float(field3) < 0.95:
				continue

			if note == "":
				mane = trans
				mane_exons = exons
				exon_start = 1
				exon_end = len(exons) 
			else:
				mane, mane_exons = manelist[genename]
				exon_start, exon_end = [int(x) for x in note[1:-1].split("-")]
			
			genesize = sum([x[1]-x[0] for x in mane_exons])
			
			if exon_end > len(exons):
				continue
			exon_range = range(exon_start-1, exon_end)
			exon_coordi = [exons[i] for i in exon_range]
			exon_overlaps, exon_index = interval_intersection(mane_exons, exons)
			
			
			if len(exon_index) == 0:
				continue
			
			exon_index = '_'.join(compress_ranges(sorted(list(exon_index))))
			totalsize = sum([x[1]-x[0] for x in exon_overlaps])
			cn = totalsize/max(1,genesize)
			gene_cn[mane] += cn
			
				
			gene_records[mane].append((c1, genesize,len(mane_exons), f"Exon_Index:{exon_index}" ,field3))
		
			
	
	return gene_cn,gene_records
		
		
		

def run(args):
	
	records = readtable(args.input)
	
	genetoanno, transtoexons, transtogenes, manelist = readgff3(args.gene, records)
	
	gene_cn,gene_records = summary(records, genetoanno, transtoexons, transtogenes, manelist)
	# Write results
	
	results = []
	
	for trans, cn in gene_cn.items():
		
		genename = transtogenes[trans]
		genesize,exonnum = gene_records[trans][0][1:3]
		
		record = sorted(gene_records[trans], key = lambda x: (int(x[3].split(":")[1].split('-')[0]),int(x[3].split(":")[1].split('-')[-1]),-x[2],x[0]))
		
		record_str = ";".join([f"{x[0]}:{x[3]}:{x[4]}" for x in record])
		results.append(f"{trans}\t{genename}\t{genetoanno[genename][0]}\t{genetoanno[genename][1]}\t{genesize}\t{exonnum}\t{cn}\t{record_str}\n")
	
	results = sorted(results, key = lambda x: (x.split('\t')[2], x))

	
	with open(args.output, "w", encoding="utf-8") as fout:
		fout.write("genename\tid\ttype\ttotal_size\ttotal_exons\taggregate_copy_number\talleles\n")
		
		for line in results:
			fout.write(line)
			
def main():
	parser = argparse.ArgumentParser(
		description="Getting aggregate copy number calling."
	)
	parser.add_argument("-i", "--input", required=True, help="Input text file")
	parser.add_argument("-g", "--gene", required=True, help="Input genecode database")
	parser.add_argument("-o", "--output", required=True, help="Output TSV file")
	args = parser.parse_args()
	run(args)
	
if __name__ == "__main__":
	main()
