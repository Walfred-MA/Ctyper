#!/usr/bin/env python3

import re
import collections as cl
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
import os 
import argparse
import colorsys
import gzip


def adjust_lightness(color, amount=0.5):
	try:
		c = mcolors.cnames[color]
	except:
		c = color
	c = colorsys.rgb_to_hls(*mcolors.to_rgb(c))
	return colorsys.hls_to_rgb(max(0, min(1, amount * c[0])), max(0, min(1, amount * c[1])), max(0, min(1, amount * c[2])))


def getcigarsize(text):
	
	allcigars = sum([int(x[:-1]) for x in re.findall(r'\d+[MIX=]', text)])
	
	return allcigars

def maptoseg(allexons, segments):
	
	l = len( segments)
	
	allexons = [x if i %2 == 0 else x -1 for i,x in enumerate(allexons)]
	
	
	allposi =  segments + allexons 
	
	
	allposi_sortindex = sorted(range(len(allposi)), key = lambda i: allposi[i])
	
	overlaps = [[] for x in segments] 
	overlap = overlaps[0]
	seg_start = 0
	for rank, index in enumerate( allposi_sortindex):
		
		posi = allposi[index]
		if index < l :
			overlap = overlaps[index]
			seg_start = allposi[index]
		else:
			overlap.append(index)
			
	overlaps = overlaps[:-1]
	for i,overlap in enumerate(overlaps):
		
		overlap_ = [allposi[x]- allposi[i] if x % 2 == 0 else allposi[x]- allposi[i]+1 for x in overlap]
		
		if len(overlap) % 2 == 1:
			if (overlap[0] - l) % 2:
				overlaps[i] = [0] + overlap_
			else:
				overlaps[i] = overlap_ + [allposi[i+1] - allposi[i]]
		else:
			overlaps[i] = overlap_
			
	return overlaps

def overlapexons(allexons):
	
	if len(allexons) and type(allexons[0]) == type([]):
		l = len(allexons)
		allposi = sum(allexons, [])
	else:
		l = len(allexons)//2
		allposi = allexons
		
	allposi_sortindex = sorted(range(len(allposi)), key = lambda i: allposi[i])
	
	allgroups = []
	current_group = []
	last_endpoint = 0
	overlap_count = 0
	
	for rank, index in enumerate( allposi_sortindex):
		
		posi = allposi[index]
		if index %2 :
			overlap_count -= 1
			
			if overlap_count == 0:
				
				last_endpoint = allposi[index]
				current_group.append(posi)
		else:
			
			overlap_count += 1
			
			if overlap_count == 1 :
				
				allgroups.append(current_group)
				current_group = [posi]
				
	allgroups.append(current_group)
	
	return allgroups[1:]

def cigar_findrposi(cigar, find_qposes):
	
		l = len(find_qposes)
		if l == 0:
			return []
	
		curr_rposi = 0
		curr_qposi = 0
		new_rposi = 0
		new_qposi = 0
	
		curr_size = 0
	
		find_qpos_index = 0
		find_qpos = find_qposes[0]
		find_rposes = []
		find_qposes = list(find_qposes)+[-1]
	
	
		for char in cigar:
			
				if char <= '9':
					
						curr_size *= 10
						curr_size += ord(char) - ord('0')
						continue
			
				if char == '=' or char ==  'M':
					
						new_rposi = curr_rposi + curr_size
						new_qposi = curr_qposi + curr_size
					
				elif char == 'I':
					
						new_qposi = curr_qposi + curr_size
					
				elif char == 'D' or char == 'H':
					
						new_rposi = curr_rposi + curr_size
					
				elif char == 'X':
					
					
						new_rposi = curr_rposi + curr_size
						new_qposi = curr_qposi + curr_size
					
					
				elif char == 'S':
					
						if curr_qposi == 0:
								curr_qposi = curr_size
							
						curr_size = 0
						continue
			
				while find_qpos_index<l and find_qpos < new_qposi:
					
						if char == '=' or char ==  'M' or char =='X' or char == 'N':
								find_rposes.append(curr_rposi+(find_qpos-curr_qposi))
						else:
								find_rposes.append(curr_rposi)
							
						find_qpos_index += 1
						find_qpos = find_qposes[find_qpos_index]
					
					
				curr_size = 0
				curr_rposi = new_rposi
				curr_qposi = new_qposi
			
		while find_qpos_index<l:
			
				find_qpos_index += 1
				find_rposes.append(curr_rposi)
			
		return find_rposes

def adjust_lightness(color, amount=0.5):
	try:
		c = mcolors.cnames[color]
	except:
		c = color
	c = colorsys.rgb_to_hls(*mcolors.to_rgb(c))
	return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


def segfiltration(cigar, usednames = ""):
	
	cigars = re.findall('[<>][0-9_A-Za-z:]+',cigar)
	
	newcigars = []
	for cigarstr in cigars:
		
			thename,thecigar = cigarstr.split(":")
			if type(usednames) != type("") and thename[1:] not in usednames:
				
					continue
		
			thesize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=]', thecigar)]+[0])
		
			if thesize >0 and thesize < 50:
				
					fullsize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XHD]', thecigar)]+[0])
				
					thecigar = "{}H".format(fullsize )
				
			newcigars.append(thename+":"+thecigar)
		
	newcigars= "".join(newcigars)
	
	return newcigars

def getgenelocation(table, gff3file, used_cigars):
	
	segment_exonalign = cl.defaultdict(lambda : cl.defaultdict(list))
	gffrecord = cl.defaultdict(lambda : cl.defaultdict(list))
	if len(gff3file):
		
		with open(gff3file, mode = 'r') as f:
			for line in f:
				line = line.strip().split()
				
				if line[2] == 'exon':
					genename = [x for x in line[-1].split(";") if x.startswith("gene_name=")][0][10:]
					
					gffrecord[line[0]][genename].append([int(line[3]), int(line[4])])
					
		gffrecords = {thechr:{gene:overlapexons(coordi) for gene, coordi in record.items()} for thechr, record in gffrecord.items()}
		
		for row in table.values.tolist():
			
			locus = row[7]
			refchr = locus.split(":")[0]
			if refchr not in gffrecords:
				continue
			
			strd = locus[-1]
			
			thepath = row[-3]
			
			cigar = row[-2]
			segment_cigars = re.findall('[<>][0-9_A-Za-z:]+',cigar)
			segment_orders = re.findall('[<>][0-9_A-Za-z:]+',thepath)
			
			segment_cigars = {x.split(":")[0][1:]:x for x in segment_cigars if x.split(":")[0][1:].split("_")[0] in used_cigars}
			segment_orders = [x for x in segment_orders if x[1:].split("_")[0] in used_cigars]
			
			segment_cigars = [segment_cigars[x[1:]] for x in segment_orders]
			
			segment_alignsizes = [(x.split(":")[0], getcigarsize(x.split(":")[1])) for x in segment_cigars ]
			
			qcoordi = 0
			segment_qranges = [0]
			segment_names = []
			for segment in segment_alignsizes:
				segment_names.append(segment[0])
				segment_qranges.append(qcoordi + segment[1])
				qcoordi += segment[1]
			segment_qranges.append(qcoordi)
			
			start, end = locus.split(":")[1][:-1].split('-')
			start, end = int(start), int(end)
			
			for gene, gffrecord in gffrecords[refchr].items():
				
				exons = [x for x in gffrecord if max(x[1],end) - min(x[0],start) - (end - start) - (x[1]-x[0]) < 0 ] 
				
				if len(exons) == 0:
					continue
				
				if strd == '+' :
					exons = sum([ [ max(0, x[0] - start), x[1] - start ] for x in exons],[])
				else:
					exons = sum([ [ max(0, end - 1 - x[1]), end - x[0] ] for x in exons[::-1]],[])
					
				exons_onsegments = maptoseg(exons, segment_qranges)
				
				
				for cigar, exons_onsegment, size in zip(segment_cigars, exons_onsegments, segment_alignsizes):
					
					if len(exons_onsegment) == 0:
						continue
					
					name, cigar = cigar.split(":")
					
					if name[0] == '<':
						exons_onsegment = [size[1]-x for x in exons_onsegment]
						
					exons_onsegment_rposi = cigar_findrposi(cigar, exons_onsegment)
					
					segment_exonalign[gene][name[1:]].extend(exons_onsegment_rposi)
					
		for gene, segment_exonalign_ in segment_exonalign.items():	
			
			for name, exonalign in segment_exonalign_.items():
				
				exonalign = overlapexons(exonalign)
				
				exonalign = [x for x in exonalign if x[1]-x[0] > 50]
				
				segment_exonalign[gene][name] = exonalign
				
	return segment_exonalign


def plotexons(fig, ax, lelement, inputfile, gff3file, used_cigars):
	
	table = pd.read_csv(inputfile, sep = '\t', header = None)
	
	row0 = table.values.tolist()[0]
	segment_cigars = re.findall('[<>][0-9_A-Za-z:]+',row0[-2])
	
	segment_order = []
	for segment_cigar in segment_cigars:
			segment_order.append(segment_cigar.split(':')[0][1:])
		
	exonlocations = getgenelocation(table, gff3file, used_cigars)
	
	segment_order = [x for x in segment_order if x in used_cigars]
	exonlocations_order = {gene: [exonlocation[segment.split('_')[0]] for segment in segment_order]  for gene,exonlocation in exonlocations.items() }
	
	max_dupnum = max([int(x.split('_')[-1]) if '_' in x else 1 for x in segment_order])
	
	c = 'red'
	l = int(lelement * 0.005 + 0.5)
	L= l * 3
	chromy = lelement + l * 10
	fullsize = sum(list(used_cigars.values()))
	dupnum = 0
	for genename, exonlocation in exonlocations_order.items():
		
		if len(sum(exonlocation, [])) == 0:
			continue
		
		lastdupnum = 0
		leftend = 0
		exons_posi = []
		gene_posi = [[]]
		
		for name, pathexon in zip(segment_order, exonlocation):
			
			dupnum_ = int(name.split("_")[-1]) if "_" in name else 0
			
			size = used_cigars[name.split('_')[0]]
			rightend = leftend + size
			
			if dupnum_ > lastdupnum :
				gene_posi.append([])
				
				
			for exon in pathexon:
				gene_posi[-1].append([ (leftend + exon[0]) / fullsize, (leftend + exon[1]) / fullsize])
				exons_posi.append([ (leftend + exon[0]) / fullsize, (leftend + exon[1]) / fullsize])
			leftend = rightend
			
			lastdupnum = dupnum_
			
			
		gene_posi_ = [[min(sum(x,[])),max(sum(x,[]))] for x in gene_posi if len(x)]
		
		exons_posi_ = sum(exons_posi, [])
		
		if len(exons_posi_) == 0:
			continue
		
		for gene in gene_posi_ :
			
			ax.plot(gene, [chromy, chromy], color=c, linewidth=l, linestyle='-')
			
		if len(gene_posi_):
			ax.text(gene_posi_[0][0] - 0.05, chromy, genename, fontsize=l*8, color='red', ha='center', va='center')
			
			
		for exon in exons_posi:
			
			ax.plot(exon, [chromy, chromy], color=c, linewidth=L, linestyle='-')
			
		chromy += l * 10
		
	return 

def segfiltration(cigar, usednames = ""):
	
	cigars = re.findall('[<>][0-9_A-Za-z:]+',cigar)
	
	newcigars = []
	for cigarstr in cigars:
		
			thename,thecigar = cigarstr.split(":")
			if type(usednames) != type("") and thename[1:] not in usednames:
				
					continue
		
			thesize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=]', thecigar)]+[0])
		
			if thesize >0 and thesize < 50:
				
					fullsize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XHD]', thecigar)]+[0])
				
					thecigar = "{}H".format(fullsize )
				
			newcigars.append(thename+":"+thecigar)
		
	newcigars= "".join(newcigars)
	
	return newcigars

def getsegments(cigarstr,highlight,exons,ifmask = 0):
	
	cigars = re.findall(r'[><0-9_:]*\d+[=A-Za-z]+', cigarstr)
	
	
	variants = []
	segments = []
	highlight_points = {}
	exon_points = {}
	
	highlight_coordinate = -1
	highlight_type = cl.defaultdict(int)
	
	oldname = ""
	variant = []
	
	localpath_offsite = 0
	coordinate = 0
	laststart = 0
	findposition = -1
	
	for cigar in  cigars:
		
		if ":" in cigar:
			
			localpath_offsite = coordinate
			
			name,cigar = cigar.split(":")
			
			name = name[1:]
			
			if name in highlight:
				
				highlight_points[name] = [[x + coordinate for x in seg] for seg in highlight[name] ]
				
			if name in exons:
				
				exon_points[name] =  [[x + coordinate for x in seg] for seg in exons[name] ]
				
		size_str = re.findall(r'\d+',cigar)[0]
		size_str_len = len(size_str)
		
		size = int(size_str)
		thetype = cigar[size_str_len]
		
		ifmasked = cigar[size_str_len:].islower() and ifmask
		
		if thetype in ['D','H']:
			
			if size > 50:
				
				if coordinate - laststart > 500:
					
					segments.append([laststart , coordinate])
					
					if not ifmasked:
						variants.append(variant)
						
				laststart = coordinate + size
				variant = []
				
			else:
				variant.extend([x for x in range(coordinate,coordinate + size,20)])
				if name in highlight_points:
					for i,seg in enumerate(highlight_points[name]):
						
						if seg[0] <= coordinate + size and seg[1] >= coordinate :
							#highlight_type[(name,seg)] = 1
							
							highlight_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
							
							
				if name in exon_points:
					for i,seg in enumerate(exon_points[name]):
						if seg[0] <= coordinate + size and seg[1] >= coordinate :
							#highlight_type[(name,seg)] = 1
							
							exon_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
							
							
			coordinate += size
			
			
		elif thetype in ['X']:
			
			
			
			if name in highlight_points:
				for i,seg in enumerate(highlight_points[name]):
					
					if seg[0] <= coordinate + size and seg[1] >= coordinate :
						#highlight_type[(name,seg)] = 1
						
						highlight_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
						
						
			if name in exon_points:
				for i,seg in enumerate(exon_points[name]):
					if seg[0] <= coordinate + size and seg[1] >= coordinate :
						#highlight_type[(name,seg)] = 1
						
						exon_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
						
			if not ifmasked:
				variant.extend([x for x in range(coordinate,coordinate + size,20)])
				
				
			coordinate += size
			
		elif thetype in ['I']:
			
			if not ifmasked:
				variant.append(coordinate)
				
				
			size = 1
			if name in highlight_points:
				for i,seg in enumerate(highlight_points[name]):
					
					if seg[0] <= coordinate + size and seg[1] >= coordinate :
						#highlight_type[(name,seg)] = 1
						
						highlight_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
						
						
			if name in exon_points:
				for i,seg in enumerate(exon_points[name]):
					if seg[0] <= coordinate + size and seg[1] >= coordinate :
						#highlight_type[(name,seg)] = 1
						
						exon_points[name][i].extend([-x for x in range(max(coordinate, seg[0]),min(coordinate + size,seg[1])+1,20)])
						
						
		else:
			coordinate += size
			
	segments.append([laststart, coordinate])
	variants.append(variant)
	
	
	return segments, variants,highlight_points,exon_points


def plotmutant(alltypes,elements,highlightnames,fullsize):
	
	light_colors = [ "light"+color for color in ['yellow','skyblue' ,'gray','coral', 'seagreen' , 'steelblue', 'cyan', 'pink', 'green', 'grey','skyblue','salmon','slategray']]
	
	color_dict = {i: light_colors[i%len(light_colors)] for i in range(max(alltypes)+1)}
	
	plotsize = len(elements)/10
	matplotlib.rcParams["path.simplify_threshold"] = 0.00001
	matplotlib.rcParams['path.simplify'] = False
	
	fig,ax = plt.subplots(figsize=(plotsize, plotsize))
	
	c = 'red'
	l = 5
	L= 20
	
	
	allvariants_x = []
	allvariants_y = []
	
	allhighlights_x = []
	allhighlights_y = []
	
	allexons_x = []
	allexons_y = []
	#elements = elements[400:500]
	
	height = 4*plotsize/len(elements)
	width = 2*plotsize/fullsize
	
	
	theminsize = min( width,height) 
	repeattime = max( width,height) / theminsize
	repeattime = int(max(1, repeattime/10))
	
	offsite = 10
	lasttype = -1
	gapsize = 10
	repeattime = 1
	for i, (thename,segments, variants,hilis,exonpts) in enumerate(elements):
		
		if i % 100 == 0:
			print(i)
			
		thetype = alltypes[i]
		
		c = color_dict[thetype]
		if c != lasttype:
			offsite += gapsize
			
		lasttype = c
		
		i += offsite
		
		ifhighlight = 0
		L = 1
		if len([x for x in  highlightnames if x in thename]):
			c = adjust_lightness(c, amount=0.7)
			L = 3
			ifhighlight = 1
			
		breaks = [x/fullsize for y in segments for x in y[:2]]
		
		
		posix = [v/fullsize for x in variants for v in x] * repeattime
		posiy = sum([ [i+10*r*theminsize]*len(posix)  for r in range(repeattime) ], [])
		
		hilix = [-v/fullsize for v in hilis if v < 0 ] * repeattime
		hiliy = sum([ [i+10*r*theminsize]*len(hilix)  for r in range(repeattime) ], [])
		
		exonx = [-v/fullsize for v in exonpts if v < 0 ] * repeattime
		exony = sum([ [i+10*r*theminsize]*len(exonpts)  for r in range(repeattime) ], [])
		
		
		allvariants_x.extend(posix)
		allvariants_y.extend( posiy)
		
		allhighlights_x.extend(hilix)
		allhighlights_y.extend(hiliy)
		
		allexons_x.extend(exonx)
		allexons_y.extend(exony)
		
		last_coordinate = 0
		for index, coordinate in enumerate(breaks):
			
			if index % 2 == 0:
				l = 0.2
				ax.plot([last_coordinate, coordinate], [i, i], color=c, linewidth=l, linestyle='-.')
			else:
				l = 2
				ax.plot([last_coordinate, coordinate], [i, i], color=c, linewidth=l, linestyle='-')
				
			if ifhighlight:
				ax.scatter(-0.05, i, marker='o', color='red', s=L*50, zorder=5)
				
				
			last_coordinate = coordinate
			
			
	ax.scatter(allvariants_x, allvariants_y, marker='.',color =  'black'  ,s = 5,zorder=3,linestyle='None')
	ax.scatter(allexons_x, allexons_y, marker='.',color =  'blue'  ,edgecolor='blue',s = 10,zorder=10,linestyle='None',alpha=1, cmap='viridis')
	ax.scatter(allhighlights_x, allhighlights_y, marker='.',color =  'red'  ,s = 10,zorder=3,linestyle='None')
	
	return fig,ax,i

def loadmsafile(inputfile, typefile, genename = ""):
	
	names = []
	cigars = []
	types = dict()
	
	with open(inputfile, mode = 'r') as f:
		
		for line in f:
			
			if len(line.strip()) == 0:
				continue
			
			if len(genename) and line.startswith(genename) == False:
				continue
			
			if line[0] == 'L':
				
				names.append(line[1]) 
				cigars.append(line[5])
				
	types = cl.defaultdict(int)
	with open(typefile, mode = 'r') as f:
		for line in f:
			if len(line.strip()) == 0:
				continue
			
			line = line.strip().split('\t')
			for name in line[-1].split(","):
				types[name] = len(typeinfo)
				
	return names, cigars, types

def loadannofile(inputfile, genename=""):
	names = []
	cigars = []
	types = {}
	
	if inputfile.endswith('.gz'):
		open_func = lambda f: gzip.open(f, mode='rt')
	else:
		open_func = open
		
	with open(inputfile) as f:
		for line in f:
			
			if len(line.strip()) == 0 or line.startswith("#"):
				continue
			
			if genename and not line.startswith(genename):
				continue
			fields = line.rstrip('\n').split('\t')
			if len(fields) < 18:
				continue  # skip malformed lines
			name = fields[0]
			thetype = fields[1]
			cigar = fields[17]
			type_int = int(thetype.split("_")[2])
			names.append(name)
			cigars.append(cigar)
			types[name] = type_int
			
	return names, cigars, types



def selectcigar(names,cigars,types):
	
	used_cigars = dict()
	for i,(name,cigar) in enumerate(zip(names,cigars)):
		
		segment_cigars = cigar.replace('<','>').split(">")[1:]
		segment_cigars = [x for x in segment_cigars if not re.match(r"[0-9H]+$", x.split(":")[1]) ]
		
		for segment_cigar in segment_cigars:
			
			thesize = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XHD]', segment_cigar)]+[0])
			
			if thesize>= 300:
				used_cigars[segment_cigar.split(":")[0].split('_')[0]] = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XHD]', segment_cigar)]+[0])
				used_cigars[segment_cigar.split(":")[0]] = sum([int(x[:-1]) for x in re.findall(r'\d+[M=XHD]', segment_cigar)]+[0])
				
		if name in types:
				pass
		elif "_".join(name.split("_")[2:]) in types:
				pass
		else:   
				continue
		
	return used_cigars

def segment_order(cigars,used_cigars):
	
	cigars0 = cigars[0]
	segment_cigars = re.findall('[<>][0-9_A-Za-z:]+',cigars0)
	
	segment_orders = []
	for segment_cigar in segment_cigars:
		segment_orders.append(segment_cigar.split(':')[0][1:])
		
	segment_orders = [x.split('_')[0] for x in segment_orders if x in used_cigars]
	
	return segment_orders


def getvariants(names,cigars,types,used_cigars,highlight,exons,ifmask):
	
	elements = []
	counter = 0 
	alltypes = []
	for i,(name,cigar) in enumerate(zip(names,cigars)):
		
			cigar = segfiltration(cigar,used_cigars)
			if name in types:
					alltypes.append(types[name])
			elif "_".join(name.split("_")[2:]) in types:
					alltypes.append(types["_".join(name.split("_")[2:])])
			else:
					print("skipping: ",name)
					continue
					#alltypes.append(0)
		
			segments, variants, highlight_points,exon_points = getsegments(cigar,highlight,exons,ifmask)
		
			if len(highlight_points):
					counter += 1
				
			highlight_points = [y for seg in list(highlight_points.values()) for x in seg for y in x[2:]]
		
			exon_points =  [y for seg in list(exon_points.values()) for x in seg for y in x[2:]]
		
			elements.append([name, segments, variants,highlight_points, exon_points] )
		
	return alltypes, elements

def extract_names(inputfile, genename):
	output = []
	with open(inputfile) as f:
		for line in f:
			line = line.rstrip('\n')
			if line.startswith("results:"):
				elements = line[len("results:"):].strip().split(',')
				for ele in elements:
					if ele.startswith(genename):
						output.append(ele)
			elif line.startswith(genename):
				first_elem = line.split('\t')[0]
				output.append(first_elem)
	return output


def makemutant(inputfile, outputfile, annotation, hnames=None, gff3file = "" ,ifintron = 1, ifmask = 0, genename = ""):
	
	exons = {}
	highlight = {}
	
	names, cigars, types = loadannofile(annotation, genename)
	
	used_cigars = selectcigar(names,cigars,types)
	
	segment_orders = segment_order(cigars,used_cigars)
	
	fullsize = sum(used_cigars.values())
	
	alltypes, elements = getvariants(names,cigars,types,used_cigars,highlight,exons,ifmask)
	
	if len(hnames):
		highlightnames = [x for x in hnames.split(",") if len(x)]
	elif len(inputfile):
		highlightnames = extract_names(inputfile, genename)
		
	fig,ax, offsite = plotmutant(alltypes,elements,highlightnames,fullsize)
	
	# Show the plot
	if len(gff3file):
		
		fullfile = inputfile
		
		plotexons(fig, ax, offsite , fullfile , gff3file, used_cigars)
		
	plt.axis('off')
	plt.savefig(outputfile)
	
	
def main(args):
	
	makemutant(args.input, args.output, args.anno, args.name,  args.gff, args.intron, args.mask, args.gene)
	
def run():
	"""
		Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="program to visualize pangenome-alleles")
	parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, default = "")
	parser.add_argument("-a", "--anno", help="path to annotation dababase file",dest="anno", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
	parser.add_argument("-g", "--gene", help="name of the plotting gene or gene group", dest="gene",type=str, required=True)
	parser.add_argument("-n", "--name", help="highlight names, separate by comma, or input the ctyper genotyping results/annotation table file", dest="name",type=str, default = "")
	parser.add_argument("-intron", "--intron", help="if plot intron duplications", dest="intron",type=int, default = 0)
	parser.add_argument("-mask", "--mask", help="if hidden variants in masked region", dest="mask",type=int, default = 0)
	parser.add_argument("-G", "--gff", help="path to gff3 file to plot GRCH38 genes", dest="gff",type=str, default = "")
	
	
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
