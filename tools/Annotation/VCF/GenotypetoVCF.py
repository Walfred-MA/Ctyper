#!/usr/bin/env python3

import argparse
import collections as cl
import re 

def makereverse(seq):
    
    tran=str.maketrans('ATCGatcg', 'TAGCtagc')
    
    return seq[::-1].translate(tran)

def CIGARbreak(CIGAR):
    
    cigar_list = []
    allcigars = re.findall(r'\d+[A-Za-z]+', CIGAR)
    for cigar in allcigars:
        
        size_str = re.findall(r'^\d+',cigar)[0]
        seq = cigar[(len(size_str) +1):]
        thetype = cigar[(len(size_str))]
        cigar_list.append((thetype,int(size_str), seq))
        
    return cigar_list

def CIGAR_polish(rstart, rend,strd, qstart, qend,qstrd, allcigars):
    
    leftsize = 0
    leftindex = 0
    
    rightsize = 0
    rightindex = 0
    
    for i,cigar in enumerate(allcigars):
        
        thetype,length, seq = cigar
        
        if thetype == 'M' and length > 30:
            leftindex = i
            break
        
    for i,cigar in enumerate(allcigars[::-1]):
        
        thetype,length, seq = cigar
        
        if thetype == 'M' and length > 30:
            rightindex = len(allcigars) - i
            break
        
    leftsize = sum([x[1] for x in allcigars[:leftindex] if x[1] in ['M','D','X']]+[0])
    rightsize = sum([x[1] for x in allcigars[rightindex:] if x[1] in ['M','D','X']]+[0])
    
    if strd == '+':
        rstart, rend = rstart+leftsize, rend-rightsize
    else:
        rstart, rend = rstart-leftsize, rend+rightsize
        
    leftsize = sum([x[1] for x in allcigars[:leftindex] if x[1] in ['M','I','X']]+[0])
    rightsize = sum([x[1] for x in allcigars[rightindex:] if x[1] in ['M','I','X']]+[0])
    
    if qstrd == '+':
        qstart, qend = qstart+leftsize, qend-rightsize
    else:
        qstart, qend = qstart-leftsize, qend+rightsize
        
        
        
    return rstart, rend, qstart, qend, allcigars[leftindex:rightindex]

def variant_totext(chr, pos, data):
    
    type, size, seq, name, qcontig, qpos,qstrd, rstrd = data
    
    if type == "X":
        REF = seq[:len(seq)-size]
        ALT = seq[-size:]
    if type == "D":
        REF = seq
        ALT = "<DEL>"
    if type == 'I':
        REF = seq[0]
        ALT = seq[-size:]
        
    if rstrd == -1:
        if type in ['X'] and size > 1:
            REF = makereverse(REF)
            ALT = makereverse(ALT)
        elif type in ['D']:
            
            REF = makereverse(REF)
        elif type in ['I']:
            REF = makereverse(seq[len(seq) - size - 1])
            ALT = makereverse(ALT)
            
    info = {}
    if type in ['D','I']:
        
        info["END"] = str(pos+size-1)
        info["SVTYPE"] = "DEL" if type == 'D' else "INS"
        
        
    info["NAME"] = name  
    info["SOURCE"] = "{}:{}{}".format(qcontig, qpos,"+" if qstrd==1 else '-')
    columns = [chr, str(pos), ".", REF, ALT, ".", ".", info, "GT"]
    
    return columns


class vcfdata:
    
    def __init__(self):
        
        self.colnames = "#CHROM POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE"
        
        self.allchrom_variants = cl.defaultdict(lambda: cl.defaultdict(list))
        
        self.coverages =  cl.defaultdict(list)
        
    def determine_cover(self, chr, pos):
        
        
        return len([x for x in self.coverages[chr] if x[1] <= pos and x[2] > pos])
    
    
    def loadCIGAR(self, name, chr, rstart, rend,rstrd, qcontig, qstart, qend,qstrd, CIGAR):
        
        allcigars = CIGARbreak(CIGAR)
        
        rstart, rend, qstart, qend, allcigars = CIGAR_polish(rstart, rend,rstrd, qstart, qend,qstrd, allcigars)
        
        variants = self.allchrom_variants[chr]
        
        qstrd = 1 if qstrd == '+' else -1
        if qstrd == -1:
            qstart_ = qend - 1
            qend = qstart-1
            qstart = qstart_
            
        rstrd = 1 if rstrd == '+' else -1
        if rstrd == -1:
            rstart_ = rend - 1
            rend = rstart-1
            rstart = rstart_
            
        rpos = rstart
        qpos = qstart
        for cigar in allcigars:
            
            thetype,length, seq = cigar
            
            
            if thetype == "I" and rstrd > 0:
                variants[rpos-1].append(( thetype,length, seq, name, qcontig, qpos, qstrd, rstrd))
            elif thetype in ['D','X'] and rstrd < 0:
                variants[rpos - length +1].append(( thetype,length, seq, name, qcontig, qpos, qstrd, rstrd))
            elif thetype != "M" and thetype != "=":
                variants[rpos].append(( thetype,length, seq, name, qcontig, qpos, qstrd, rstrd))
                
                
            if thetype in ['M','D','X']:
                rpos += rstrd * length
                
            if thetype in ['M','I','X']:
                qpos += qstrd * length
                
                
                
    def output(self, outputfile):
        
        self.allchrom_variants = {chr:sorted([(rpos, vars) for rpos,vars in variants.items()]) for chr, variants in self.allchrom_variants.items()}
        
        
        with open(outputfile, mode = 'w') as f:
            
            f.write(self.colnames+"\n")
            
            for chr, variants in self.allchrom_variants.items():
                
                
                for variant in variants:
                    
                    text = ['','','','','']
                    names = set()
                    totalcounter = len(set([x[3] for x in variant[1]]))
                    for var in variant[1]:
                        newtext = variant_totext(chr, variant[0], var)
                        
                        if text[:5] == newtext[:5] :
                            text[-2]["NAME"]+= ","+newtext[-2]["NAME"]
                            text[-2]["SOURCE"]+= ","+newtext[-2]["SOURCE"]
                            names.add(var[3])
                        else:
                            
                            if  len(text[0]):
                                samplecover = self.determine_cover(chr, variant[0])
                                
                                if max(2,samplecover)>len(names):
                                    sample = "0/1"
                                else:
                                    sample = "1/1"
                                    
                                    
                                text[-2] = ";".join([k+"="+x for k,x in text[-2].items()])
                                f.write("\t".join(text+[sample])+"\n")
                                
                            text = newtext
                            names = set([var[3]])
                            
                            
                    if  len(text[0]):
                        samplecover = self.determine_cover(chr, variant[0])
                        if max(2,samplecover)>len(names):
                            sample = "0/1"
                        else:
                            sample = "1/1"
                            
                        text[-2] = ";".join([k+"="+x for k,x in text[-2].items()])
                        f.write("\t".join(text+[sample])+"\n")
                        
                        
                        
        coverages = [ sorted(coverage, key = lambda x: (len(x[0]), x[0], int(x[1]))) for chr, coverage in self.coverages.items()]
        
        with open(outputfile+".bed",mode = 'w') as f:
            
            for coverage in coverages:
                
                for cover in coverage:
                    
                    f.write("\t".join([str(x) for x in cover]) + "\n")
                    
                    
def vcf(inputfile, annofile, outputfile):
    
    allnames = set()
    with open(inputfile, mode = 'r') as f:
        
        for line in f:
            
            line = line.split("\t")
            allnames.add(line[0])
            
    fullcigars = cl.defaultdict(str)
    with open(annofile, mode = 'r') as f:
        
        for line in f:
            
            line = line.split("\t")
            
            if line[0] in allnames:
               
                if line[-1].count(":") > 1: 
                    data = line[-1].split(":")
                    if min([len(x.strip()) for x in data]) > 0:
                        fullcigars[line[0]] = (data[-3],data[-2],data[-1])
                
                
    vcf_record = vcfdata()
    
    name_counter = cl.defaultdict(int)
    results = set()
    with open(inputfile, mode = 'r') as f:
        
        for line in f:
            
            line = line.split("\t")
            name = line[0]
            name_counter[name] += 1
            
            if name_counter[name] > 1:
                name += "#" + str(name_counter[name])
                
            alignment = fullcigars[name]
            if len(alignment) == 0:
                continue
            
            qlocation = line[5]
            chr, location, CIGAR = alignment
            strd = location[-1]
            rstart, rend = location[:-1].split('-')
            rstart, rend = int(rstart), int(rend)
            
            
            vcf_record.coverages[chr].append([chr, rstart, rend,name])
            
            qcontig, qlocation = qlocation.split(":")
            qstrd = qlocation[-1]
            qstart, qend = qlocation[:-1].split('-')
            qstart, qend = int(qstart), int(qend)
            
            vcf_record.loadCIGAR(name, chr, rstart, rend, strd, qcontig, qstart, qend, qstrd, CIGAR)
            
    vcf_record.output(outputfile)
    
    
    
    
def main(args):
    
    vcf(args.input, args.anno, args.output)
    
def run():
    """
            Parse arguments and run
    """
    parser = argparse.ArgumentParser(description="program to convert ctyper genotyping results to VCF")
    parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
    parser.add_argument("-a", "--anno", help="path to annotation data file",dest="anno", type=str, required=True)
    parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
    
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
    
    
if __name__ == "__main__":
    run()
