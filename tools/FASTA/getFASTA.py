#!/usr/bin/env python3

import argparse
import gzip
import re


def read_fasta(fasta_path, seq_dict):
    read = []
    header = None
    with open(fasta_path, mode='r') as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith(">"):
                if header:
                    seq_dict[header] = "".join(read)
                header = line[1:].split()[0]
                read = []
            else:
                read.append(line.strip())
        if header:
            seq_dict[header] = "".join(read)
    return seq_dict


def make_reverse(seq, strand="-"):
    if strand == "+":
        return seq
    tran = str.maketrans("ATCGatcg", "TAGCtagc")
    return seq[::-1].translate(tran)


def cigar_break(cigar):
    parsed = []
    for match in re.finditer(r'(\d+)([A-Za-z]+)([^0-9A-Za-z]*)', cigar):
        length, op, extra = match.groups()
        parsed.append((op, int(length), extra))
    return parsed


def cigar_to_seq(rlocation, cigar, ref_seq):
    cigar_ops = cigar_break(cigar)
    strand = 1 if rlocation[-1] == '+' else -1
    start, end = map(int, rlocation[:-1].split('-'))

    rpos = end if strand == -1 else start - 1
    segments = []

    for op, length, seq in cigar_ops:
        if op in {"M", "="}:
            if strand > 0:
                segments.append(ref_seq[rpos:rpos + length])
            else:
                segments.append(make_reverse(ref_seq[rpos - length:rpos]))
        elif op == "I" and strand > 0:
            segments.append(seq[2:])  # assumes extra is >2 chars, else fix
        elif op == "X":
            segments.append(seq[length:])
        if op in {"M", "D", "X"}:
            rpos += strand * length

    return "".join(segments)


def get_fasta(input_path, anno_path, ref_files, output_path):
    target_names = set()
    with open(input_path, mode='r') as f:
        for line in f:
            if line.strip():
                target_names.add(line.strip().split("\t")[0])

    ref_seqs = {}
    for ref_file in ref_files.split(","):
        read_fasta(ref_file, ref_seqs)

    if anno_path.endswith(".gz"):
        anno_file = gzip.open(anno_path, mode='rt')
    else:
        anno_file = open(anno_path, mode='r')

    with open(output_path, mode='w') as out:
        for line in anno_file:
            line = line.strip().split("\t")
            if not line:
                continue
            name = line[0]
            if name not in target_names:
                continue

            qlocation, alignment = line[7], line[-1]
            try:
                chrom, rlocation, cigar = alignment.split(":")[-3:]
            except ValueError:
                out.write(f">{name} {qlocation} {alignment}\nINVALID_ALIGNMENT\n")
                continue

            if rlocation == "NA" or chrom == "NA":
                out.write(f">{name} {qlocation} {alignment}\nNA\n")
                continue

            if chrom not in ref_seqs:
                print(f"[WARN] Missing reference contig '{chrom}' — check for alternative contigs or use -r HG38_main.fa,$Ctyper_PATH/data/Allalters.fa")
                continue

            seq = cigar_to_seq(rlocation, cigar, ref_seqs[chrom])
            out.write(f">{name} {qlocation} {chrom}:{rlocation}\n{seq}\n")

    anno_file.close()


def main(args):
    get_fasta(args.input, args.anno, args.ref, args.output)


def run():
    parser = argparse.ArgumentParser(description="Convert Ctyper genotyping results to FASTA format.")
    parser.add_argument("-i", "--input", required=True, type=str, help="Input tab-delimited file with target names.")
    parser.add_argument("-r", "--ref", required=True, type=str, help="Comma-separated reference FASTA files.")
    parser.add_argument("-a", "--anno", required=True, type=str, help="Annotation file (can be .gz).")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output FASTA file.")

    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    run()
