#! /usr/bin/env python

from __future__ import print_function
from readfq import readfq
import cyvcf2
import argparse
import collections
import numpy
import sys
import csv
import gzip


def rev(c):
    if c == 'A':
        return 'T'
    elif c == 'C':
        return 'G'
    elif c == 'G':
        return 'C'
    elif c == 'T':
        return 'A'
    else:
        return 'N'

# Parse command line
parser = argparse.ArgumentParser(description='SNV statistics for BCF')
parser.add_argument('-b', '--bcf', metavar='snv.bcf', required=True, dest='bcf', help='BCF file (required)')
parser.add_argument('-r', '--ref', metavar='ref.fa', required=True, dest='ref', help='input reference (required)')
args = parser.parse_args()

# Get all the positions
bcf = cyvcf2.VCF(args.bcf)
selectedPos=collections.defaultdict(set)
for record in bcf:
    selectedPos[record.CHROM].add(record.POS)

# Store the prefix and suffix reference nucleotide
nuclPrevDict = collections.defaultdict(set)
nuclPostDict = collections.defaultdict(set)
f_in = gzip.open(args.ref) if args.ref.endswith('.gz') else open(args.ref)
for seqName, seqNuc, seqQuals in readfq(f_in):
    for pos in selectedPos[seqName]:
        nuclPrevDict[ (seqName, pos) ] = seqNuc[(pos-2):(pos-1)]
        nuclPostDict[ (seqName, pos) ] = seqNuc[pos:(pos+1)]


# Build mutation dictionary    
mt = dict()
for i in ['A', 'C', 'G', 'T']:
    for j in ['A', 'C', 'G', 'T']:
        for k in ['A', 'C', 'G', 'T']:
            for l in ['A', 'C', 'G', 'T']:
                if j != k:
                    if (j == 'C') or (j == 'T'):
                        mt[(i, j, k, l)] = (i, j, k, l)
                    else:
                        mt[(i, j, k, l)] = (rev(l), rev(j), rev(k), rev(i))
    
# VCF
vcf = cyvcf2.VCF(args.bcf)
for record in vcf:
    if (record.CHROM == "Y") or (record.CHROM == "X"):
        continue
    if len(record.ALT)==1:
        if record.is_snp:
            pass
        elif record.is_indel:
            continue
        else:
            continue
    else:
        continue
    if record.FILTER is None:
        m = mt[(nuclPrevDict[(record.CHROM, record.POS)], record.REF, record.ALT[0], nuclPostDict[(record.CHROM, record.POS)])]
        mtype = str(m[1]) + ">" + str(m[2])
        col = "puor-6-div-1"
        if mtype == "C>A":
            col = "puor-6-div-1"
        elif mtype == "C>G":
            col = "puor-6-div-2"
        elif mtype == "C>T":
            col = "puor-6-div-3"
        elif mtype == "T>A":
            col = "puor-6-div-4"
        elif mtype == "T>C":
            col = "puor-6-div-5"
        elif mtype == "T>G":
            col = "puor-6-div-6"
        offset = 0
        prec = str(m[0])
        if prec == 'A':
            offset = 0
        elif prec == 'C':
            offset = 1
        elif prec == 'G':
            offset = 2
        elif prec == 'T':
            offset = 3
        offset *= 4
        succ = str(m[3])
        if succ == 'A':
            offset += 0
        elif succ == 'C':
            offset += 1
        elif succ == 'G':
            offset += 2
        elif succ == 'T':
            offset += 3
        print(record.CHROM, record.POS, record.POS, offset, "color=" + col + ",mut=" + prec + mtype + succ, sep=" ")
