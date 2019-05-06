#! /usr/bin/env python

from __future__ import print_function
import cyvcf2
import argparse
import collections
import numpy
import sys
import csv
import gzip


# Parse command line
parser = argparse.ArgumentParser(description='Convert BCF file to link file')
parser.add_argument('-i', '--input', metavar='delly.bcf', required=True, dest='svinput', help='input SV file (required)')
parser.add_argument('-f', '--format', metavar='delly', required=False, default='delly', dest='svformat', help='SV format [delly|bedpe]')
parser.add_argument('-m', '--min', metavar='20000000', type=int, required=False, default=20000000, dest='minSize', help='min. SV size')
args = parser.parse_args()

# Arguments
minSize = args.minSize

# Parse Delly BCF file
if args.svformat == "delly":
    vcf = cyvcf2.VCF(args.svinput)
    samples = numpy.array(vcf.samples)
    for record in vcf:
        if (record.CHROM == "Y") or (record.CHROM == "chrY"):
            continue
        if len(record.ALT)==1:
            if record.is_snp:
                continue
            elif record.is_indel:
                continue
            else:
                pass
        else:
            continue
        chr1 = record.CHROM.replace("chr","")
        pos1 = record.POS
        chr2 = record.INFO["CHR2"].replace("chr","")
        pos2 = record.INFO["END"]
        svt = record.INFO["SVTYPE"]
        ct = record.INFO["CT"]
        if (chr1 != chr2) or (pos1 + minSize < pos2):
            print(chr1, pos1, pos1 + 1, chr2, pos2, pos2 + 1, "svtype=" + svt + ",ct=" + ct, sep=" ")
elif args.svformat == "bedpe":
    with gzip.open(args.svinput) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            chr1 = row['chrom1']
            pos1 = int(row['start1'])
            chr2 = row['chrom2']
            pos2 = int(row['end1'])
            svt = row['svclass']
            if svt == "DEL":
                ct = "3to5"
            if svt == "DUP":
                ct = "5to3"
            if svt == "h2hINV":
                svt = "INV"
                ct = "5to5"
            if svt == "t2tINV":
                svt = "INV"
                ct = "3to3"
            if svt == "TRA":
                svt = "BND"
                ct = "NtoN"
            if (chr1 != chr2) or (pos1 + minSize < pos2):
                print(chr1, pos1, pos1 + 1, chr2, pos2, pos2 + 1, "svtype=" + svt + ",ct=" + ct, sep=" ")

            

        
