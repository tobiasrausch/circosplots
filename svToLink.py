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
parser.add_argument('-b', '--bcf', metavar='indel.bcf', required=True, dest='bcf', help='BCF file (required)')
args = parser.parse_args()

minSize = 20000000

vcf = cyvcf2.VCF(args.bcf)
samples = numpy.array(vcf.samples)
count = 1
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
    count += 1
