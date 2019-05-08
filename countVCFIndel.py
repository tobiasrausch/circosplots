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
parser = argparse.ArgumentParser(description='InDel statistics for BCF')
parser.add_argument('-b', '--bcf', metavar='indel.bcf', required=True, dest='bcf', help='BCF file (required)')
parser.add_argument('-o', '--outprefix', metavar='out', required=False, default='out', dest='outprefix', help='output prefix')
args = parser.parse_args()

fdel = args.outprefix + ".small.del.txt"
with open(fdel, 'w') as fd:
    fins = args.outprefix + ".small.ins.txt"
    with open(fins, 'w') as fi:
        # Phased calls BCF1
        bcf = dict()
        vcf = cyvcf2.VCF(args.bcf)
        for record in vcf:
            if record.CHROM == "Y":
                continue
            if len(record.ALT)==1:
                if record.is_snp:
                    continue
                elif record.is_indel:
                    pass
                else:
                    continue
            else:
                continue
            if record.FILTER is not None:
                continue
            if len(record.REF) > len(record.ALT[0]):
                size = len(record.REF) - len(record.ALT[0])
                col = "rdylgn-5-div-1"
                if size <= 3:
                    off = 0
                elif size <= 9:
                    off = 1
                else:
                    off = 2
                print(record.CHROM, record.POS, record.POS, off, "color=" + col, sep=" ", file=fd)
            else:
                size = len(record.ALT[0]) - len(record.REF)
                col = "rdylgn-5-div-4"
                if size <= 3:
                    off = 0
                elif size <= 9:
                    off = 1
                else:
                    off = 2
                print(record.CHROM, record.POS, record.POS, off, "color=" + col, sep=" ", file=fi)
