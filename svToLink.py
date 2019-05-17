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
parser.add_argument('-o', '--outprefix', metavar='out', required=False, default='out', dest='outprefix', help='output prefix')
parser.add_argument('-p', '--pos', metavar='59283206', type=int, required=False, default=-1, dest='posHighlight', help='highlight position')
parser.add_argument('-m', '--min', metavar='10000000', type=int, required=False, default=10000000, dest='minSize', help='min. SV size')
args = parser.parse_args()

# Arguments
minSize = args.minSize

fout = args.outprefix + ".links.txt"
with open(fout, 'w') as flinks:
    foutDel = args.outprefix + ".tiles.del.txt"
    with open(foutDel, 'w') as fdel:
        foutDup = args.outprefix + ".tiles.dup.txt"
        with open(foutDup, 'w') as fdup:
            foutInv3 = args.outprefix + ".tiles.inv3to3.txt"
            with open(foutInv3, 'w') as finv3:
                foutInv5 = args.outprefix + ".tiles.inv5to5.txt"
                with open(foutInv5, 'w') as finv5:
                    # Parse Delly BCF file
                    if args.svformat == "delly":
                        vcf = cyvcf2.VCF(args.svinput)
                        samples = numpy.array(vcf.samples)
                        for record in vcf:
                            highlight=0
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
                            if pos1 == args.posHighlight:
                                highlight = 1
                            if (chr1 == "Y") or (chr1 == "chrY") or (chr2 == "Y") or (chr2 == "chrY"):
                                continue
                            svt = record.INFO["SVTYPE"]
                            ct = record.INFO["CT"]
                            if (chr1 != chr2) or (pos1 + minSize < pos2):
                                print(chr1, pos1, pos1 + 1, chr2, pos2, pos2 + 1, "svtype=" + svt + ",ct=" + ct + ",hl=" + str(highlight), sep=" ", file=flinks)
                            else:
                                if svt == 'DEL':
                                    print(chr1, pos1, pos2, "svtype=" + svt + ",ct=" + ct + ",hl=" + str(highlight) + ",color=dark2-5-qual-2", sep=" ", file=fdel)
                                elif svt == 'DUP':
                                    print(chr1, pos1, pos2, "svtype=" + svt + ",ct=" + ct + ",hl=" + str(highlight) + ",color=dark2-5-qual-3", sep=" ", file=fdup)
                                elif (svt == 'INV') and (ct == '3to3'):
                                    print(chr1, pos1, pos2, "svtype=" + svt + ",ct=" + ct + ",hl=" + str(highlight) + ",color=dark2-5-qual-4", sep=" ", file=finv3)
                                elif (svt == 'INV') and (ct == '5to5'):
                                    print(chr1, pos1, pos2, "svtype=" + svt + ",ct=" + ct + ",hl=" + str(highlight) + ",color=dark2-5-qual-5", sep=" ", file=finv5)
                    elif args.svformat == "bedpe":
                        with gzip.open(args.svinput) as f:
                            reader = csv.DictReader(f, delimiter="\t")
                            for row in reader:
                                highlight = 0
                                chr1 = row['chrom1']
                                pos1 = (int(row['start1']) + int(row['end1'])) / 2
                                chr2 = row['chrom2']
                                pos2 = (int(row['start2']) + int(row['end2'])) / 2
                                if pos1 == args.posHighlight:
                                    highlight = 1
                                if (chr1 == "Y") or (chr1 == "chrY") or (chr2 == "Y") or (chr2 == "chrY"):
                                    continue
                                svt = row['svclass']                                
                                if svt == "DEL":
                                    ct = "3to5"
                                    col = "dark2-5-qual-2"
                                elif svt == "DUP":
                                    ct = "5to3"
                                    col = "dark2-5-qual-3"
                                elif svt == "h2hINV":
                                    svt = "INV"
                                    ct = "5to5"
                                    col = "dark2-5-qual-5"
                                elif svt == "t2tINV":
                                    svt = "INV"
                                    ct = "3to3"
                                    col = "dark2-5-qual-4"
                                elif svt == "TRA":
                                    svt = "BND"
                                    ct = "NtoN"
                                    col = "dark2-5-qual-1"
                                if (chr1 != chr2) or (pos1 + minSize < pos2):
                                    print(chr1, pos1, pos1 + 1, chr2, pos2, pos2 + 1, "svtype=" + svt + ",ct=" + ct + ",hl=" + str(highlight) + ",color=" + col, sep=" ", file=flinks)
                                else:
                                    if svt == 'DEL':
                                        print(chr1, pos1, pos2, "svtype=" + svt + ",ct=" + ct + ",hl=" + str(highlight) + ",color=" + col, sep=" ", file=fdel)
                                    elif svt == 'DUP':
                                        print(chr1, pos1, pos2, "svtype=" + svt + ",ct=" + ct + ",hl=" + str(highlight) + ",color=" + col, sep=" ", file=fdup)
                                    elif (svt == 'INV') and (ct == '3to3'):
                                        print(chr1, pos1, pos2, "svtype=" + svt + ",ct=" + ct + ",hl=" + str(highlight) + ",color=" + col, sep=" ", file=finv3)
                                    elif (svt == 'INV') and (ct == '5to5'):
                                        print(chr1, pos1, pos2, "svtype=" + svt + ",ct=" + ct + ",hl=" + str(highlight) + ",color=" + col, sep=" ", file=finv5)
