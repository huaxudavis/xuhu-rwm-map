#!/usr/bin/env python

# !!! requires python 2.6 or 2.7 !!!

# This script creates a "BED file" with SNP locations from a table of SNPs exported from CLC as a tab-delimited text file.
# Please export only "SNV"s from CLC (and "MNV"s if newer version of CLC is used i.e. starting from CLC 6.05).
# CLC exports the variants as "SNV" and since version CLC 6.05 as both "SNV" and "MNV" (the latter variants are up to 5 bases long).
#
# By default the script will use only SNVs and the first base of MNVs.
# If you want to use the bases of MNVs as separate SNPs and split them into SNPs please provide the max. length as
# an integer number (smaller than 6) as an aditional argument for the script. Please see below.
# For SNP location tables exported from CLC 6.04 or earlier all records (SNVs as well as the first bases of MNVs) will be retained as SNPs. 
#
# input arguments:
#   1.input file.
#   2. optional: Max length of MNV to use (integer  < 6)  - for example use "2" if you do not want to retain MNVs at all.
#                The default is: 2
#   3. optional: Number of MNV bases to convert into SNPs  (integer < 6)
#
# Usage:         python 01-CLC_Table_Converter.py  tab-separated-CLC-table.txt &
# Optional:      python 01-CLC_Table_Converter.py  tab-separated-CLC-table.txt   Number-of-MNV-bases-to-use   Max-Length-of-MNV  &
# E.g.           python 01-CLC_Table_Converter.py  tab-separated-CLC-table.txt   2  3  &
# The resulting tab-delimited BED file has threee columns: chromosome-name, coordinate-before-SNP, coordinate-of-the-SNP

import csv, sys
from os.path import basename, splitext

infile = sys.argv[1]
mnvb = 1
max = 1

if len(sys.argv) == 3 or len(sys.argv) == 4:
    max =  int(sys.argv[2]) + 1

if len(sys.argv) == 4:
    mnvb = int(sys.argv[3])
    if mnvb > 5:
        mnvb = 5


infbase = splitext(basename(infile))[0]
outfile = 'SNP-table_' + infbase + '.bed.tsv'

inf = open(infile,'rb')
outf = open(outfile, 'wb')

tsvinreader = csv.reader(inf, delimiter='\t')
tsvoutwriter = csv.writer(outf, delimiter='\t')

next(tsvinreader)

for row in tsvinreader:
# print int(row[1])
    length = 1
    startpos = 0
    endpos= 0
    position = 0

    if ".." in row[1]:
        startpos = int(row[1].split("..")[0])
        if startpos > 0:
            position = startpos -1
        endpos = int(row[1].split("..")[1])
        length = endpos - startpos + 1
    else:
        startpos = int(row[1])
        if startpos > 0:
            position = startpos -1
        length = 1
        
    if length > max:
        pass
    else:
        if position > 0:
            tsvoutwriter.writerow([row[0], str(position), str(position + 1)])
        
        if length > 1:
            if mnvb > 1:
                tsvoutwriter.writerow([row[0], str(position + 1), str(position + 2)])
            else:
                pass
        
        if length > 2:
            if mnvb > 2:
                tsvoutwriter.writerow([row[0], str(position + 2), str(position + 3)])
            else:
                pass
        
        if length > 3:
            if mnvb > 3:
                tsvoutwriter.writerow([row[0], str(position + 3), str(position + 4)])
            else:
                pass
        
        if length > 4:
            if mnvb > 4:
                tsvoutwriter.writerow([row[0], str(position + 4), str(position + 5)])
            else:
                pass



