#!/usr/bin/env python

import csv, os, sys
from os.path import basename, splitext

infile = sys.argv[1]
infbase = splitext(basename(infile))[0]
outfile = 'convCLCtable_' + infbase + '_.tsv'

inf = open(infile,'r')
outf = open(outfile, 'w')

tsvinreader = csv.reader(inf, delimiter='\t')
tsvoutwriter = csv.writer(outf, delimiter='\t')
for row in tsvinreader:
    position = (int(row[1])) -1
    tsvoutwriter.writerow([row[0], str(position), row[1]])

