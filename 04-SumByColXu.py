#!/usr/bin/python
##################################################################################
# Author: Huaqin Xu (huaxu@ucdavis.edu)
# Date: Aug.16 2013
# Description:
#
# This python script sum the occurence of 'A', 'B','-' for each scaffold by defined lines  
#
# =================================================================================
# input arguments:
#	1.input file.
#	2.number of lines for each scaffold to sum
#		  
# Output: sum files.
#
######################################################################################

import csv, os, sys
from os.path import basename, splitext
from itertools import islice
import timeit

######################################################

start = timeit.default_timer()

# ----- get options and file names and open files -----
if len(sys.argv) == 3:
    infile = sys.argv[1]
    cnt= int(sys.argv[2])
else: 
    print len(sys.argv)
    print 'Usage: [1]infile, [2]num of lines to sum'
    sys.exit(1) 

infbase = splitext(basename(infile))[0]
outfile = 'sum_' + infbase + '.tsv'

# ----- count the number of lines for each scaffold -------
tsvid = open(infile,'rb')
tsvidreader = csv.reader(tsvid, delimiter='\t')

idlist = []
first = 1
idcount = 0
for id in tsvidreader:
    if id == [] or id == '\n':
        continue
    if first == 1:
        curid = id[0]
        first = 0

    if id[0] != curid:
        idlist = idlist + [[curid, idcount]]
        idcount = 1
        curid = id[0]
    else:
        idcount = idcount + 1
idlist = idlist + [[curid, idcount]]

tsvid.close()

# ----- main: count the occurance of 'A','B','-' --------

tsvin = open(infile,'rb')
tsvout = open(outfile, 'wb')
tsvinreader = csv.reader(tsvin, delimiter='\t')
tsvoutwriter = csv.writer(tsvout, delimiter='\t')

cntA = 0
cntB = 0
cntN = 0

for alist in idlist:
    id = alist[0]
    scaflines = alist[1]
    print id + ":" + str(scaflines)
    for m in range(scaflines/cnt): # loop by each cnt-100 lines
        total = []
        cnt_lines = list(islice(tsvinreader, cnt))  # get cnt-100 lines
        pos = cnt_lines[0][1]
        for x in zip(*cnt_lines)[2:]:               # convert to list of columns and count
            cnt_N = list(x).count('-')
            cnt_A = list(x).count('A')
            cnt_B = list(x).count('B')
            cnt_U = list(x).count('U')
            total = total + [cnt_A, cnt_B, cnt_N, cnt_U]        
        tsvoutwriter.writerow([id]+ [pos] + total)

    if scaflines%cnt != 0: # get the rest of lines 
        total = []
        cnt_lines = list(islice(tsvinreader, scaflines%cnt))
        pos = cnt_lines[0][1]
        for x in zip(*cnt_lines)[2:]:
            cnt_N = list(x).count('-')
            cnt_A = list(x).count('A')
            cnt_B = list(x).count('B')
            cnt_U = list(x).count('U')
            total = total + [cnt_A, cnt_B, cnt_N, cnt_U]        
        tsvoutwriter.writerow([id]+ [pos] + total)
    tsvoutwriter.writerow([])
    tsvoutwriter.writerow([])    

stop = timeit.default_timer()
print stop - start 

