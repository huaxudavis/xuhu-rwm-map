#!/usr/bin/python
##################################################################################
# Author: Lutz Froenicke(lfroenicke@ucdavis.edu) and Huaqin Xu (huaxu@ucdavis.edu)
# Date: Aug.16 2013; last update: Aug.29 2013
# Description:
#
# This python script sum the occurence of 'A', 'B','-','U' for each scaffold by defined lines  
#
# =================================================================================
# input arguments:
#	1.input file.
#	2.if opt is l: it is the number of lines for each scaffold to sum, if num set to 0, scaffold will sum by ID name;
#         if opt is b, it is the position range to sum.
#       3.opt: l - sum by lines or b - sum by position range
#       4.format: f - first run or s - second run(has two extra blank lines between groups and two extra cols after pos)
#
# For example: python 03-SumByGroup.py x-haplo-split-sum.test.tab 0 l s
#  
# Output: sum files.
# Format: SNP group | start position | end position | # of SNPs in group | # of 'A' | # of 'B' | # of '-' | # of 'U' |......
#
######################################################################################

import csv, os, sys, timeit, math
from os.path import basename, splitext
from itertools import islice


######################################################
#count genotype
def genotypesum(x):
    cnt_N = list(x).count('-')
    cnt_A = list(x).count('A')
    cnt_B = list(x).count('B')
    cnt_U = list(x).count('U')
    return [cnt_A, cnt_B, cnt_N, cnt_U] 
    
######################################################

start = timeit.default_timer()

# ----- get options and file names and open files -----
if len(sys.argv) == 5:
    infile = sys.argv[1]
    cnt= int(sys.argv[2])
    opt = sys.argv[3]
    format = sys.argv[4]
else: 
    print len(sys.argv)
    print 'Usage: [1]infile, [2]num of lines/bases to sum, [3]opt: l or b, [4]format: f or s'
    sys.exit(1) 

if format == 's' and opt == 'b':
    print 'b option can not use with s option!'
    sys.exit(0)

infbase = splitext(basename(infile))[0]
outfile = 'sum_' + infbase + '.tsv'

# ----- count the number of lines for each scaffold -------
tsvid = open(infile,'rb')
tsvidreader = csv.reader(tsvid, delimiter='\t')

rowlist = []
first = 1
SNPcount = 0 # num of SNPs per scaffold or num of SNPs in position range
cutoff = cnt

for row in tsvidreader:
    if row == [] or row == '\n':
        continue
    if len(row)<2 and len(row)>0:
        continue
    if first == 1:
        curid = row[0]
        first = 0
        rowlen = len(row)
        
    if row[0] != curid:  # next scaffold
        rowlist = rowlist + [[curid, SNPcount]]
        SNPcount = 0
        curid = row[0]
        if opt == 'b':    
            cutoff = cnt

    if opt == 'l':
        SNPcount = SNPcount + 1     # count num of SNPs in scaffold group
    else:
        if int(row[1]) <= cutoff:
            SNPcount = SNPcount +1  # count num of SNPs in position range
        else: 
            rowlist = rowlist + [[curid, SNPcount]]
            cutoff = cutoff+cnt
            SNPcount = 0
            while int(row[1]) > cutoff:
                rowlist = rowlist + [[curid, SNPcount]]
                cutoff = cutoff+cnt                    
            SNPcount = 1  
                
rowlist = rowlist + [[curid, SNPcount]]

tsvid.close()

# ----- main: count the occurance of 'A','B','-' --------

tsvin = open(infile,'rb')
tsvout = open(outfile, 'wb')
tsvinreader = csv.reader(tsvin, delimiter='\t')
tsvoutwriter = csv.writer(tsvout, delimiter='\t')

if format == 'f': 
    s = 2   # id and position columns
else:
    s = 4   # id, start, end, and count columns

if opt == 'l':
    for alist in rowlist:
        id = alist[0]
        SNPlines = alist[1]
        print id + ":" + str(SNPlines)
        if cnt > 0:  #sum by the predefined number
            for m in range(SNPlines/cnt): # loop by each cnt-100 lines
                total = []
                cnt_lines = list(islice(tsvinreader, cnt))  # get cnt-100 lines
                posStart = cnt_lines[0][1]
                if format == 'f':
                    posEnd = cnt_lines[cnt-1][1]
                    sumcnt = cnt
                else:
                    posEnd = cnt_lines[cnt-1][2]
                    sumcnt = sum(int(z) for z in zip(*cnt_lines)[s-1])
                for x in zip(*cnt_lines)[s:]:            # convert to list of columns and count
                    total = total + genotypesum(x)      
                tsvoutwriter.writerow([id]+ [posStart] + [posEnd] + [sumcnt] + total)
            
            if SNPlines%cnt != 0: # get the rest of lines 
                total = []
                cnt_lines = list(islice(tsvinreader, SNPlines%cnt))
                posStart = cnt_lines[0][1]
                if format == 'f':
                    posEnd = cnt_lines[SNPlines%cnt-1][1]
                    sumcnt = SNPlines%cnt
                else:
                    posEnd = cnt_lines[SNPlines%cnt-1][2]
                    sumcnt = sum(int(z) for z in zip(*cnt_lines)[s-1])

                for x in zip(*cnt_lines)[s:]:
                    total = total + genotypesum(x)       
                tsvoutwriter.writerow([id]+ [posStart] + [posEnd] + [sumcnt] + total)
            tsvoutwriter.writerow([])
            tsvoutwriter.writerow([])
            if format != 'f':
                blank_lines = list(islice(tsvinreader, 2))
            
        else: # sum by each scaffold
            total = []
            cnt_lines = list(islice(tsvinreader, SNPlines))
            posStart = cnt_lines[0][1]
            if format == 'f':
                posEnd = cnt_lines[SNPlines-1][1]
                sumcnt = SNPlines
            else:
                posEnd = cnt_lines[SNPlines-1][2]
                sumcnt = sum(int(z) for z in zip(*cnt_lines)[s-1])

            for x in zip(*cnt_lines)[s:]:
                total = total + genotypesum(x)       
            tsvoutwriter.writerow([id]+ [posStart] + [posEnd] + [sumcnt] + total)
            if format != 'f':
                blank_lines = list(islice(tsvinreader, 2))
            
else:     
    id = 'first' # first scaffold

    for rlist in rowlist: # loop by position range for each scaffold
        if rlist[0] != id: # New scaffold
            posStart = 0
            posEnd = 0
            if id != 'first': # print two extra lines at the end of previous scaffold if it is not the first one
                tsvoutwriter.writerow([])
                tsvoutwriter.writerow([])
    
        id = rlist[0]
        SNPlines = rlist[1]
        print id + ":" + str(SNPlines)        
        total = []
        if SNPlines != 0:            
            cnt_lines = list(islice(tsvinreader, SNPlines))  # get lines by range
            posStart = (int(cnt_lines[0][1])/cnt)*cnt+1      
            posEnd = ((int(cnt_lines[SNPlines-1][1])/cnt)+1)*cnt
            sumcnt = SNPlines
            for x in zip(*cnt_lines)[s:]:               # convert to list of columns and count
                total = total + genotypesum(x)
                
        else:      # no data in this range
            posStart = posEnd +1
            posEnd = posEnd + cnt
            sumcnt = 0
            total = [0]*(rowlen-s)*4
        tsvoutwriter.writerow([id]+ [posStart] + [posEnd] + [sumcnt] + total)    

stop = timeit.default_timer()
print stop - start 

