#!/usr/bin/python
##################################################################################
# Atuthor: Lutz Froenicke(lfroenicke@ucdavis.edu) and Huaqin Xu (huaxu@ucdavis.edu)
# Date: Aug.16 2013, last update: Aug.29 2013
# Description:
#
# This python script split the SNP chunk files by the turning points 'A'<->'B'  
#
# =================================================================================
# input arguments:
#	1.input file.
#	2. number of samples
#       3. cutoff (for example: 25)
#		  
# Output: splited file.
#
######################################################################################

import csv, os, sys
from os.path import basename, splitext


######################################################

# ----- get options and file names and open files -----
if len(sys.argv) == 4:
    infile = sys.argv[1]
    samples = int(sys.argv[2])
    cutoff = int(sys.argv[3])
else: 
    print len(sys.argv)
    print 'Usage: [1]infile, [2]# of smaples, [3]cutoff'
    sys.exit(1) 

infbase = splitext(basename(infile))[0]
outfile = 'split_' + infbase + '.tsv'

tsvin = open(infile,'rb')
tsvout = open(outfile, 'wb')

tsvinreader = csv.reader(tsvin, delimiter='\t')
tsvoutwriter = csv.writer(tsvout, delimiter='\t')

# ----- main: count the occurance of 'A','B','-','U'--------

cntU = 0
start = 1
prerow =[]
urow = []
rowlist = []
scafcnt = 0
splitcnt = 0

for row in tsvinreader:
    if len(row) == 0 or row[0] == '':
        if rowlist != []: 
            for r in rowlist: # print the last part of the previous scaffold
                tsvoutwriter.writerow([curid] + r[1:])
            scafcnt = scafcnt + 1
            splitcnt = splitcnt + 1
        rowlist=[]
        tsvoutwriter.writerow([])
        start = 1
    else:
        id = row[0]
        mismatches = 0
        
        if start == 1:   # initial values if it is the first line of scaffold
            curid = id
            prerow = row[4:]
            start = 0
            sub = 1
            rowlist.append(row)

        else:
            cntU = row[4:].count('U') # count the occurence of 'U'
            mismatches = len([x for x, y in zip(prerow, row[4:]) if x != y]) # find mismatch items in previous and current rows         
            
            if cntU >= cutoff:
                urow = row
                row = tsvinreader.next()
                mismatches = len([x for x, y in zip(prerow, row[4:]) if x != y])
            prerow = row[4:]
     
            if mismatches > cutoff: # split the mismatch parts
                curid = id + '_' + str(sub)
                for r in rowlist: # print part of scaffold 
                    tsvoutwriter.writerow([curid] + r[1:])
                tsvoutwriter.writerow([])
                tsvoutwriter.writerow([])
                splitcnt = splitcnt + 1
                
                sub = sub + 1 # prepare for the next part of scaffold
                curid = id + '_' + str(sub)
                rowlist = []
                rowlist.append(row)                    
            else:
                if cntU >= cutoff: 
                    rowlist.append(urow)
                rowlist.append(row)
                    
                
for r in rowlist: # last part of last scaffold
    tsvoutwriter.writerow([curid] + r[1:])
scafcnt = scafcnt + 1
splitcnt = splitcnt + 1    

print "Processed Scaffolds:" + str(scafcnt)
print "New Scaffold segments:" + str(splitcnt)
