#!/usr/bin/python
##################################################################################
# Author: Lutz Froenicke(lfroenicke@ucdavis.edu) and Huaqin Xu (huaxu@ucdavis.edu)
# Date: Aug.16 2013, last update: Nov.4 2013
# Description:
#
# This python script split the SNP chunk files by the turning points 'A'<->'B'
#
# =================================================================================
# input arguments:
#   1.input file.
#   2. number of samples
#   3. cutoff for mismatches (ratio of mismatches/matches between A and B genotype callls - please provide a value in the range > 0 and < 1 )
#   4. cutoff for 'U' (for example: 25)
#
# Output: split file.

# usage: python 05-ScaffoldSplit.py genotypetable.tsv 88 0.1 10
#
######################################################################################

import csv, os, sys
from os.path import basename, splitext

######################################################

# ----- get options and file names and open files -----
if len(sys.argv) == 5:
    infile = sys.argv[1]
    samples = int(sys.argv[2])
    cutoff = float(sys.argv[3])
    cutoffU = int(sys.argv[4])
#    cutoffM = int(sys.argv[5])
else:
    print len(sys.argv)
    print "Usage: [1]infile, [2]# of samples, [3]cutoff of missmatches, [4]cutoff of 'U'" #, [5]cutoff of '-'"
    sys.exit(1)

infbase = splitext(basename(infile))[0]
outfile = 'split_AB_' + sys.argv[3] + '-' + sys.argv[4] + infbase + '.tsv'  #  + '-' + sys.argv[5]

tsvin = open(infile,'rb')
tsvout = open(outfile, 'wb')
log = open("log_split_AB_" + sys.argv[3] + '-' + sys.argv[4] + infbase + '.txt','w') # + '-' + sys.argv[5]

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
                if(len(row)!= ''):
                    tsvoutwriter.writerow([curid] + r[1:])
                else:
                    tsvoutwriter.writerow(r)
            scafcnt = scafcnt + 1
            splitcnt = splitcnt + 1
        rowlist=[]
        tsvoutwriter.writerow([])
        start = 1
    else:
        id = row[0]
        matches = 0
        mismatches = 0

        if start == 1:   # initial values if it is the first line of scaffold
            curid = id
            prerow = row

            start = 0
            sub = 1
            rowlist.append(row)

        else:
            cntU = row[4:].count('U') # count the occurence of 'U'
        #    cntM = row[4:].count('-') # count the occurence of '-'
            matches = len([x for x, y in zip(prerow[4:], row[4:]) if ( x == 'A' and y == 'A' ) or ( x == 'B' and y == 'B' ) or ( x == 'U' and y == 'U' )]) # find match items in previous and current rows
            mismatches = len([x for x, y in zip(prerow[4:], row[4:]) if ( x == 'A' and y == 'B' ) or ( x == 'B' and y == 'A' )]) # find mismatch items in previous and current rows

            if cntU >= cutoffU:
                urow = row
                row = tsvinreader.next()
                if len(row) != 0 and row[0] != '':
                    matches = len([x for x, y in zip(prerow[4:], row[4:]) if ( x == 'A' and y == 'A' ) or ( x == 'B' and y == 'B' ) or ( x == 'U' and y == 'U' )])
                    mismatches = len([x for x, y in zip(prerow[4:], row[4:]) if ( x == 'A' and y == 'B' ) or ( x == 'B' and y == 'A' )])

            if  float(mismatches)/(matches+mismatches) > cutoff and len(row) != 0 and row[0] != '': # split the mismatch parts
                curid = id + '__' + str(sub)
                log.write(id + "\t" + str(prerow[2]) + '\n') # write the split position to log file
                for r in rowlist: # print part of scaffold
                    tsvoutwriter.writerow([curid] + r[1:])

                breakline = ['']+ [''] +['I']*(samples+2)
                tsvoutwriter.writerow(breakline)
                tsvoutwriter.writerow(breakline)
                splitcnt = splitcnt + 1

                sub = sub + 1 # prepare for the next part of scaffold
                curid = id + '__' + str(sub)
                rowlist = []
                rowlist.append(row)
            else:
                if cntU >= cutoffU:
                    rowlist.append(urow)
                rowlist.append(row)
            prerow = row

for r in rowlist: # last part of last scaffold
    tsvoutwriter.writerow([curid] + r[1:])
scafcnt = scafcnt + 1
splitcnt = splitcnt + 1

print "Processed Scaffolds:" + str(scafcnt)
print "New Scaffold segments:" + str(splitcnt)

log.write("\n\n")
log.write("Filename:" + sys.argv[1] + '\n')
log.write("Number Of Samples:" + sys.argv[2] + '\n')
log.write("Cutoff Threshold:" + sys.argv[3] + '-' + sys.argv[4] +'\n')  # + '-' + sys.argv[5]
log.write("Processed Scaffolds:" + str(scafcnt) + '\n')
log.write("New Scaffold segments:" + str(splitcnt)+ '\n')
