#!/usr/bin/python
import csv, os, sys
from os.path import basename, splitext
##################################################################################
# Author: Huaqin Xu (huaxu@ucdavis.edu)
# Date: Aug.15 2013
# Description:
#
# This python script sum the occurence of [ATCG|atcg] or [,.] for each scaffold,
# and also convert count result to genotype table
# and the genotype table can be cleaned if the option is 1,
# =================================================================================
# input arguments:
#	1.input file.
#	2.number of sample
#       3.opt: cleaned (1) or not cleaned(0)
#		  
# Output: count files and genotype files, if option is 1: cleaned genotype file and cleanCount file.
#
######################################################################################

# count in cells
def countcell(cnt, pu):
    Vcnt = pu.count("a") + pu.count("A") + pu.count("c") + pu.count("C") +pu.count("t") + pu.count("T") + pu.count("g") + pu.count("G")
    Ccnt = pu.count(",") + pu.count(".")
    if int(cnt) != Vcnt + Ccnt:
        Vcnt = 0
        Ccnt = 0
    C = [Vcnt, Ccnt]    
    return C

#convert count number to genotype
def genotype(C):    
    Vcnt = float(C[0])
    Ccnt = float(C[1])

    GT = ""
    if Vcnt == Ccnt == 0:
        GT = "-"
    elif Vcnt == 0:
        GT = "A"
    elif Ccnt == 0:
        GT = "B"
    elif 0.2 <= Vcnt / Ccnt <= 5:
        GT = "U"
    elif Vcnt / Ccnt < 0.2:
        GT = "A"
    elif Vcnt / Ccnt > 5:
        GT = "B"    
    return [GT]

#clean up genotype table by genotype count
def genotypecount(row, samples):
    U = row.count("U")
    A = row.count("A")
    B = row.count("B")
    M = row.count("-")
    R = ""

    F = ""
    if M > 0.8 * samples :
        F = F + "M"
    if U > 0.05 * samples:
        F = F + "U"
    if A == 0 or B == 0:
        F = F + "Y"
    if A != 0 and B != 0 :
        R = round((A/float(B)), 2)
    if R < 0.1 or R > 9:
        F = F + "X"

    C = [U, A, B, M, R, F]
    return C

######################################################

# ----- get options and file names and open files -----
if len(sys.argv) == 4:
    infile = sys.argv[1]
    samples = int(sys.argv[2])
    opt = int(sys.argv[3])
else: 
    print len(sys.argv)
    print 'Usage: [1]infile, [2]num of samples, [3]opt: 0 or 1'
    sys.exit(1)

infbase = splitext(basename(infile))[0]
outfile1 = 'parsed_' + infbase + '.tsv'
outfile2 = 'genotyped_' + infbase + '.tsv'

tsvin = open(infile,'rb')
tsvout1 = open(outfile1, 'wb')
tsvout2 = open(outfile2, 'wb')

tsvinreader = csv.reader(tsvin, delimiter='\t')
tsvoutwriter1 = csv.writer(tsvout1, delimiter='\t')
tsvoutwriter2 = csv.writer(tsvout2, delimiter='\t')

if opt == 1:
    outfile3 = 'cleanedGTs_' + infbase + '.tsv'
    outfile4 = 'GTcounts_' + infbase + '.tsv'
    tsvout3 = open(outfile3, 'wb')
    tsvout4 = open(outfile4, 'wb')
    tsvoutwriter3 = csv.writer(tsvout3, delimiter='\t')
    tsvoutwriter4 = csv.writer(tsvout4, delimiter='\t')

# ------ loop through each line to count ocurrence -------

for row in tsvinreader:
    cntresult = [row[0]] + [row[1]]
    gtresult = [row[0]] + [row[1]]
    for x in range(1,samples+1):
        Cnt = countcell(row[x*3] , row[x*3+1])
        cntresult = cntresult + Cnt   # countcell has two input values: cnt and teh cell with the genotypes
        gtresult = gtresult + genotype(Cnt)
    tsvoutwriter1.writerow(cntresult)
    tsvoutwriter2.writerow(gtresult)
    
    if opt == 1:
        gtcnt = [row[0]] + [row[1]]
        gtcnt = gtcnt + genotypecount(gtresult, samples) #get genotype count
        tsvoutwriter4.writerow(gtcnt)
        if gtcnt[7] == '':             # write to cleanGT file if genotype count is good
            tsvoutwriter3.writerow(gtresult)


    
