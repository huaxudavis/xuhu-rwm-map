#!/usr/bin/python
import csv, os, sys
from os.path import basename, splitext

#####################################################
''' 
This script converts the group counts into group haplotypes with the GroupHaplotyper script.
The haplotyping thresholds are hardcoded in the script and should be edited there at the top of the script.
If the B/A call ratio is between 0.2 and 5, the group is typed as heterozygous;
If it is lower than 0.2 the haplotype is set to A; if it is bigger than 5, the haplotype is set to B.
If the rate of missing data is greater, the 90% the group is typed as missing data.
If the number of het calls is equal or greater than the sum of A calls and B calls, the haplotype is set to heterozygous.

Usage: python  04-GroupHaplotyper.py  chunked-clean-genotype-table.tsv  number-of-samples &
e.g:   python  04-GroupHaplotyper.py  genotype-table.tsv  88 &
'''
#####################################################
# B/A ratio between "min" and "max" --> resulting call: heterozygous
# B/A ratio between smaller than or equal to "min" --> resulting call: A
# B/A ratio between bigger than or equal to "max" --> resulting call: B

min = 0.2
max = 5

######################################################
#count in cells
def genotype(groupsize, Vcnt, Ccnt, Mcnt, Ucnt):
    Vcnt = float(Vcnt)
    Ccnt = float(Ccnt)
    Mcnt = float(Mcnt)
    Ucnt = float(Ucnt)
    gs = int(groupsize)
	
    GT = ""
    if Vcnt == Ccnt == 0:
        GT = "-"
    elif Vcnt == 0:
        GT = "A"
    elif Ccnt == 0:
        GT = "B"
    elif min <= Vcnt / Ccnt <= max:    # thresholds for typing as heterozygous
        GT = "U"
    elif Vcnt / Ccnt <  min:
        GT = "A"
    elif Vcnt / Ccnt > max:
        GT = "B"
    if Vcnt < Ucnt and Ccnt < Ucnt:      # previously: "if Vcnt + Ccnt < Ucnt:"  or  "if Vcnt + Ccnt <= Ucnt:"  
        GT = "U"
    if Mcnt > 0.9 * gs:
        GT = "-"
    GTL = [GT]
    return GTL
######################################################

infile = sys.argv[1]
samples = int(sys.argv[2])
infbase = splitext(basename(infile))[0]
outfile = 'hpt-' + infbase + '.tsv'

tsvin = open(infile,'rb')
tsvout = open(outfile, 'wb')
tsvinreader = csv.reader(tsvin, delimiter='\t')
tsvoutwriter = csv.writer(tsvout, delimiter='\t')

for row in tsvinreader:
    if row == []:
        gtresult = row
    else:
        gtresult = [row[0]] + [row[1]] + [row[2]] + [row[3]]
        for x in range(1,samples + 1):
            gtresult = gtresult + genotype(row[3], row[x*4+1], row[x*4], row[x*4+2], row[x*4+3])   #
    tsvoutwriter.writerow(gtresult)
