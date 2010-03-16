#!/usr/bin/python
##################################################################################
# Author: Huaqin Xu (huaxu@ucdavis.edu)
# Supervisor: Alexander Kozik (akozik@atgc.org)
# Date: Oct.12. 2009
# Last update: Mar.1. 2010
# Description:
#
# This python script splits file by prefix and count the match information. 
#
# =================================================================================
# input arguments:
#	1.input file.
#	2.windowsize
#	3.windowstep.
#		  
# Output: split files and log file.
#
######################################################################################

import sys
import re
import array
import os

from math import *
from os.path import exists, join, basename, splitext


# ================================ functions ==========================================
# ======================== Open and read file functions ===============================
def open_file(file_name, mode):
	if file_name == "":
		print 'Empty input file name!'
		raw_input("\nPress the enter key to exit.")
		sys.exit(0)

	try:
		the_file = open(file_name, mode)
	except(IOError), e:
		print "Unable to open the file", file_name, "Ending program.\n", e
		raw_input("\nPress the enter key to exit.")
		sys.exit(0)
	else:
		return the_file

def read_file(afile):
	try:
		flines = afile.readlines()
	except:
		print 'Failed to read from: ', afile
		sys.exit(0)
	else:
		return flines


#=============================== function getData  ===================================
# This function reads data from spp files
# 1. Get the header of the file -> varaible: header
# 2. Get different spps 4 letter codes -> array: queue 
# 3. For each 4_letter_code group, get the spps and seqences -> dictionary: content
# 4. Put all spps, get seqences -> dictionary: seq

def getData(lines):
	#------------------------ variables -----------------------------------------
	global header				# header of the input file
	global queue				# spps groups
	global content				# spps and seqs stored by each group
	global seq				# seqs for all spps


	prefix=4				# prefix of spp names are 4 letters
	delimiter="\t"
	datalen=len(lines)
	if(datalen == 0):
		print "Empty Data file!"
		sys.exit(0)
		
	colcount=lines[0].count(delimiter)+1 	# how many columns in the file
	header = lines[0] 			# get the header of the file
	
	queue=[]				
	content={}				
	seq={}
	aq = 'XXXX'				# 4_letter_code group name
	
	#--------------- read lines from file ----------------------------------------
	for l in range(1, datalen):

		# check for empty lines and incorrect field numbers
		if lines[l] != '\n':
			if lines[l].count(delimiter)==colcount-1:
				arow=lines[l].rstrip().split(delimiter)
				if(arow[0][0:prefix] != aq):		# find a differnt group
					aq = arow[0][0:prefix]
					queue.append(aq)		# append the different group to queue
					content[aq] = {}		# initialize the content[groupname]
					print aq			# debugging
				
				content[aq][arow[0]] = arow[1:] 	# content[groupname][sppname] = seq
				seq[arow[0]] = arow[1:]			# seq[sppname] = seq
			
			else:
				print "Error: Line #%s has inconsistent number of columns.\n" %(l+1)
				sys.exit(0)
		else:
			print "Skip an empty line at #%s.\n" %(l+1)        


#=========================== function findgroup ======================================
# This function will find spps groups, generate consensus seq and print output file
# This function will call function getconsensus and printtable

def findgroup(q):
	
	# ---------------------------- variables  -----------------------------------
	global windowsize
	global windowstep
	global content
	global header
	global seq

	misscutoff = 180 			# Cut off value total scores
	radiolistA={}				# store the radio of 'A' for each spp
	groupx={}				# store spps in each group
	
	#=== calculate the minimum allow distance ===
	maxgroupnum = int((1-windowsize)/windowstep)+1
	i=1
	mad=0					# min allow distence
	while(windowstep*i < windowsize and i<maxgroupnum):
		mad = i
		i=i+1

	#=== initialize all groups ===		
	for i in range(0, maxgroupnum+1):
		groupx[i] = []		

	
	#------------------- calculate radioA and locate spp into the group -----------
	for (id, line) in content[q].items():
		unisum = dict(map(lambda k:(k,1), line)).keys() # get unique letters in the spp
		unisum.sort()
		countsum = map(lambda t: (float(line.count(t))/float(len(line))), unisum) # count No. of each letter
		
		#---------- calculate A/(A+B) ratio ----------
		if(len(unisum) == 1):
			if(unisum[0]== 'A'):
				radiolistA[id] = 1.0
			else:
				radiolistA[id] = 0
		elif(len(unisum) == 2):
			if(unisum[0]!= '-'):
				radiolistA[id] = countsum[0]/(countsum[0]+countsum[1])
			elif(unisum[0]== '-' and countsum[0] < misscutoff):
				if(unisum[1]== 'A'):
					radiolistA[id] = 1.0
				else:
					radiolistA[id] = 0
			else:
				radiolistA[id] = 0
		else:
			if(countsum[0] < misscutoff):
				radiolistA[id] = countsum[1]/(countsum[1]+countsum[2])
			else:
				radiolistA[id] = 0
						
		#---------- locate spps into group ----------
		for i in range(0, maxgroupnum+1):
			if(radiolistA[id] >= windowstep*i and radiolistA[id] <= windowsize+windowstep*i):
				groupx[i].append(id)
	
	#------------------- find the first two largest groups ------------------------

	max1 = 0 	# No. of spps in the first large group
	order1 = 0	# Group No. of the first large group 
	subgroup1 ={}	# spps and seqs of the first large group
	max2 = 0 	# No. of spps in the second large group
	order2 = 0	# Group No. of the second large group 
	subgroup2 ={}	# spps and seqs of the second large group 
	
	for(group,item) in groupx.items():
		l = len(item)
		if( l > max1):
			max1 = l
			order1 = group
		elif(l > max2):
			max2 = l
			order2 = group
	
	# ---------------- print each group and get consensus --------------------

	if(abs(order1 - order2) <= mad or max2 == 0):	# have one group
		suffix = '0'			# suffix for subgroup
		for id in groupx[order1]:	# get spps and seqs in the group
			subgroup1[id] = seq[id]
		spp = id[0:4]
		radioA = str(windowstep*order1)+"-"+str(windowsize+windowstep*order1) # calculate radioA range for the group
		matchseq1 = getconsensus(subgroup1)	# get consensus   
		printtable(spp, suffix, subgroup1, matchseq1, radioA)	# print output
	else:					#have two groups
		suffix = '1'			
		for id in groupx[order1]:	
			subgroup1[id] = seq[id]
		spp = id[0:4]
		radioA = str(windowstep*order1)+"-"+str(windowsize+windowstep*order1)	
		matchseq1 = getconsensus(subgroup1)  
		printtable(spp, suffix, subgroup1, matchseq1, radioA)
		
		suffix = '2'			
		for id in groupx[order2]:	
			subgroup2[id] = seq[id]
		radioA = str(windowstep*order2)+"-"+str(windowsize+windowstep*order1)	
		matchseq2 = getconsensus(subgroup2)  
		printtable(spp, suffix, subgroup2, matchseq2, radioA)

	
		
#=========================== function getconsensus =====================================
# Input: a set of spp seqs
# Function: gets the number of mismatches and generates the consensus seq
# Output:  consensus seq and No. of mismatches

def getconsensus(groupseq):
	global header
	
	totallen = header.count("\t")
	matchseq = []
	mismatch = 0
	
	for j in range(0, totallen):
		sum = map(lambda i:i[j], groupseq.values())  # get all letters in j position
		unisum = dict(map(lambda k:(k,1), sum)).keys() # get unique letters in j position
		unisum.sort()
		
		# get the number of mismatch in the group
		if(len(unisum)==2 and unisum[0] != '-') or (len(unisum)>2):
			mismatch = mismatch +1
			
		# get No. of each letter	
		countsum = map(lambda t: (float(sum.count(t))/float(len(sum))), unisum)


		# get the consensus seq
		if(len(unisum)==1): 
			matchseq.append(unisum[0])	                   
		elif(len(unisum)==2):
			if(countsum[0] >= countsum[1]):
				matchseq.append(unisum[0])
			else:
				matchseq.append(unisum[1])
		else:
			if(countsum[0] >= 0.5):
				matchseq.append(unisum[0])
			elif(countsum[1] >= countsum[2]):
				matchseq.append(unisum[1])
			else:
				matchseq.append(unisum[2])
				
	matchseq.append(mismatch) # append No. of mismatches at the end of consensus seq
	return matchseq


#=========================== function printtable =====================================
# Input: spp: 4_letter_code,
#	 suffix: suffix of the subgroup
#	 subgroup: spp and seqs of this subgroup
#	 matchseq: consensus seq wtih number of mismatch 
#	 radioA: radioA range
# Function: write to outfile and logfile

def printtable(spp, suffix, subgroup, matchseq, radioA):
	global header
	global totalspps
	global totalmismatch
	global logf
	global option	
	
	maxmismatch = 10	# max No. of mismatches to generate consensus seq
	
	#---------------- write output file --------------------	
	outfile= spp + "_" + suffix +".out"
	outf=open_file(outfile,'w')
		
	outf.write(header)			# write header
	for (key,line) in subgroup.items():	# write spps
		outf.write(key+"\t"+"\t".join(line)+"\n")
	
	mismatch = matchseq.pop()
	if(mismatch < maxmismatch): 		# write consensus seq
		consensus_success = "_Y_"
		outf.write("CONSENSUS_"+spp+"_"+suffix+"\t"+"\t".join(matchseq))
	else:
		consensus_success = "_N_"
	outf.close()
	
	#---------------- write log file --------------------
	logf.write(spp+"\t"+suffix+"\t"+str(totalspps)+"\t"+str(totalmismatch)+"\t"+str(len(subgroup))+"\t"+str(mismatch)+"\t"+consensus_success+"\t"+radioA + "\n")

#----------------------------- main ------------------------------------------------------

# ----- get options and file names and open files -----
if len(sys.argv) == 4:
	infile=sys.argv[1]
	windowsize = float(sys.argv[2])		# Ratio distance between groups
	windowstep = float(sys.argv[3])		# sliding window size
elif len(sys.argv) == 2:
	infile=sys.argv[1]
	windowsize = 0.1 			# default windowsize			
	windowstep = 0.05			# default windowstep
else:
	print len(sys.argv)
	print 'Usage: [1]infile ([2]windowsize [3]windowstep: optional)'
	sys.exit(1)
	
if(windowsize < windowstep):
	print 'Error: windowsize < windowstep'
	sys.exit(1)

# ----- read infile ---------------------------
inf=open_file(infile,'r')
inlines = read_file(inf)
infilebase = splitext(basename(infile))[0]
getData(inlines)

# ----- open log file to write -----------------
logfile=infilebase+".log"
logf =open_file(logfile,'w')
logf.write("spp\tgroup_id\ttotal_spps\ttotal_mismatch\tspps_in_subgroup\tmismatch_in_subgroup\tconsensus_success\tradio_group\n")

# ------ generate output file -------------------	
for q in queue:
	totalspps = len(content[q])			# total No. of spp
	mismatchseq = getconsensus(content[q])		# get total mismatches.
	totalmismatch = mismatchseq.pop()
	findgroup(q)					# find group, write output file and log file
logf.close()

print 'Please find log file: '+ logfile + ".\n"
#-------------------------- end of the program ----------------------------------------

		
		