#!/usr/bin/python
##################################################################################
# Author: Huaqin Xu (huaxu@ucdavis.edu)
# Supervisor: Alexander Kozik (akozik@atgc.org)
# Date: Oct.12. 2009
# Last update: April.20. 2010
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
	global colcount
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
	header = ""			# get the header of the file
	
	queue=[]				
	content={}				
	seq={}
	aq = 'XXXX'				# 4_letter_code group name
	
	#--------------- read lines from file ----------------------------------------
	for l in range(0, datalen):

		# check for empty lines and incorrect field numbers
		if lines[l] != '\n':
			if lines[l][0:1] == ";":
				header += lines[l]
			elif lines[l].count(delimiter)==colcount-1:
				arow=lines[l].rstrip().split(delimiter)
				if(arow[0][0:prefix] != aq):		# find a differnt group
					aq = arow[0][0:prefix]
					queue.append(aq)		# append the different group to queue
					content[aq] = {}		# initialize the content[groupname]
					#print aq			# debugging
				
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
	ratiolistA={}				# store the ratio of 'A' for each spp
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

	
	#------------------- calculate ratioA and locate spp into the group -----------
	for (id, line) in content[q].items():
		unisum = dict(map(lambda k:(k,1), line)).keys() # get unique letters in the spp
		unisum.sort()
		countsum = map(lambda t: (float(line.count(t))/float(len(line))), unisum) # count No. of each letter
		
		#---------- calculate A/(A+B) ratio ----------
		if(len(unisum) == 1):
			if(unisum[0]== 'A'):
				ratiolistA[id] = 1.0
			else:
				ratiolistA[id] = 0
		elif(len(unisum) == 2):
			if(unisum[0]!= '-'):
				ratiolistA[id] = countsum[0]/(countsum[0]+countsum[1])
			elif(unisum[0]== '-' and countsum[0] < misscutoff):
				if(unisum[1]== 'A'):
					ratiolistA[id] = 1.0
				else:
					ratiolistA[id] = 0
			else:
				ratiolistA[id] = 0
		else:
			if(countsum[0] < misscutoff):
				ratiolistA[id] = countsum[1]/(countsum[1]+countsum[2])
			else:
				ratiolistA[id] = 0
						
		#---------- locate spps into group ----------
		for i in range(0, maxgroupnum+1):
			if(ratiolistA[id] >= windowstep*i and ratiolistA[id] <= windowsize+windowstep*i):
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
		ratioA = str(windowstep*order1)+"-"+str(windowsize+windowstep*order1) # calculate ratioA range for the group
		matchseq1 = getconsensus(spp, suffix, subgroup1)	# get consensus   
		printtable(spp, suffix, subgroup1, matchseq1, ratioA)	# print output
	else:					#have two groups
		suffix = '1'			
		for id in groupx[order1]:	
			subgroup1[id] = seq[id]
		spp = id[0:4]
		ratioA = str(windowstep*order1)+"-"+str(windowsize+windowstep*order1)	
		matchseq1 = getconsensus(spp, suffix, subgroup1)  
		printtable(spp, suffix, subgroup1, matchseq1, ratioA)
		
		suffix = '2'			
		for id in groupx[order2]:	
			subgroup2[id] = seq[id]
		ratioA = str(windowstep*order2)+"-"+str(windowsize+windowstep*order1)	
		matchseq2 = getconsensus(spp, suffix, subgroup2)  
		printtable(spp, suffix, subgroup2, matchseq2, ratioA)

	
		
#=========================== function getconsensus =====================================
# Input: a set of spp seqs
# Function: gets the number of mismatches and generates the consensus seq
# Output:  consensus seq and No. of mismatches

def getconsensus(spp, suffix, groupseq):
	global colcount
	
	matchseq = []
	countstr = ""
	consensusstr = ""
	matchresult = []
	countseq={}
	fractseq={}

	mismatch = 0
	countseq['A']={}
	countseq['B']={}
	countseq['M']={}
	countseq['AB']={}
	countseq['ALL']={}
	fractseq['A']={}
	fractseq['B']={}
	fractseq['M']={}
	
	for j in range(0, colcount-1):
		sum = map(lambda i:i[j], groupseq.values())  # get all letters in j position
		unisum = dict(map(lambda k:(k,1), sum)).keys() # get unique letters in j position
		unisum.sort()
		
		# get the number of mismatch in the group
		if(len(unisum)==2 and unisum[0] != '-') or (len(unisum)>2):
			mismatch = mismatch +1
			
		# get No. of each letter	
		letter = ['A', 'B', '-']
		countsum = map(lambda t: (sum.count(t)), letter)		
		
		countseq['A'][j] = str(countsum[0])
		countseq['B'][j] = str(countsum[1])
		countseq['M'][j] = str(countsum[2])
		countseq['AB'][j] = str(countsum[0]+countsum[1])
		countseq['ALL'][j] = str(countsum[0]+countsum[1]+countsum[2])
		print spp+"\t"+str(countsum[0])+"\t"+str(countsum[1])+"\t"+str(countsum[2])		
		
		if(countseq['AB'][j] == "0"):
			fractseq['A'][j] = "0"
			fractseq['B'][j] = "0"
		else:
			fractseq['A'][j] = str(int(round((float(countsum[0])/float(countsum[0]+countsum[1]))*100)))
			fractseq['B'][j] = str(int(round((float(countsum[1])/float(countsum[0]+countsum[1]))*100)))
		fractseq['M'][j] = str(int(round((float(countsum[2])/float(countsum[0]+countsum[1]+countsum[2]))*100)))	
		

		
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
				
	countstr += "COUNT_A_"+spp+"_"+suffix+"\t"+"\t".join(countseq['A'].values()) + "\n"
	countstr += "COUNT_B_"+spp+"_"+suffix+"\t"+"\t".join(countseq['B'].values()) + "\n"
	countstr += "COUNT_M_"+spp+"_"+suffix+"\t"+"\t".join(countseq['M'].values()) + "\n"
	countstr += "COUNT_AB_"+spp+"_"+suffix+"\t"+"\t".join(countseq['AB'].values()) + "\n"
	countstr += "COUNT_ALL_"+spp+"_"+suffix+"\t"+"\t".join(countseq['ALL'].values()) + "\n"
	countstr += "FRACT_A_"+spp+"_"+suffix+"\t"+"\t".join(fractseq['A'].values()) + "\n"
	countstr += "FRACT_B_"+spp+"_"+suffix+"\t"+"\t".join(fractseq['B'].values()) + "\n"
	countstr += "FRACT_M_"+spp+"_"+suffix+"\t"+"\t".join(fractseq['M'].values()) + "\n"
	
	consensusstr = "CONSENSUS_"+ spp+"_"+ str(len(groupseq))+ "_" +suffix+"\t"+"\t".join(matchseq)+ "\n"
	
	matchresult.append(countstr) # count informaiton
	matchresult.append(consensusstr) # consensus sequence
	matchresult.append(mismatch) # No. of mismatches
	return matchresult


#=========================== function printtable =====================================
# Input: spp: 4_letter_code,
#	 suffix: suffix of the subgroup
#	 subgroup: spp and seqs of this subgroup
#	 matchresult: consensus seq wtih number of mismatch 
#	 ratioA: ratioA range
# Function: write to outfile and logfile

def printtable(spp, suffix, subgroup, matchresult, ratioA):
	global header
	global totalspps
	global totalmismatch
	global logf
	global option
	global maxmismatch
	
	#---------------- write output file --------------------	
	outfile= spp + "_" + suffix +".out"
	outf=open_file(outfile,'w')
		
	outf.write(header)			# write header
	for (key,line) in subgroup.items():	# write spps
		outf.write(key+"\t"+"\t".join(line)+"\n")
	outf.write(matchresult[0])
	
	mismatch = matchresult.pop()
	if(mismatch < maxmismatch): 		# write consensus seq
		consensus_success = "_Y_"
		outf.write(matchresult[1]) # CONSENSUS_XXXX_
	else:
		consensus_success = "_N_"
	outf.close()
	
	#---------------- write log file --------------------
	logf.write(spp+"\t"+suffix+"\t"+str(totalspps)+"\t"+str(totalmismatch)+"\t"+str(len(subgroup))+"\t"+str(mismatch)+"\t"+consensus_success+"\t"+ratioA + "\n")

#----------------------------- main ------------------------------------------------------

# ----- get options and file names and open files -----
if len(sys.argv) == 5:
	infile=sys.argv[1]
	windowsize = float(sys.argv[2])		# Ratio distance between groups
	windowstep = float(sys.argv[3])		# sliding window size
	maxmismatch = int(sys.argv[4])
elif len(sys.argv) == 2:
	infile=sys.argv[1]
	windowsize = 0.1 			# default windowsize			
	windowstep = 0.05			# default windowstep
	maxmismatch = 10			# max No. of mismatches to generate consensus seq
else:
	print len(sys.argv)
	print 'Usage: [1]infile ([2]window size [3]window step [4]max mismatch: optional)'
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
logf.write("spp\tgroup_id\ttotal_spps\ttotal_mismatch\tspps_in_subgroup\tmismatch_in_subgroup\tconsensus_success\tratio_group\n")

# ------ generate output file -------------------	
for q in queue:
	print q
	totalspps = len(content[q])			# total No. of spp
	mismatchseq = getconsensus(q, '0', content[q])		# get total mismatches.
	totalmismatch = mismatchseq.pop()
	findgroup(q)					# find group, write output file and log file
logf.close()

print 'Please find log file: '+ logfile + ".\n"
#-------------------------- end of the program ----------------------------------------

		
		