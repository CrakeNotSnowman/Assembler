#!/usr/bin/env python


import JBGenSuffTree
import KM_TreeAlign
import matplotlib.pyplot as plt
from time import gmtime, strftime, sleep
import gc
import sys

def jbGenSuffTreeAlign(infile, outfile):
    JBGenSuffTree.GenTreeAlign(infile, outfile)

def sfxTreeBasicAlign(filePath, overlapThreshold, outAlignedFragReads):

	#get files: send to dictionary 
	#  if already in dictionary, entry += 1
	gc.enable()

	outFile 	= str(outAlignedFragReads)

	rawfile 	= open(str(filePath), 'r')
	fragReadsA 	= {}				#dictionary
	#fragReadsB 	= {}				#dictionary
	fragList	= []
	unMerged 	= {}
	unMergedList 	= []
	startLen	= []
	finalLen 	= []
	repeats 	= 0
	stringA 	= ''
	stringAcount 	= 0
	stringB 	= ''
	stringBcount 	= 0
	dictSize 	= 0
	i 		= 0
	merge 		= False
	mergeCount 	= 0

	#print "File Name:\t\t", filePath

	while True:
	    line = rawfile.readline().strip()
	    if (line == ''):
		# end of file
		rawfile.close()
		break
	    if (line in fragReadsA) :
		# frag is a repeat
		fragReadsA[line] += 1
		repeats += 1
		startLen.append(len(line))
	    else :
		# fragment is new
		fragReadsA[line.upper()] = 1
		fragList.append(line.upper())
		dictSize += 1
		startLen.append(len(line))



	rawfile.close()

	#print "repeats: \t\t", repeats
	#print "fragments: \t\t", len(fragList)
	dictSize = dictSize -1


	#print stringA
	#print stringB
	#stringA = "TTGGCTATCCTTTATATTTTAAGGGTTATTAGGATATTTTTTATTATGACTACATGGGATAAATGTTTAAAAAAAATAAAAAAAAACCTTTCTACGTTTG"
	#dictSize = 3

	#one pass through, align string A
	unMergedCount = 0
	loopCount = 0
	merge = False

	while (len(fragList) > 0) :
	    #print "Loop Count", unMergedCount
	    loopCount += 1
	    i = 0
	    if (merge == False) :
		stringA = fragList.pop(0)

	    merge = False
	    for i in range(len(fragList) ) :
		#print "\t", i, len(fragList)
		#print fragList
		#print unMergedList
		stringB = fragList[i]
		alignVals = KM_TreeAlign.Align(stringA, len(stringA), stringB, len(stringB))
		#print "\t", alignVals[0], "\t", alignVals[1], "\t", alignVals[2], "\t", alignVals[3]
		if (alignVals[1] + 1 >= overlapThreshold):
		    #print "YES"
		    if (alignVals[0] + alignVals [1] + 1 == len(stringA) ):
			#stich
			stringA = ''.join([stringA, stringB[(alignVals [1] + 1 ):len(stringB)  ]])
			#print "OVERLAP", len(stringA)
			fragList.remove(stringB)
			#print dictSize
			merge = True
			mergeCount += 1
			i = 0
			break
		    elif  ((alignVals[1] +1 >= len(stringB)) and ( len(stringB) > 500)):
			# Check if B is a substring of A
			#print alignVals[1], len(stringB)
			fragList.remove(stringB)	
			merge = True
			mergeCount += 1
			i = 0
			break		
		elif (alignVals[3] + 1 >= overlapThreshold):
		    #print "ALSO"
		    if (alignVals[2] + alignVals [3] + 1 == len(stringB) ):
			#print "OVERLAP"
			#stich
			stringA = ''.join([stringB, stringA[(alignVals [3] + 1 ):len(stringA) ]])
			fragList.remove(stringB)
			#print dictSize
			merge = True
			mergeCount += 1
			i = 0
			break
		    elif  ((alignVals[3] +1  >= len(stringA)) and (len(stringA) > 500)):
			# Check if A is a substring of B
			#print alignVals[3], len(stringA)
			fragList.remove(stringA)
			stringA = stringB	
			merge = True
			mergeCount += 1
			i = 0
			break	
		
	    if (merge == False) :		
		unMergedList.append(stringA)
		unMergedCount += 1
	
	if (merge == True) :
	    unMergedList.append(stringA)

	#print "Overlap Threshold\t", overlapThreshold
	#print "Merge Count: \t\t", mergeCount


	i = 0
	with open(outFile,'a') as file_Write:
		for i in range(len(unMergedList) ):
			finalLen.append(len(unMergedList[i]))
			file_Write.write(str(unMergedList[i]))
			file_Write.write("\n")
	gc.disable()
	
	#plotAlignment(finalLen)
	return mergeCount, finalLen, startLen


