#!/usr/bin/env python
'''
Keith Murray

Sequence Alignment Program:
  Mr. C's v1.22
  I need a better name for this




This program is intended to take fragment reads and assemble them 
over multiple itterations.

TODO:

STARTED: make text/email UI cleaner/not a text bomb
STARTED: Let user choose which emails to send
*  in email: check size of file, if too large, don't send
* comment program so it's not a mess
* add a garbage collection: Delete internal files
* add kmer option: open to user
* Consider making program keep a runtime log: what would go in it?
* Is stats page good enough?
* Is A_fspecs = "A_Output_Vectors/file_Specs" necessary? Can it be deleted
* plot 36 against source file
* Add notation explaining what is going on, my methods, motivations



TOTEST:
* AMI
* non-default
    Outfile



'''
import os.path
import sys
import argparse, textwrap
import math
import gc
import timeit
import datetime
import numpy
from time import gmtime, strftime

import A_AMIGenVects
import A_fragmentCluster
import B_ClusterAll
import alignMain
import sendMessage

# User Dependent Defaults
inFragReads_Default 		= 'KM_Hold/inputFile.txt'
outAlignedFragReads_Default 	= 'Outputs/OutputAlignments/outputPass0.txt'
metric_Default 			= 'N'
refGenomeFile_Default 		= 'KM_Hold/referenceFasta.txt'
window_Default 			= 25
kmerSize_Default 		= 7
numOfClusters_Default 		= 128
clusterDistMet_Default		= 'E'
warning_Default 		= False
AlignmentMethod_Default 	= 0
STOverlapThreshold_Default 	= 30
userFreqContact_Default 	= 'X@vtext.com'
userEmail_Default 		= 'X@gmail.com'
recurse_Default			= 0
stageText_Default 		= False
completeText_Default 		= False

# Internal Defaults 
clustCentFile = "KM_Hold/clustCents"
A_vectorFile = "A_Output_Vectors/Vector_File"
A_fspecs = "A_Output_Vectors/file_Specs"
prioriFile = "KM_Hold/priori"
clusterDistFolder = "B_ClusterDists/"
clusterFolder = "B_Clusters/"
statsFile = "Outputs/StatsFiles/StatsRun.txt"

# User values
toaddrsSMS  = 'X@vtext.com' 
toaddrsMMS  = 'X@vzwpix.com'


def interface():
    args = argparse.ArgumentParser(
        prog='KM_Assembler01.py', 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Assemble a genetic sequence from the input fragreads.',
        epilog=textwrap.dedent('''\
            * Fragread input file format : 
		Each line should be a new read, consisting of an 
		  aplhabet of A,C,G,T, and if necessary, N

              Example Input file:
		GAGACAT
		TAGGACTATTCANA
		ATTATAAGCCAG
		CATAGAG
		ATGTCTC
		
		

        '''))
    args.add_argument('-i', '--input-file', default= inFragReads_Default, help='Input Fragread file*')
    args.add_argument('-o', '--output-file', default= outAlignedFragReads_Default, help='Output Fragread filename')
    args.add_argument('-m', '--cluster-metric', default= metric_Default , help='Method to generate clusters: "D": Dictionary, "A" AMI')
    args.add_argument('-rf', '--reference-file', default= refGenomeFile_Default, help='Reference FASTA file for Dictionary Method')
    args.add_argument('-sss', '--stats', default= statsFile, help='Does Nothing')
    args.add_argument('-aw', '--ami-window', type=int, default= window_Default , help='Window Size for AMI Method')
    args.add_argument('-ks', '--kmer-size', type=int, default= kmerSize_Default , help='k value for kmer Method')
    args.add_argument('-r', '--recurse', type=int, default= recurse_Default , help='How many times the program should re-align the output file: stored in a sepperate file')
    args.add_argument('-c', '--cluster-count', type=int, default= numOfClusters_Default, help='Number of Clusters for alignment: MUST BE BASE 2')
    args.add_argument('-w', '--show-warning', type=bool, default= warning_Default, help='If anything other than "False", warnings will be shown')
    args.add_argument('-ts', '--send-text', type=bool, default= stageText_Default, help='Send a text after stage A and B are complete')
    args.add_argument('-tc', '--send-textemail', type=bool, default= completeText_Default, help='Send a text for completion, and an email with attached stats page on completion')
    args = args.parse_args()
    return args

def clusterSpecs(clustNum, readInClust, distortion, mergeCount, finalLen, startLen):
    clustInfo 		= []
    score 		= 0
    mino 		= 0
    minf 		= 0
    avgo 		= 0
    avgf 		= 0
    mediano 		= 0 
    medianf 		= 0
    maxo		= 0
    maxf 		= 0

    distortion 		= '%.10f' % float(distortion)
    score 		= mergeCount / float(float(distortion) * readInClust)
    score 		= '%.10f' % score
    mediano 		= numpy.median(startLen)
    medianf 		= numpy.median(finalLen)
    avgo 		= numpy.mean(startLen)
    avgo 		= '%.10f' % avgo
    avgf 		= numpy.mean(finalLen)
    avgf 		= '%.10f' % avgf
    mino 		= min(startLen)
    maxo 		= max(startLen)
    minf 		= min(finalLen)
    maxf 		= max(finalLen)


    clustInfo.append(clustNum)
    clustInfo.append(score)
    clustInfo.append(distortion)
    clustInfo.append(readInClust)
    clustInfo.append(mergeCount)
    clustInfo.append(mino)
    clustInfo.append(minf)
    clustInfo.append(avgo)
    clustInfo.append(avgf)
    clustInfo.append(mediano)
    clustInfo.append(medianf)
    clustInfo.append(maxo)
    clustInfo.append(maxf)



    return clustInfo




def main(args, t):
    gc.enable()
    # assign inputs to varibles 
    inFragReads 		= args.input_file
    outAlignedFragReads 	= args.output_file
    metric 			= args.cluster_metric
    refGenomeFile 		= args.reference_file
    window 			= args.ami_window
    kmerSize 			= args.kmer_size
    numOfClusters 		= args.cluster_count
    warning 			= args.show_warning
    recurse 			= args.recurse
    stageText 			= args.send_text
    completeText 		= args.send_textemail
    statsFile			= args.stats

    # TODO: Make these open to user
    AlignmentMethod		= AlignmentMethod_Default
    STOverlapThreshold 		= STOverlapThreshold_Default
    clusterDistMet		= clusterDistMet_Default
    userFreqContact		= userFreqContact_Default
    userEmail	 		= userEmail_Default

    # assign varibles for later use
    trainingSetSize 		= 0
    dimensions 			= 0
    timeStart			= 0
    timeStartTotal		= 0
    totalEnd			= 0
    timeStageA			= 0
    timeStageB			= 0
    timeStageC			= 0
    attchedFile 		= []
    recip 			= []



    # check that input are valid
    #   Input File:
    infileCheck = ''
    if (warning == True):
	if (inFragReads == inFragReads_Default):
	    infileCheck = (raw_input("WARNING:\n\tThe default input Fragread file is selected: '"+ str(inFragReads_Default) + "'\n    If this is ok, press enter.\n    If not, please enter the desired input file: \n"))
	if (infileCheck != ''):
	    inFragReads = infileCheck

    while (os.path.isfile(inFragReads) == False):
	print "ERROR:\n\tThe Fragment Reads file was invalid:\n\t    It does not appear to exist."
	inFragReads = (raw_input("Please enter the desired input file: \n"))
    
    #   Output File:
    if (warning == True):
	outfileCheck = ''
	while (os.path.isfile(outAlignedFragReads) == True):
	    print "WARNING:\n\tThe Output File Already Exists:\n\t    The contents of the file will be deleted unless you change the Output File"
	    outfileCheck = (raw_input("    If this is ok, press enter.\n    If not, please enter the desired output file:  \n"))
	    if (outfileCheck == ''):
		outfileCheck = outAlignedFragReads
		break
	outAlignedFragReads = outfileCheck

    #   Metric
    metricD = {"D": 0, "DICT": 0,  "DIC": 0,  "DICTIONARY": 0,  "A": 1,   "AMI": 1}
    metric = metric.upper()
    while (metric not in metricD) :
	
	metric = (raw_input("Please Enter the Desired Metric\n    D :  Dictionary\n    A :  AMI\n"))
	metric = metric.upper()

    #       Dictionary
    if (metricD[metric] == 0):
	if (warning == True):
	    if (refGenomeFile == refGenomeFile_Default):
		reffileCheck = ''
		refFileCheck = (raw_input("WARNING:\n\tThe default Reference Genome file is selected: '"+ str(refGenomeFile_Default) + "'\n    If this is ok, press enter.\n    If not, please enter the desired Reference Genome file: \n"))
		if (refFileCheck != ''):
		    refGenomeFile = refFileCheck

	while (os.path.isfile(refGenomeFile) == False):
	    print "ERROR:\n\tThe Reference Genome file file was invalid:\n\t    It does not appear to exist."
	    refGenomeFile = (raw_input("Please enter the desired input file: \n"))


    #       AMI
    elif (metricD[metric] == 1):
	if (warning == True):
	    if (window == window_Default):
		windowCheck = ''
		windowCheck = (raw_input("WARNING:\n\tThe default AMI Window Size is selected: '"+ str(window_Default) + "'\n    If this is ok, press enter.\n    If not, please enter the desired AMI Window Size: \n"))
		if (windowCheck != ''):
		    window = int(windowCheck)


    #       kmer ***NOT YET SUPPORTED***
    elif (metricD[metric] == 2):
	if (warning == True):
	    if (kmerSize == kmerSize_Default):
		kmerSizeCheck = ''
		kmerSizeCheck = (raw_input("WARNING:\n\tThe default Kmer Size is selected: '"+ str(kmerSize_Default) + "'\n    If this is ok, press enter.\n    If not, please enter the desired Kmer Size: \n"))
		if (kmerSizeCheck != ''):
		    kmerSize = int(kmerSizeCheck)


    #   Clusters
    if (warning == True):
	if (numOfClusters == numOfClusters_Default):
	    numOfClustersCheck = ''
	    numOfClustersCheck = (raw_input("WARNING:\n\tThe default Number of Clusters is selected: '"+ str(numOfClusters_Default) + "'\n    If this is ok, press enter.\n    If not, please enter the desired Number of Clusters. \n    **IT MUST BE BASE 2:** \n"))
	    if (numOfClustersCheck != ''):
		numOfClusters = int(numOfClustersCheck)

    clusterBaseCheck = math.log(numOfClusters, 2)
    while ( clusterBaseCheck != int(clusterBaseCheck) ):
	print "ERROR: \n\tNumber of Clusters was not Base 2"
	numOfClusters = int(raw_input("Please Enter the Desired Number of Clusters: \n"))
	clusterBaseCheck = math.log(numOfClusters, 2)
    '''
    # assign inputs to varibles 
    args.input_file		= inFragReads
    outAlignedFragReads 	= args.output_file
    metric 			= args.cluster_metric
    refGenomeFile 		= args.reference_file
    window 			= args.ami_window
    kmerSize 			= args.kmer_size
    numOfClusters 		= args.cluster_count
    warning 			= args.show_warning
    recurse 			= args.recurse
    '''




    stats = open(statsFile, 'w')
    stats.write(strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()) )
    stats.write('\nInput File:\t\t\t' + str(inFragReads) + '\nOutput File:\t\t\t' + str(outAlignedFragReads) )    
    if (metricD[metric] == 0):
	# Dictionary Metric
	stats.write("\nVectorization Method:\t\tDictionary\nReference Sequence File:\t" + str(refGenomeFile) )
    elif (metricD[metric] == 1):
	# Run AMI Metric
	stats.write("Vectorization Method:\tAverage Mutual Information\n\tSpecs:\tWindow Size:\t" + str(window) )
    elif (metricD[metric] == 2):
	# Run Kmer Metric
	print 'Waldo: \n"Phooey.\n"How on earth did you manage to get here?\n"' 
	stats.write("Vectorization Method:\tKmer\n\tSpecs:\tKmer Size:\t" + str(kmerSize) )
    stats.write("\nNumber of Clusters:\t\t" + str(numOfClusters) )
    if (clusterDistMet == 'E'):
	stats.write("\nCluster Assignment:\t\tEuclidean" )
    elif (clusterDistMet == 'C'):
	stats.write("\nCluster Assignment:\t\tCorrelational" )
    if (AlignmentMethod == 0):
	stats.write("\nAlignment Method:\t\tSuffix Tree based Overlap Search" )
	stats.write("\n\tOverlap Threshold:\t" + str(STOverlapThreshold) )
    stats.close()

    '''
 
    inFragReads 
    outAlignedFragReads 
    metric 
    refGenomeFile 
    window 
    kmerSize 
    numOfClusters
    warning 
    '''
    # msg = "TEST"
    # attchedFile = []
    # attchedFile.append(prioriFile)
    # recip = []
    # recip.append(userEmail)
    # sendMessage.message_Send_Full_Email(recip, "TEST ATTACHMENT", msg, attchedFile)

    
    #*********#
    # Stage A
    print "\t\t******** STARTING PROGRAM ********\n"
    if (numOfClusters > 1):
	timeStartTotal 	= timeit.default_timer()
	timeStart 		= timeStartTotal
	
	if (metricD[metric] == 0):
	    # Run Dictionary Metric
	    trainingSetSize, dimensions = A_fragmentCluster.main(refGenomeFile, inFragReads, A_vectorFile, A_fspecs)
	elif (metricD[metric] == 1):
	    # Run AMI Metric
	    trainingSetSize, dimensions = A_AMIGenVects.main(inFragReads, window, A_vectorFile, A_fspecs)
	elif (metricD[metric] == 2):
	    # Run Kmer Metric
	    print 'Waldo: \n"Phooey.\n"How on earth did you manage to get here?\n"' 

	totalEnd 		= timeit.default_timer()
	timeStageA		= totalEnd - timeStart
	timeStageA		= str(datetime.timedelta(seconds = timeStageA))
	msg 		= ("Stage A: building Vectors is complete.\nTime: " + str(timeStageA))
	try:
	    if (stageText == True):
		sendMessage.sms_message_Send(userFreqContact, msg)
	except:
	    print "Could not send text message"
	print msg
	stats = open(statsFile, 'a')
	stats.write("\nNumber of Reads:\t\t" + str(trainingSetSize) )
	stats.write("\nVector Dimensions:\t\t" + str(dimensions) )
	stats.write("\n\nStage A time:\t\t\t" + str(timeStageA) )
	stats.close()


	#*********#
	# Stage B
	timeStart 		= timeit.default_timer()

	#distort = B_trvqsplit.startClustering(A_vectorFile, clustCentFile, numOfClusters, dimensions, trainingSetSize)
	#B_clusterAssignment.mainClustAssignment(clustCentFile, inFragReads, A_vectorFile, clusterDistFolder, clusterFolder, clusterDistMet)
	#clustPrioriArray, stdDevArray = B_clusterPrioritization.clusterPriori(clusterDistFolder, clustCentFile, prioriFile)
	clustStats = B_ClusterAll.startClustering(A_vectorFile, clustCentFile, numOfClusters, dimensions, trainingSetSize, inFragReads, clusterFolder)
	# Clust Number, number of Contigs, Distortion Value

	totalEnd 		= timeit.default_timer()
	timeStageB		= totalEnd - timeStart
	timeStageB		= str(datetime.timedelta(seconds = timeStageB))
	msg 			= ("Stage B: Clustering Vectors is complete.\nTime: " + str(timeStageB))
	try:
	    if (stageText == True):
		sendMessage.sms_message_Send(userFreqContact, msg)
	except:
	    print "Could not send text message"
	print msg
	stats = open(statsFile, 'a')
	#stats.write("\nCluster Distortion:\t\t" + str(distort) )
	stats.write("\nStage B time:\t\t\t" + str(timeStageB) )
	stats.close()
    else :
	clustStats	= []
	temp 			= []
	temp.append(0)
	temp.append(2)
	temp.append(1)
	clustStats.append(temp)
	#clustPrioriArray[0][0] = 0
	#clustPrioriArray[0][1] = 2
	
    stats = open(statsFile, 'a')
    stats.write("\n\nCluster Specs\n")
    stats.write("01: Cluster Number\n02: KM score\n03: Distortion\n04: Number of Reads\n05: Number of Mergers\n06: Starting Smallest Read\n07: Final Smallest Read\n08: S. Average Size Read\n09: F. Average Size Read\n10: S. Median Size Read\n11: F. Median Size Read\n12: S. Largest Read\n13: F. Largest Size Read\n\n")
    stats.write("01\t02\t\t03\t\t04\t05\t06\t07\t08\t\t09\t\t10\t11\t12\t13\n\n")
    stats.close()

    #*********#
    # Stage C
    timeStart 		= timeit.default_timer()
    outfile 		= open(outAlignedFragReads, 'w')
    outfile.close()
    totalMerged 	= 0
    for i in range(len(clustStats)):
	stats = open(statsFile, 'a')
	if (clustStats[i][1] == 1):
	    # Move single seq to output file: no use doing anything else
	    filePath = str(clusterFolder) + str(clustStats[i][0])
	if (clustStats[i][1] > 1):
	    # Valid Cluster of more than 1 sequence
	    # Use Alignement Method Preselected
	    if (AlignmentMethod == 0):
		# Suffix Tree Direct Overlap
		# print "Valid: " + str(clustPrioriArray[i])
		try: 
		    if (numOfClusters > 1) :
			filePath = str(clusterFolder) + str(clustStats[i][0])
		    else:
			filePath = inFragReads
		    mergeCount, finalLen, startLen = alignMain.sfxTreeBasicAlign(filePath, STOverlapThreshold, outAlignedFragReads)
		    clustInfo = clusterSpecs(clustStats[i][0], clustStats[i][1], clustStats[i][2], mergeCount, finalLen, startLen)
		    totalMerged += mergeCount
		    for j in range(len(clustInfo)):
			stats.write(str(clustInfo[j]) + "\t")
		except:
		    print ("Error in cluster " + str(clustStats[i][0]))
		    stats.write("Error in cluster " + str(clustStats[i][0]))
	stats.write("\n")
	stats.close()


    totalEnd 		= timeit.default_timer()
    timeStageC		= totalEnd - timeStart
    timeStageC		= str(datetime.timedelta(seconds = timeStageC))
    msg 		= ("Stage C: Contig Alignment is complete.\nTime: " + str(timeStageC))
    print msg
    msg 		= ("Program Pass " + str(t + 1) + " complete.\nTotal Time:   " + str(datetime.timedelta(seconds = (totalEnd - timeStartTotal) )) + "\nTotal Merged: " + str(totalMerged))
    try:
	if (completeText == True):
	    sendMessage.sms_message_Send(userFreqContact, msg)
    except:
	print "Could not send text message"
    
    stats 		= open(statsFile, 'a')
    stats.write("\n\nTotal Reads Merged:\t\t" + str(totalMerged) )
    stats.write("\nStage C Time:\t\t" + str(timeStageC) )
    stats.write("\nTotal Run Time:\t\t" + str(totalEnd - timeStartTotal) )
    stats.close()

    msg 		= "Attached is the stats for the most recent alignment run.\n\tEnjoy!"
    attchedFile 	= []
    attchedFile.append(statsFile)
    recip 		= []
    recip.append(userEmail)
    #recip.append('X@gmail.com')
    try :
	if (completeText == True):
	    sendMessage.message_Send_Full_Email(recip, "Alignment Results", msg, attchedFile)
    except:
	print "Could not send stats page"
    t += 1
    if ( (recurse > 0) and (args.cluster_count > 1) ):
	args.cluster_count 	= args.cluster_count / 2
	recurse 		= recurse - 1
	args.input_file 	= args.output_file
	tempStr 		= args.output_file.split('.')[0]
	tempStr 		= tempStr.split('Pass')[0]
	args.output_file 	= str(tempStr + 'Pass' + str(t) + ".txt" )
	tempStr 		= args.stats.split('.')[0]
	tempStr			= tempStr.split('Run')[0]
	args.stats 		= str(tempStr + 'Run' + str(t) + ".txt" )
	args.show_warning 	= False
	main(args, t)
	
	

    gc.disable()
    return

t = 0
args = interface()
main(args, t)












