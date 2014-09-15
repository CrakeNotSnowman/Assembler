#!/usr/bin/env python
'''
Keith Murray

Sequence Alignment Program:
  La Casa v1.00 March 06 2014

  Previous Versions:
    Venice Inn v1.00 Forked off March 06 2014
    Mr. C's v1.22 Forked off Feb 20 2014
  I need a better name for this




This program is intended to take fragment reads and assemble them 
over multiple itterations.

TODO:
* Prep Buffer for large files 
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

CODE FROM THE WEBS:
    FILE COMPRESSION
	import tarfile
	tar = tarfile.open("sample.tar.gz", "w:gz")
	for name in ["file1", "file2", "file3"]:
	    tar.add(name)
	tar.close()
    If you want to create a tar.bz2 compressed file, 
    just replace file extension name with ".tar.bz2" 
    and "w:gz" with "w:bz2".

    DELETING FILES
	import os

	filelist = [ f for f in os.listdir(".") if f.endswith(".bak") ]
	for f in filelist:
	    os.remove(f)

TOTEST:
* AMI
* non-default
    Outfile



'''
import os.path
import os
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
#import MAS			# and pray that it works 

# User Dependent Defaults
inFragReads_Default 		= 'KM_Hold/inputFile.txt'
outAlignedFragReads_Default 	= 'Outputs/OutputAlignments/outputPass0.txt'
metric_Default 			= 'N'
refGenomeFile_Default 		= 'KM_Hold/referenceFasta3.txt'
window_Default 			= 25
kmerSize_Default 		= 7
numOfClusters_Default 		= 128
clusterDistMet_Default		= 'E'
warning_Default 		= False
AlignmentMethod_Default 	= 0
STOverlapThreshold_Default 	= 30
userFreqContact_Default 	= 'X@vtext.com'
userEmail_Default 		= 'kX@gmail.com'
recurse_Default			= 0
stageText_Default 		= False
completeText_Default 		= False

# Internal Defaults 
clustCentFile 			= "temp/clustCents"
A_vectorFile 			= "A_Output_Vectors/Vector_File"
A_fspecs 			= "A_Output_Vectors/file_Specs"
prioriFile 			= "KM_Hold/priori"
clusterDistFolder 		= "B_ClusterDists/"
clusterFolder 			= "B_Clusters/"
statsFile 			= "Outputs/StatsFiles/StatsRun.txt"

# User values
toaddrsSMS  			= 'X@vtext.com' 
toaddrsMMS  			= 'X@vzwpix.com'


def getArgsfromFile():
    print "Eventually I'll finish this section"
    return

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


#def alignChunk(args, t) :


#def preClassIClusters(
'''
OUTLINE:
INFILE:
  Read in chunks of INFILE
    VECTORIZE
    CLUSTER -> Save cluster Centers in own file
	Cluster count == line of Cluster Center
    
'''

def master(args) :
     
    CRpPreCISizeLimit		= 10000		# Cluster Read per Class I size limit (Contig Count)
    DRpCISize			= 500		# Approx. Reads per Cluster in Class I
    DRpCIISize			= 80000		# Approx. Reads per Cluster in Class II: 
    MASfileSizeGoal		= 10		# MB
    MASAddedBufferLimit		= 2		# MB
    MASfileSizeTrimmed 		= 8		# MB
    bytesPerMB			= 1024*1024
    MASfileSizeGoal		= long(MASfileSizeGoal) * bytesPerMB		# Bytes
    MASAddedBufferLimit		= long(MASAddedBufferLimit) * bytesPerMB	# Bytes
    MASfileSizeTrimmed 		= long(MASfileSizeTrimmed) * bytesPerMB		# Bytes
    

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
    tsSize			= []
    metricD 			= {"D": 0, "DICT": 0,  "DIC": 0,  "DICTIONARY": 0,  "A": 1,   "AMI": 1}


    classIClustCount 		= 0 			# Class 1 clusters: ~78 contigs per file
    totalContigCount 		= 0
    preCIClustFolder		= "temp/preCIfiles/"	# ~10,000 contigs per file
    preCIClustCount		= 0
    preCIVectorsFolder		= "temp/preCIVectors/"
    classICenters		= "temp/classICenters.txt"
    classIClustersFolder	= "temp/classIclusters/"
    classIIClusterFolder	= "temp/classIIclusters/"
    clustCentFile		= "temp/clusterCenters.txt"
    cIClusterCenters		= "temp/cIClustCenters.txt"
    cIIClusterCenters		= "temp/cIIClustCenters.txt"
    



    

    cleanUpTemp("temp/")

    # Generate 10,000 reads in pre class 1 files
    inFile = open(inFragReads, 'r')
    while True:
	contigCount = 0					# Setting the counter for contigs/file
	preCIClustFile = open(str(preCIClustFolder) + str(preCIClustCount), 'w' )
	while contigCount < CRpPreCISizeLimit:
	    contig = inFile.readline()
	    if (contig == ''):
		break
	    preCIClustFile.write(contig)
	    contigCount += 1	
	preCIClustFile.close()
	preCIClustCount += 1
	totalContigCount += contigCount
	if (contig == ''):
	    break
    inFile.close()

    # Generate PreCI Vectors
    #cent = open(classICenters, 'w')
    #cent.close()
    clustCenters = open(cIClusterCenters, 'w')
    clustCenters.close()
    for i in range(preCIClustCount):
	print i/float(preCIClustCount)
	pCIClustFile = str(preCIClustFolder) + str(i)	
	A_vectorFile = str(preCIVectorsFolder) + str(i)
	# Run Dictionary Metric
	trainingSetSize, dimensions = A_fragmentCluster.main(refGenomeFile, pCIClustFile, A_vectorFile, A_fspecs)
	tsSize.append(trainingSetSize)	
	# Find out size and the such
        tempVar = int(math.log(trainingSetSize/float(DRpCISize), 2))
	
	numOfClusters = int(math.pow(2, tempVar))
	clustStats = B_ClusterAll.startClusteringVTwo(A_vectorFile, clustCentFile, classIClustCount, numOfClusters, dimensions, trainingSetSize, pCIClustFile, classIClustersFolder, cIClusterCenters)
	classIClustCount += numOfClusters



    # Cluster/Align Cluster chunks
    tempVar = int(math.log(classIClustCount * DRpCISize/float(DRpCIISize), 2))
    numOfCIIClusters = int(math.pow(2, tempVar))
    #numOfCIIClusters = 2
    # Group/Merge Clusters
    clustStats = B_ClusterAll.classIIClustering(cIClusterCenters, cIIClusterCenters, numOfCIIClusters, dimensions, classIClustCount, classIClustersFolder, classIIClusterFolder)

    # Refine CII Clusters: Ensure filesize does not exceed (    MASfileSizeGoal + MASAddedBufferLimit )
    # MASfileSizeTrimmed
    fullCIIClusterCount = numOfCIIClusters
    i = 0
    MASsizeLimit = MASfileSizeGoal + MASAddedBufferLimit

    for i in range(numOfCIIClusters):
	# Looking at a specific file
	currFile = str(classIIClusterFolder) + str(i)
	statinfo = os.stat(currFile)
	fSize = statinfo.st_size		# Bytes
	byteCounter = [0, 0, 0]			# MB, kB, B
	bytesTotal = long(0)
	if (fSize > MASsizeLimit):
	    # Making sure clusters are smaller than the limit
	    remainA = fSize % MASfileSizeGoal
	    remainC = MASfileSizeGoal - remainA
	    remainB = fSize % MASfileSizeTrimmed
	    remainD = MASfileSizeTrimmed - remainD

	    minValue = min(remainA, remainB, remainC, remainD)


	    if (remainA == minValue):
		# best files will be files close to 10 MB limit
		#    And the remainder will be distributed in each file
		remainA = fSize / float(MASfileSizeGoal)
		# remainA is now the limit, in MB for each file. 
		if (int(remainA) > 1) :
		    sourceFile = open(currFile, 'r')
		    #for j in range(0,
		    line = sourceFile.readline()
		    lineSize = len(line)
		    
		    '''
		    ok so I need to see if I need to make more than one file
		    Here I do, so next I need to write the new file,
		    incriment the counter
		    open the new file, and source file,
		    copy line by line 

		    '''
		    
	    elif (remainB == minValue):
		print "HI"
	    elif (remainC == minValue):
		print "HI"
	    elif (remainD == minValue):
		print "NI"

    # ReCluster/Align Groups

	
    '''
    Read INFILE into 10,000 line chunks
    Save each chunk to file N
    For file i in N
       cluster (64)
    	    For j in 64:
    		Save Center
    		Align cluster: save alignment in cluster
     Cluster the cluster centers 

    Number of Lines 			= 5,466,888
    Number of chunks 			= 547
    Start subCluster (Class 1) Count	= 64
    number of total clusters 		= Start subCluster Count * Number of chunks
					= 547 * 64
					= 35,008
    number of Class 2 Clusters 		= number of total clusters / subCluster (Class 1) Count
					= 35,008 / 64
					= 547
					= 512 (base 2)
    
    ((Number of lines/ 10,000) + 1 ) / subClusterCount)

    '''
    return

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
    #recip.append('jbohac61@gmail.com')
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

def prepPreCIfiles(inFragReads, startingIndex, CRpPreCISizeLimit, preCIClustFolder):
    totalContigCount = 0
    preCIClustCount = startingIndex
    inFile = open(inFragReads, 'r')
    while True:
	contigCount = 0
	preCIClustFile = open(str(preCIClustFolder) + str(preCIClustCount), 'w' )
	while contigCount < CRpPreCISizeLimit:
	    contig = inFile.readline()
	    if (contig == ''):
		break
	    preCIClustFile.write(contig)
	    contigCount += 1	
	preCIClustFile.close()
	preCIClustCount += 1
	totalContigCount += contigCount
	if (contig == ''):
	    break
    inFile.close()
    return preCIClustCount, totalContigCount

def generateVectsAndClusts(preCIClustFileCount, cIClusterCenters, refGenomeFile, preCIClustFolder, preCIVectorsFolder, DRpCISize, classIClustersFolder):
    tsSize = []
    clustCenters = open(cIClusterCenters, 'w')
    clustCenters.close()
    classIClustCount = 0
    for i in range(preCIClustFileCount):
	print i/float(preCIClustFileCount)
	pCIClustFile = str(preCIClustFolder) + str(i)	
	A_vectorFile = str(preCIVectorsFolder) + str(i)
	# Run Dictionary Metric
	trainingSetSize, dimensions = A_fragmentCluster.main(refGenomeFile, pCIClustFile, A_vectorFile, A_fspecs)
	tsSize.append(trainingSetSize)	
	# Find out size and the such
        tempVar = int(math.log(trainingSetSize/float(DRpCISize), 2))
	print tempVar
	numOfClusters = int(math.pow(2, tempVar))
	print A_vectorFile, clustCentFile, classIClustCount, numOfClusters, dimensions, trainingSetSize, pCIClustFile, classIClustersFolder, cIClusterCenters
	clustStats = B_ClusterAll.startClusteringVTwo(A_vectorFile, clustCentFile, classIClustCount, numOfClusters, dimensions, trainingSetSize, pCIClustFile, classIClustersFolder, cIClusterCenters)
	print "BYE"
	classIClustCount += numOfClusters

    return classIClustCount, dimensions
    

def groupClusters(classIClustCount, DRpCISize, DRpCIISize, cIClusterCenters, cIIClusterCenters, dimensions, classIClustersFolder, classIIClusterFolder, MASfileSizeGoal, MASAddedBufferLimit):
    tempVar = int(math.log(classIClustCount * DRpCISize/float(DRpCIISize), 2))
     
    print "Log[ " + str(classIClustCount) + " * " + str(DRpCISize) + " / " + str(DRpCIISize) + " ]"
    print tempVar
    numOfCIIClusters = int(math.pow(2, tempVar))
    numOfCIIClusters = 4
    # ** I'm being lazy here: go back and thance this later **
    if (numOfCIIClusters < 2):
	numOfCIIClusters = 1
    origNumCIIClusts = numOfCIIClusters
    print "APES"
    print "(cIClusterCenters, cIIClusterCenters, numOfCIIClusters, dimensions, classIClustCount, classIClustersFolder, classIIClusterFolder)"
    print (cIClusterCenters, cIIClusterCenters, numOfCIIClusters, dimensions, classIClustCount, classIClustersFolder, classIIClusterFolder)
    clustStats = B_ClusterAll.classIIClustering(cIClusterCenters, cIIClusterCenters, numOfCIIClusters, dimensions, classIClustCount, classIClustersFolder, classIIClusterFolder)
    
    # TESTING TIME OHH YEAH
    MASfileSizeGoal = 1024*10
    MASAddedBufferLimit = 1024*2
    
    MASsizeLimit = MASfileSizeGoal + MASAddedBufferLimit
    print numOfCIIClusters

    for i in range(origNumCIIClusts):
	# Looking at a specific file
	currFile = str(classIIClusterFolder) + str(i)
	statinfo = os.stat(currFile)
	fSize = statinfo.st_size		# Bytes
	print fSize
	fSize = os.path.getsize(currFile)
	print fSize
	byteCounter = [0, 0, 0]			# MB, kB, B
	bytesTotal = long(0)
	bytesRunningTotal = long(0)
	print "I EXIST YOU DAMN EXISTENTIALISTS" + "!"*i + " Oh, and this file is " + str(fSize / float(1024*1024) ) + " Mb"
	if (fSize > MASsizeLimit):
	    # Check the remainder of the cut size,
	    cuts = int(fSize / MASfileSizeGoal) - 1
	    remainder = fSize % (MASfileSizeGoal)
	    print "\t---Started with ", numOfCIIClusters, " Clusters---"
	    print "\tfSize:\t", fSize
	    print "\tCuts:\t", cuts
	    print "\tRemain:\t", remainder
	    print "\ti:\t", i
	    # If it's > ~2 MB, make it into it's own file,
	    if (remainder > MASAddedBufferLimit):
		cuts += 1
	    print "\tCuts:\t", cuts
	    # ----------- White Board Code (WBC) ----------- #
	    tofile = open(str(classIIClusterFolder) + "temp", "w")
	    fromfile = open(currFile, "r")
	    newFileCount = 0
	    while True:
		line = fromfile.readline()
		if (line == ""):
		    print "\t\tFinalSize:     ", bytesTotal
		    print "\t\tRemainderSize: ", fSize - bytesRunningTotal
		    tofile.close()
		    break
		bytesTotal = bytesTotal + len(line)#sys.getsizeof(line)
		#if (newFileCount == cuts):
		#    #print "\t\t\t\t", bytesTotal
		if ((bytesTotal > MASfileSizeGoal) and (newFileCount < cuts)):
		    #print "\t\t\tBytesCount", bytesTotal, " / ", MASfileSizeGoal
		    #print "\t\t\tBytesReal", sys.getsizeof(line)
		    tofile.close()
		    tofile = open(str(classIIClusterFolder) + str(numOfCIIClusters), "w")
		    numOfCIIClusters += 1
		    newFileCount += 1
		    bytesRunningTotal += bytesTotal
		    bytesTotal = len(line)#sys.getsizeof(line)
		    #print bytesTotal, len(line), len(line)*4
		tofile.write(line)
	    try:
		tofile.close()
	    except:
		print "\t\t\tFile already closed"
	    print 
	    fromfile.close()
	    # Time to move the temp file back to its original file
	    tofile = open(currFile, "w")
	    fromfile = open(str(classIIClusterFolder) + "temp", "r")
	    while True:
		line = fromfile.readline()
		if (line == ""):
		    break
		tofile.write(line)
	    tofile.close()
	    fromfile.close()
	    # Now remove temp
	    try:
		os.remove(str(classIIClusterFolder) + "temp")
	    except OSError:
		print "Could not delete " + str(filename) + " during cleanup"
		pass
	    print "\t----Ended with ", numOfCIIClusters, " Clusters----"
	   
	    print "APPLES"
	    

    return numOfCIIClusters
    

def assembleCIIClusts(numOfCIIClusters, infolder, outfolder):
    for i in range(numOfCIIClusters):
	alignMain.jbGenSuffTreeAlign(str(infolder) + str(i), str(outfolder) + str(i))

    return
    


def cleanUpTemp(start_Folder):
    print "Cleaning up " + str(start_Folder)
    name_OF_Files = []
    for dirname, dirnames, filenames in os.walk(start_Folder):
	for filename in filenames:
	    name_OF_Files.append(os.path.join(dirname, filename))
    i = 0
    for i in range(len(name_OF_Files)):
	filename = name_OF_Files[i]
	try:
	    os.remove(filename)
	except OSError:
	    print "Could not delete " + str(filename) + " during cleanup"
	    pass
    return

def assembleHome(args):
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


    CRpPreCISizeLimit		= 10000		# Cluster Read per Class I size limit (Contig Count)
    DRpCISize			= 500		# Approx. Reads per Cluster in Class I
    DRpCIISize			= 80000		# Approx. Reads per Cluster in Class II: 
    MASfileSizeGoal		= 10		# MB
    MASAddedBufferLimit		= 2		# MB
    MASfileSizeTrimmed 		= 8		# MB
    bytesPerMB			= 1024*1024
    MASfileSizeGoal		= long(MASfileSizeGoal) * bytesPerMB		# Bytes
    MASAddedBufferLimit		= long(MASAddedBufferLimit) * bytesPerMB	# Bytes
    MASfileSizeTrimmed 		= long(MASfileSizeTrimmed) * bytesPerMB		# Bytes
    classIClustCount 		= 0 			# Class 1 clusters: ~78 contigs per file
    totalContigCount 		= 0
    preCIClustFolder		= "temp/preCIfiles/"	# ~10,000 contigs per file
    preCIClustCount		= 0
    preCIVectorsFolder		= "temp/preCIVectors/"
    classICenters		= "temp/classICenters.txt"
    classIClustersFolder	= "temp/classIclusters/"
    classIIClusterFolder	= "temp/classIIclusters/"
    cIIClusterOutputFolder	= "temp/outputCIIclusters/"
    clustCentFile		= "temp/clusterCenters.txt"
    cIClusterCenters		= "temp/cIClustCenters.txt"
    cIIClusterCenters		= "temp/cIIClustCenters.txt"
    cIIClusterOutputFolder	= "temp/outputCIIclusters/"
    

    '''
    Go through all of the temp folders:
	Delete everything 

    1.  On startup: Generate files of 10,000 contigs in 
	temp/preCIfiles files will be labeled as ints 0-N
    2.  Generate files of vectors for each file 0-N in a seperate
	folder, temp/preCIVectors/
    3.	For each contig&vector file pair, break it up into k clusters
	and save each cluster center in a center file, cIClusterCenters.txt
	there will be P centers, where P = N * k
    4.	Cluster the centers into Q clusters, and generate the coorosponding
	files in temp/classIIclusters/
    5.  If file i in temp/classIIclusters/ is >10Mb, split the file, name
	one file i, and the other Q. set Q = Q + 1.
    6.  Assemble the classIIclusters, output the results into preCIfiles
	under the same name, then merge files into one flatfile in Output
    7.  if the final files are too large, ensure the files are about 10,000
	contigs in length, and repeat from step 2
    '''
    # Delete everything in temp
    cleanUpTemp("temp/")



    #1.	On startup: Generate files of 10,000 contigs in 
    #	temp/preCIfiles files will be labeled as ints 0-N
    # IN:	Input file, starting index, contig count limit, folderLocation
    # OUT:	total number of files (N)
    startingIndex = 0
    preCIClustFileCount, totalContigCount = prepPreCIfiles(inFragReads, startingIndex, CRpPreCISizeLimit, preCIClustFolder)

    
    # 2. Generate files of vectors for each file 0-N in a seperate
    #    folder, temp/preCIVectors/
    # 3. For each contig&vector file pair, break it up into k clusters
    #    and save each cluster center in a center file, cIClusterCenters.txt
    #    there will be P centers, where P = N * k
    # IN: preCIClustFileCount, cIClusterCenters file, refGenomeFile
    classIClustCount, dimensions = generateVectsAndClusts(preCIClustFileCount, cIClusterCenters, refGenomeFile, preCIClustFolder, preCIVectorsFolder, DRpCISize, classIClustersFolder)


    # 4. Cluster the centers into Q clusters, and generate the coorosponding
    #    files in temp/classIIclusters/
    groupClusters(classIClustCount, DRpCISize, DRpCIISize, cIClusterCenters, cIIClusterCenters, dimensions, classIClustersFolder, classIIClusterFolder, MASfileSizeGoal, MASAddedBufferLimit)

    assembleCIIClusts(numOfCIIClusters)
    # Delete everything in temp
    #cleanUpTemp()
    gc.disable()
    return

t = 0
getArgsfromFile()
args = interface()
#main(args, t)
#master(args)
assembleHome(args)











