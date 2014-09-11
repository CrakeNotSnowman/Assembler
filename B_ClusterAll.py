#!/usr/bin/env python
'''
Keith Murray
in: (A_vectorFile, clustCentFile, numOfClusters, dimensions, trainingSetSize)(inFragReads, clusterFolder)

Out: stuff


-> run kamivq
   output file format:
     line 1: "Codebook centroid:"
     line 2: <Floats: centroid vector>
     line 3: "count[" <int: clusterNumb> "]=" <int: number of files in cluster>
     line 4: "Cluster membership:"
     line 5: <int values refering to the specific line in the A_vectorFile/inFragReads: those lines are members of that clust>
     line 6: "Cluster distortion: " <float: distortion value>
     line 7: "\n"

-> Clean B_Clusters files
	for i in range numOfClusters:
	    cluster = open(str(clusterFolder) + str(i), 'w')
	    cluster.close()

-> Parse file
-> Assign Clusters
-> Prioritize
  
'''


## import c++MAS
import kvqsplitForClust
import sys

##def assemble(filename, outfilename){
	#c++.MAS_main(filename, outfilename)
    #return

	

def startClusteringVTwo(A_vectorFile, clustCentFile, classIClustCount, numOfClusters, dimensions, trainingSetSize, pCIClustFile, clusterFolder, cIClusterCenters):
    # run kamivq
    kvqsplitForClust.cluster(A_vectorFile, clustCentFile, numOfClusters, dimensions, trainingSetSize)
    clustAssgn = [0]*trainingSetSize
    clusterStats = []
    distort = []

    # Clean B_Clusters Files
    for i in range (int(numOfClusters)):
	cluster = open(str(clusterFolder) + str(i + classIClustCount), 'w')
	cluster.close

    # Parse File
    clustSpecs = open(clustCentFile, 'r')
    while True:
	statsList = []
	line = clustSpecs.readline().strip()		# Line 1
	if (line == ''):
	    break
	line = clustSpecs.readline().strip()		# Line 2: Cluster Center: add to file, center line should corrosponde with file number
	clustCenters = open(cIClusterCenters, 'a')
	clustCenters.write(line)
	clustCenters.write("\n")
        clustCenters.close()
	line = clustSpecs.readline().strip()		# Line 3: Clust Number, number of contigs in Clust
	temp = line.split("]=")
	contigCount = int(temp[1])
	temp = temp[0].split("[")
	clustNum = int(temp[1])
	line = clustSpecs.readline().strip()		# Line 4
	line = clustSpecs.readline().strip()		# Line 5: All lines in A_VectorFile/inFragReads which belong to this cluster
	temp = line.split(' ')
	#print temp
	if (contigCount > 0) :
	    for i in range (len(temp)):
		clustAssgn[int(float(temp[i]))] = clustNum		#	Could also use an dictionary
	line = clustSpecs.readline().strip()		# Line 6: Distortion Value
	temp = line.split(":")
	temp = temp[1].strip()
	distortVal = float(temp)
	line = clustSpecs.readline().strip()		# Line 7

	distort.append(distortVal)
	statsList.append(clustNum)
	statsList.append(contigCount)
	statsList.append(distortVal)
	clusterStats.append(statsList)
    clustSpecs.close()

    # Assign Clusters
    contigs = open(pCIClustFile, 'r')
    i = 0
    while True:
	contig = contigs.readline().strip()
	clusterFile = open( str(clusterFolder) + str(clustAssgn[i] + classIClustCount), 'a')
	clusterFile.write(contig)
	clusterFile.write("\n")
	clusterFile.close()
	i += 1
	if ( i >= trainingSetSize):
	    break
    contigs.close()

    # Prioritize 
    distort, clusterStats = (list(t) for t in zip(*sorted(zip(distort, clusterStats))))
    #print clusterStats


    return clusterStats # Clust Number, number of Contigs, Distortion Value

def classIIClustering(A_vectorFile, clustCentFile, numOfClusters, dimensions, trainingSetSize, classIClusterFolder, classIIClusterFolder):
    print "\tA_vectorFile:\t\t", A_vectorFile
    print "\tclustCentFile:\t\t", clustCentFile
    print "\tnumOfClusters:\t\t", numOfClusters
    print "\tdimensions:\t\t", dimensions
    print "\ttrainingSetSize:\t", trainingSetSize
    print "\tclassIClusterFolder:\t", classIClusterFolder
    print "\tclassIIClusterFolder:\t", classIIClusterFolder

    # run kamivq
    kvqsplitForClust.cluster(A_vectorFile, clustCentFile, numOfClusters, dimensions, trainingSetSize)
    clustAssgn = [0]*trainingSetSize
    clusterStats = []
    distort = []

    # Clean B_Clusters Files
    #for i in range (int(numOfClusters)):
	#cluster = open(str(classIIClusterFolder) + str(i), 'w')
	#cluster.close

    # Parse File
    clustSpecs = open(clustCentFile, 'r')
    while True:
	statsList = []
	line = clustSpecs.readline().strip()		# Line 1
	if (line == ''):
	    break
	line = clustSpecs.readline().strip()		# Line 2: Cluster Center: Currently Do nothing with 
	line = clustSpecs.readline().strip()		# Line 3: Clust Number, number of contigs in Clust
	temp = line.split("]=")
	contigCount = int(temp[1])
	temp = temp[0].split("[")
	clustNum = int(temp[1])
	line = clustSpecs.readline().strip()		# Line 4
	line = clustSpecs.readline().strip()		# Line 5: All lines in A_VectorFile/inFragReads which belong to this cluster
	temp = line.split(' ')
	#print temp
	cII = open(str(classIIClusterFolder) + str(clustNum), 'a')
	for i in range(len(temp)):
	    if (temp[i] != ""):
		cI = open(str(classIClusterFolder) + str(temp[i]), 'r')
		while True:
		    line = cI.readline()
		    if (line == ""):
			break
		    cII.write(line)
	        cI.close()
	cII.close()
	line = clustSpecs.readline().strip()		# Line 6: Distortion Value
	temp = line.split(":")
	temp = temp[1].strip()
	distortVal = float(temp)
	line = clustSpecs.readline().strip()		# Line 7

	distort.append(distortVal)
	statsList.append(clustNum)
	statsList.append(contigCount)
	statsList.append(distortVal)
	clusterStats.append(statsList)
    clustSpecs.close()
    print line
   

    # Prioritize 
    distort, clusterStats = (list(t) for t in zip(*sorted(zip(distort, clusterStats))))
    #print clusterStats

    return clusterStats # Clust Number, number of Contigs, Distortion Value

def startClustering(A_vectorFile, clustCentFile, numOfClusters, dimensions, trainingSetSize, inFragReads, clusterFolder):
    # run kamivq
    kvqsplitForClust.cluster(A_vectorFile, clustCentFile, numOfClusters, dimensions, trainingSetSize)
    clustAssgn = [0]*trainingSetSize
    clusterStats = []
    distort = []

    # Clean B_Clusters Files
    for i in range (int(numOfClusters)):
	cluster = open(str(clusterFolder) + str(i), 'w')
	cluster.close

    # Parse File
    clustSpecs = open(clustCentFile, 'r')
    while True:
	statsList = []
	line = clustSpecs.readline().strip()		# Line 1
	if (line == ''):
	    break
	line = clustSpecs.readline().strip()		# Line 2: Cluster Center: Currently Do nothing with 
	line = clustSpecs.readline().strip()		# Line 3: Clust Number, number of contigs in Clust
	temp = line.split("]=")
	contigCount = int(temp[1])
	temp = temp[0].split("[")
	clustNum = int(temp[1])
	line = clustSpecs.readline().strip()		# Line 4
	line = clustSpecs.readline().strip()		# Line 5: All lines in A_VectorFile/inFragReads which belong to this cluster
	temp = line.split(' ')
	#print temp
	if (contigCount > 0) :
	    for i in range (len(temp)):
		clustAssgn[int(float(temp[i]))] = clustNum		#	Could also use an dictionary
	line = clustSpecs.readline().strip()		# Line 6: Distortion Value
	temp = line.split(":")
	temp = temp[1].strip()
	distortVal = float(temp)
	line = clustSpecs.readline().strip()		# Line 7

	distort.append(distortVal)
	statsList.append(clustNum)
	statsList.append(contigCount)
	statsList.append(distortVal)
	clusterStats.append(statsList)
    clustSpecs.close()

    # Assign Clusters
    contigs = open(inFragReads, 'r')
    i = 0
    while True:
	contig = contigs.readline().strip()
	clusterFile = open( str(clusterFolder) + str(clustAssgn[i]), 'a')
	clusterFile.write(contig)
	clusterFile.write("\n")
	clusterFile.close()
	i += 1
	if ( i >= trainingSetSize):
	    break
    contigs.close()

    # Prioritize 
    distort, clusterStats = (list(t) for t in zip(*sorted(zip(distort, clusterStats))))
    #print clusterStats


    return clusterStats # Clust Number, number of Contigs, Distortion Value
#startClustering('NULL', "testOutputCluster.txt", 128, 490, 10000, "/home/keith/Documents/filesForProgramming/Clustering/GoldenDataSet01/frags.txt", "B_Clusters/")

