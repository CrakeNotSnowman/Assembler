#!/usr/bin/env python
'''
This generates the AMI vector for a specific file passed into it,
and returns the vector

window size
length
vector


'''
import os.path
from Bio import SeqIO
import numpy
import math
import sys


def dnaDataConvert(name):
#Converts string of letters (ATCG) into an
#array of numbers and returns it as "dataConverted"
        # Base A becomes 0
        # Base T becomes 1
        # Base G becomes 2
        # Base C becomes 3
    
    name = name.replace('A', '0')
    name = name.replace('a', '0')
    name = name.replace('T', '1')
    name = name.replace('t', '1')
    name = name.replace('G', '2')
    name = name.replace('g', '2')
    name = name.replace('C', '3')
    name = name.replace('c', '3')
    name = name.replace('N', '4')
    name = name.replace('R', '4')
    name = name.replace('K', '4')    
    name = name.replace('Y', '4')
    name = name.replace('W', '4')
    name = name.replace('M', '4')
    name = name.replace('S', '4')    
    name = name.replace('B', '4')
    name = name.replace('V', '4')
    name = name.replace('D', '4')
    name = name.replace('E', '4')
    name = name.replace('F', '4')
    name = name.replace('H', '4')
    name = name.replace('I', '4')    
    name = name.replace('J', '4')
    name = name.replace('L', '4')
    name = name.replace('O', '4')
    name = name.replace('P', '4')  
    name = name.replace('Q', '4')
    name = name.replace('U', '4')
    name = name.replace('X', '4')
    name = name.replace('Z', '4')
    name = name.replace('b', '4')
    name = name.replace('d', '4')
    name = name.replace('e', '4')
    name = name.replace('f', '4')
    name = name.replace('h', '4')
    name = name.replace('i', '4')
    name = name.replace('j', '4')
    name = name.replace('k', '4')
    name = name.replace('l', '4')
    name = name.replace('m', '4')
    name = name.replace('n', '4')
    name = name.replace('o', '4')
    name = name.replace('p', '4')
    name = name.replace('q', '4')
    name = name.replace('r', '4')
    name = name.replace('s', '4')
    name = name.replace('u', '4')
    name = name.replace('v', '4')
    name = name.replace('w', '4')
    name = name.replace('x', '4')
    name = name.replace('y', '4')
    name = name.replace('z', '4')
    
    baseList = list(name)
    baseList = [ int(a) for a in baseList ] 
    length = len(baseList)                    
                                                      

    dataConverted = baseList
    return dataConverted



def amiMatrixProbability(dataConverted) :
    fullBaseList = dataConverted
    length = len(fullBaseList)
    indexLoopCount = 0
    independentProb = [[0]*5]


    while indexLoopCount < length :
        if fullBaseList[indexLoopCount] == 0 :
            independentProb[0][0] += 1
        elif fullBaseList[indexLoopCount] == 1 :
            independentProb[0][1] += 1
        elif fullBaseList[indexLoopCount] == 2 :
            independentProb[0][2] += 1
        elif fullBaseList[indexLoopCount] == 3 :
            independentProb[0][3] += 1
            
        indexLoopCount +=1
    
    i = 0
    a = 0
    for i in range (5) :
        a = independentProb[0][i]
        if ((length != 0) and (a != 4)):
            a = float(a) / length
        independentProb[0][i] = a
        a = 0
        i += 1

           

    return independentProb
    


def amiJointProbability(dataConverted, window) :
    fullBaseList = dataConverted
    length = len(fullBaseList)
    independentProb = amiMatrixProbability(dataConverted)
    #window
    k = 0
    amiProfileList = []
    #window = 64

    for k in range(1, int(window) + 2) :
        i = 0
        j = 0
        jointProb =  [ [ 0. for a in range(5) ] for b in range(5) ]
        
        for window in range(length - k) :
            i = fullBaseList[window]
            j = fullBaseList[window + k]
            jointProb[i][j] += 1

        i = 0
        j = 0
        for i in range(4) :
            for j in range(4) :
                if ((length-k) != 0) :
                    jointProb[i][j] = (jointProb[i][j]) / (length - k) 

        

        
        #Log Math time

        iOfK = 0
        logNum = 0
        for i in range(4) :
            for j in range(4) :
                if ((independentProb[0][i] != 0)and(independentProb[0][j] !=0)) : 
                    logNum = (float(jointProb[i][j]) / (float(independentProb[0][i]) * independentProb[0][j]) )                
                if (logNum != 0) :
                    iOfK += (float(jointProb[i][j]) * math.log(logNum))
                
        amiProfileList.append(iOfK)        
        #print(iOfK)

    #*****#End for loop
    
    
    
    return amiProfileList


def amiComputeMain(seq, window) :
    seq = seq.strip()
    dnaConvertedData = dnaDataConvert( seq)
    amiProfileList = amiJointProbability(dnaConvertedData, window)
    return amiProfileList



def main(fragFile, window, outputFile, fspecs):
    #get user input:
    #	Reference Genome, Fragread file
    

    trainingSetSize = 0
    with open(fragFile, 'r') as fragReads, open(outputFile, 'w') as output:
	#fragReads = open(fragFile, 'r')
	#output = open(outputFile, 'a')
	for line in fragReads:
	    ami = amiComputeMain(line, window)
	    i = 0
	    for i in range(len(ami)):
		temp = ami[i]
		output.write(str(float(temp)))
		if (i < len(ami)-1):
		    output.write("\t")	    
	    output.write("\n")
	    trainingSetSize += 1
    fragReads.close
    output.close    


    specs = open(fspecs, 'w')
    specs.write("Dimension:\t\t" + str(len(ami)) + "\n")
    specs.write("Training Set Size:\t" + str(trainingSetSize))
    specs.write("\n\nTraining Set File Path:\n    ")
    specs.write("/home/keith/Documents/filesForProgramming/Clustering/A_Generate_Vectors/Output_Vectors/Vector_File") 
    specs.close()
    print "Training Set Size: ", trainingSetSize
    return trainingSetSize, len(ami)









