#!/usr/bin/env python
import sys
import fasta
import os.path

'''
Keith Murray
011414

Hypothesis: Using the language based LZ78 algorithm, fragment
    assembly is improved by achiving more accurate clusters.


Build Master Dictionary(refrence genome)
    word = Modified LZ78(refrence genome)
    if (word not in Master Dictionary):
	word count += 1
	add word to Master Dictionary, where word is a key
	add word count as value to key
	grow/append vector by one, value 1
	

    return (Master Dictionary, vector)

	* Structure of Master Dictionary/Vector Relationship	*
	* -word: as determined by LZ78				*
	* -word = key to Dictionary				*
	*     Dictionary[word]					*
	* -key to value in Dictionary				*
	*     Value = Dictionary[word]				*
	* -Value points to the position in the vector		*
	*     vector[Dictionary[word]]				*

For Each Fragment:
    initialize fragment vector to length of vector
    find word by LZ78
    if (word in Dictionary)
	fragment word count += 1
	fragment vector[Dictionary[word]] += 1
    normalize vector:
	fragment vector = fragment vector / fragment word count

    return (fragment vector)

'''

def build_Master_Dictionary_LZ78_av(s): # Modified LZ78 algorithm looking for words at each position
    s=s.upper()     # make the sequence uppercase
    D={'A':0,'C':1,'G':2,'T':3, 'N':4}     # initial dictionary, value point to vector location
    vector = [0, 0, 0, 0, 0]
    wordLoc = 0
    wordCount = 5
    refD = {0:'A', 1:'C', 2:'G', 3:'T', 4:'N'} 

    errorCount = 0
    for i in xrange(len(s)):        # do for the whole sequence
        word_len=1                  # start with 1 letter word
        while s[i:i+word_len] in D: #if we have this word in the dictionary            
            wordLoc = D[s[i:i+word_len]]  # reassigning 
	    vector[wordLoc] += 1  #incrementing vectors value
            word_len +=1            # we had the word in the dictionary, now look for a longer one with an extra letter at the end
            if (i+word_len)>=len(s):    #if we reached the end stop the process
                break
        if (i+word_len)<len(s):  # if we are at the end of the sequence, don't to anything
            D[s[i:i+word_len]] = wordCount   # this is a new word, add it to the dictionary and say it we encountered it once
	    refD[wordCount] = s[i:i+word_len]
	    wordCount += 1
	    vector.append(1)
	    if ( (len(D.keys()) != len(vector)) and (errorCount < 2) ):
		print '***********'
		print "ERROR OCCURED"
		print wordCount
		print D[s[i:i+word_len]]
		print s[i:i+word_len]
		print len(s)
		print i+word_len
		print '***********'
		errorCount += 1
    return D, vector, refD    # report the dictionary to the program


def build_Fragment_Dictionary_LZ78_av(s, D, vecLen): # Modified LZ78 algorithm looking for words at each position
    s=s.upper()     # make the sequence uppercase
    #D={'A':0,'C':1,'G':2,'T':3, 'N':4}     # initial dictionary, value point to vector location
    vector = [0] * vecLen
    wordLoc = 0
    wordCount = 5

    for i in xrange(len(s)):        # do for the whole sequence
        word_len=1                  # start with 1 letter word
        while s[i:i+word_len] in D: #if we have this word in the dictionary            
            wordLoc = D[s[i:i+word_len]]  # reassigning 
	    vector[wordLoc] += 1  #incrementing vectors value
            word_len +=1            # we had the word in the dictionary, now look for a longer one with an extra letter at the end
	    wordCount += 1
            if (i+word_len)==len(s):    #if we reached the end stop the process
                break
    
    return vector, wordCount    # report the fragVector to the program


def main(refGenomeFile, fragFile, outputFile, fspecs):
    #get user input:
    #	Reference Genome, Fragread file
    
    s = fasta.fna_read(refGenomeFile)
    masterD, masterVector, refD = build_Master_Dictionary_LZ78_av(s)
    trainingSetSize = 0
    with open(fragFile, 'r') as fragReads, open(outputFile, 'w') as output:
	#fragReads = open(fragFile, 'r')
	#output = open(outputFile, 'a')
	for line in fragReads:
	    fragVect, wordCount = build_Fragment_Dictionary_LZ78_av(line, masterD, len(masterVector))
	    i = 0
	    #tempD = refD
	    for i in range(len(fragVect)):
		count = fragVect[i]
		word = refD[i]
		temp02 = len(line) / float(len(word))
		temp = float(count/temp02) 
		output.write(str(float(temp)))
		if (i < len(fragVect)-1):
		    output.write("\t")	    
	    output.write("\n")
	    trainingSetSize += 1
    fragReads.close
    output.close    

    
    specs = open(fspecs, 'w')
    specs.write("Dimension:\t\t" + str(len(masterVector)) + "\n")
    specs.write("Training Set Size:\t" + str(trainingSetSize))
    specs.write("\n\nTraining Set File Path:\n    ")
    specs.write("/home/keith/Documents/filesForProgramming/Clustering/A_Generate_Vectors/Output_Vectors/Vector_File") 
    specs.close()
    #print "Training Set Size: ", trainingSetSize
    return trainingSetSize, len(masterVector)


























