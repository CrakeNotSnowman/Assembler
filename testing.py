#!/usr/bin/env python
import sys
import os.path

import alignMain

inputFolder = "temp/classIIclusters/"
outputFolder = "temp/outputCIIclusters/"


def main():
    name_OF_Files = []
    for dirname, dirnames, filenames in os.walk(inputFolder):
	for filename in filenames:
	    name_OF_Files.append(os.path.join(dirname, filename))
    i = 0
    print len(name_OF_Files)
    for i in range(len(name_OF_Files)):
	print i
	filename = name_OF_Files[i]
	print filename
	outfileName = str(outputFolder) + str(i)
	alignMain.jbGenSuffTreeAlign(filename, outfileName)
'''
	try:
 	    alignMain.jbGenSuffTreeAlign(filename, outfileName)
	    print "HERROOOS"
	except:
	    print "Caught something... Maybe a cold"
'''

def bashCall(inarg):
    # error with file 11
    i = int(inarg)
    STOverlapThreshold = 50
    infile = "temp/classIclusters/" + str(i)
    outfile = "temp/classIclustersAssem/" + str(i)
    try:
	#alignMain.jbGenSuffTreeAlign(infile, outfile)
	alignMain.sfxTreeBasicAlign(infile, STOverlapThreshold, outfile)
    except:
	print "Error on file number ", i
    print "\t", i/float(67)*100, "% of groups have been assembled" 
	
    return


def difficultJob():
    # error with file 11
    startf = 41
    endf = 46
    for i in range(startf, endf):
	infile = "temp/classIIclusters/" + str(i)
	outfile = "temp/outputCIIclusters/" + str(i)
	try:
	    alignMain.jbGenSuffTreeAlign(infile, outfile)
	    #print "infile: ", 
	except:
	    print "Error on file number ", i
	print "\t", i/float(67)*100, "% of groups have been assembled" 
	
    return


#main()
#difficultJob()
inarg = sys.argv[1]
bashCall(inarg)
