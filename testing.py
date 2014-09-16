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


main()
