#!/usr/bin/env python

def removeACerrors(fileName, probE):
    oldFile = open(fileName, 'r')
    newFileName = fileName.split('.')
    newFileName[0] = str(newFileName[0]) + "Cleaned."
    newFileName = str(newFileName[0]) + str(newFileName[1])
    newFile = open(newFileName, 'w')
    #for line in oldFile:
    while True:
	header = oldFile.readline() #fasta format file
	if (header == ""):
	    break
	line = oldFile.readline()
	a = 0
	t = 0
	g = 0
	c = 0
	for i in range(len(line)):
	    if ((line[i] == "A") or (line[i] == "a")):
		a += 1
	    elif ((line[i] == "T") or (line[i] == "t")):
		t += 1
	    elif ((line[i] == "G") or (line[i] == "g")):
		g += 1
	    elif ((line[i] == "C") or (line[i] == "c")):
		c += 1
	if ((a+c)/(float(len(line))) < probE) :
	    newFile.write(header)
	    newFile.write(line)

    oldFile.close()
    newFile.close()


def main():
    fileName = "Outputs/OutputAlignments/UNMCPass03GreaterThan500baseContigs.iaf"
    removeACerrors(fileName, .95)

main()
