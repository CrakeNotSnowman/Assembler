#!/usr/bin/env python

def fna_read(name):
    fid = open(name)
    dna = fid.readline()
    dna = fid.read().replace("\n","")
    fid.close()
    return dna

def multifna_read(name):
    fragList = []
    fid = open(name)
    a = fid.readline()
    s = ''
    for line in fid:
        if line[0]== '>':
            fragList.append(s)
            s=''
        else:
            s = s+line.strip()
    fid.close()
    return fragList
