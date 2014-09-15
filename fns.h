//
//  fns.h
//  CompressBases
//
//  Created by Austin Riffle on 11/28/13.
//  Copyright (c) 2013 Austin Riffle. All rights reserved.
//

#ifndef __CompressBases__fns__
#define __CompressBases__fns__

#define ENDL_HEX  0x05 //0000-0 101
#define COMP_BITS 3
#define CONT_FLAG 0x00
#define STOP_FLAG 0xFF

#define A         0x00 //0000-0 000
#define T         0x01 //0000-0 001
#define G         0x02 //0000-0 010
#define C         0x03 //0000-0 011
#define N         0x04 //0000-0 100
#define ENDL      0x05 //0000-0 101

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdint>
using namespace std;

struct MatchedIndices//alligned indexes for reference string to read
{
    int read_start;
    int read_stop;

    int ref_start;
    int ref_stop;

    bool operator < (const MatchedIndices &rhs) const { return read_start < rhs.read_start; }//overload for sorting
};

struct OpenSub {
    int start;
    int stop;
    bool isFull;
    OpenSub(int start, int stop, bool isFull) {
        this->start  = start;
        this->stop   = stop;
        this->isFull = isFull;
    }
};

void bitCompress(unsigned char value, int bits, int endFlag, ofstream& out_stream);
unsigned char getHex(const char letter);
vector<string> getVecReads(ifstream& in_stream, int n_reads);
MatchedIndices getLCSpair(string& read, string& ref, OpenSub os);
OpenSub getOpenStr(vector<MatchedIndices>& vecMi, string& read);

void setNodeDepths(shared_ptr<Node> root, int parent_depth);


//need to define these:
//void printAndCompressIndices(ofstream& out_stream, vector<MatchedIndices> vecMi);

#endif /* defined(__CompressBases__fns__) */
























