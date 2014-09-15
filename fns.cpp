//
//  fns.cpp
//  CompressBases
//
//  Created by Austin Riffle on 11/28/13.
//  Copyright (c) 2013 Austin Riffle. All rights reserved.
//

#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdint>
#include "GeneralizedSuffixTree.h"
#include "fns.h"

using namespace std;

void setNodeDepths(shared_ptr<Node> root, int parent_depth){
    shared_ptr<EdgeBag> edges = root->getEdges();
    vector<shared_ptr<Edge> > e = edges->getValues();

    for(unsigned int i=0; i<e.size(); i++){
        e[i]->getDest()->parent = root;
        e[i]->getDest()->nodeDepth = e[i]->getLabel().length() +parent_depth;
        //
        setNodeDepths(e[i]->getDest(), e[i]->getLabel().length()+parent_depth);
    }

}



vector<string> getVecReads(ifstream& in_stream,int n_reads) {
    int i=0;
    string tmp  = "";
    string tmp2 = "";
    string read = "";
    vector<string> vec_reads = vector<string>(0);
    vec_reads.push_back(tmp);

    while (i<n_reads && getline(in_stream, read)) {


        if (read.empty())
            break;
        else {
            //read = read.substr(0, read.size()-1);
            vec_reads.push_back(read);
            read.clear();
        }
    i++;
    }

    return vec_reads;
}


void bitCompress(unsigned char value, int bits, int endFlag, ofstream& out_stream)//fix output to ofstream -- fix output ?encoding? -> output data as single char in .txt file // from sayood
{
    static uint8_t data  = 0;
    static int bits_used = 0;

    int byte_size = 8;//8bits = 1 byte
    int shift     = byte_size - bits_used - bits;

    if(shift >= 0) {
        data         |= (value << shift);
        bits_used    += bits;
        if(bits_used == byte_size) {
            out_stream.put(data);
            data      = 0;
            bits_used = 0;
        }
    }

    else {
        data |= (value >> -shift);
        out_stream.put(data);
        data  = 0;
        shift = byte_size + shift;

        if(shift >= 0) {
            data     |= (value << shift);
            bits_used = byte_size - shift;
        } else {
            data     |= (value >> -shift);
            out_stream.put(data);
            data      = 0;
            shift     = byte_size + shift;
            data     |= (value << shift);
            bits_used = byte_size - shift;
        }
    }

    if(endFlag && bits_used != 0)
        out_stream.put(data);
}

unsigned char getHex(const char letter) {//changes base to encoded value
    if (letter == 'A')
        return A;
    else if (letter == 'T')
        return T;
    else if (letter == 'G')
        return G;
    else if (letter == 'C')
        return C;
    else
        return N;
}

MatchedIndices getLCSpair(string& read, string& ref, OpenSub os) {
	MatchedIndices indices;
    string tmp = read.substr(os.start, (os.stop - os.start + 1));
	GeneralizedSuffixTree gst(tmp, ref);

	string substr = gst.getLCS_k_strs(2);//*** This is returning the incorrect substring ***

	indices.read_start = (int)tmp.find(substr) + os.start;//double check to make sure correct --- may be able to just use os.start
	indices.read_stop = indices.read_start + (int)substr.size() - 1;
	indices.ref_start = (int)ref.find(substr);
	indices.ref_stop = indices.ref_start + (int)substr.size() - 1;

	return indices;
}

OpenSub getOpenStr(vector<MatchedIndices>& vecMi, string& read) {
    sort(vecMi.begin(), vecMi.end());

    if (vecMi.front().read_start > 0) {
        OpenSub os(0, vecMi.front().read_start - 1, false);
        return os;
    } else if (vecMi.back().read_stop != (int)(read.length()-1)) {
        OpenSub os(vecMi.back().read_stop + 1, (int)read.length() - 1, false);
        return os;
    } else
        for (unsigned int i=0; i<vecMi.size() - 1; i++)
            if ((vecMi[i].read_stop + 1) != (vecMi[i+1].read_start) ) {
                OpenSub os(vecMi[i].read_stop + 1, vecMi[i+1].read_start - 1, false);
                return os;
            }

    OpenSub os(0,0,true);
    return os;
}

//bool getOpenSubStr(string& read, vector<line_ref>& line_read_i, sub_stuct& read_substr) {//can probably be made more effeciently
//    vector<int> filler (read.size(), 0);//initializes vector with 0's
//
//    if (line_read_i.empty()) {//checks for references, if none sets subStr equal whole line
//        read_substr.substr = read;
//        read_substr.start  = 0;
//        return true;
//    }
//
//    for (int i=0; i<line_read_i.size(); i++)//places ones where substrings have been taken
//        for (int j=line_read_i[i].line_i[0]; j<=line_read_i[i].line_i[1]; j++)
//            filler[j]=1;
//
//    int start = -1, stop = 0;
//    for (int k=0; k<filler.size(); k++) {//gets start and stop point of first available substring
//        if (filler[k]==0 && start == -1)
//            start = k;
//        if (start != -1 && filler[k]==1) {
//            stop  = k-1;
//            break;
//        }
//    }
//
//    if (start == -1)//returns false if no substring available
//        return false;
//
//    int length=(stop-start)+1;//generates substring
//    read_substr.substr = read.substr(start, length);
//    read_substr.start  = start;
//
//    return true;
//}


//void getLCS(const string& ref_str, sub_stuct& read_substr, vector<MatchedIndices>& line_read_i) {// adapted code -- based off dynamic programming -- from wikipedia -- NOT WORKING LIKE EXPECTED!!!
//    int substr_size  = (int)read_substr.substr.size();
//    int *curr = new int[substr_size];
//    int *prev = new int[substr_size];
//    int *swap = nullptr;
//    int ref_end=0, lcs_end=0, max_substr_length=0;
//
//    for(int i=0; i<ref_str.size(); i++) {
//        for(int j=0; j<substr_size; j++) {
//            if(ref_str[i] != read_substr.substr[j])
//                curr[j] = 0;
//            else {
//                ref_end=i;
//                lcs_end=j;
//                if(i == 0 || j == 0)
//                    curr[j] = 1;
//                else
//                    curr[j] = 1 + prev[j-1];
//                max_substr_length = max(max_substr_length, curr[j]);
//            }
//        }
//        swap=curr;
//        curr=prev;
//        prev=swap;
//    }
//
//    delete [] curr;
//    delete [] prev;
//
//    int ref_beg = ref_end - (max_substr_length-1);
//    int lcs_beg = lcs_end - (max_substr_length-1);
//
//    int offset = read_substr.start;
//    insertData(ref_beg, ref_end, lcs_beg+offset, lcs_end+offset, line_read_i);
//}

//void insertData(const int ref_beg, const int ref_end, const int substr_beg, const int substr_end, vector<MatchedIndices>& line_read_i) {
//    line_read_i.push_back(line_ref() );
//    int size = (int)line_read_i.size() - 1;
//
//    line_read_i[size].read.push_back(substr_beg);
//    line_read_i[size].read.push_back(substr_end);
//    line_read_i[size].ref.push_back(ref_beg);
//    line_read_i[size].ref.push_back(ref_end);
//}


//void substrDriver(string& ref_str, string& read, vector<uint8_t>& index_data) {//!!!MAY HAVE TO CHANGE DATA TYPE FOR INDEXDATA!!!
//    vector<MatchedIndices> line_read_i;
//    sub_stuct read_substr;//new
//
//    while (getOpenSubStr(read, line_read_i, read_substr))//finds open string in read based on linei and places into subStr -
//        getLCS(ref_str, read_substr, line_read_i);//finds subStr's LCS in refStr, places indexes in refi and linei
//
//    //GROW REF STRING HERE!!!
//
//    sort(line_read_i.begin(), line_read_i.end() );//sorts line_ref_i so they go numerically according to line indexes
//
//    for (int i=0; i<line_read_i.size(); i++) {//place refi into index data and append end
//        //!!!USING TYPE CASTING FOR TESTING THIS, MAY WANT TO GO BACK AND CHANGE LINE_REF STRUCT!!!
//        index_data.push_back((uint8_t)line_read_i[i].ref[0]);
//        index_data.push_back((uint8_t)line_read_i[i].ref[1]);
//    }
//}



