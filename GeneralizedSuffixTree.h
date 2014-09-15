//
//  GeneralizedSuffixTree.h
//  CompressBases
//
//  Created by Austin Riffle on 1/9/14.
//  Copyright (c) 2014 Austin Riffle. All rights reserved.
//

#ifndef __CompressBases__GeneralizedSuffixTree__
#define __CompressBases__GeneralizedSuffixTree__

#include <utility>
#include <stack>
#include "Node.h"
#include "VkObject.h"

struct MatchedIndices;

using std::pair;
using std::stack;

const string SP_CHAR[] = {
    "0", "1", "2", "3",
    "4", "5", "6", "7",
    "8", "9", "!", "@",
    "#", "$", "%", "^",
    "&", "*", "<", ">"
};

class GeneralizedSuffixTree{
    int last;


    vector<VkObject> Vk;
    vector<vector<Node>> Li;

   pair<bool, shared_ptr<Node> > testAndSplit(shared_ptr<Node>& inputNode, const string& stringPart, const char t, const string& remainder, const int value);
	pair<shared_ptr<Node>, string> canonize(shared_ptr<Node> inputNode, const string& inputstr);
	pair<shared_ptr<Node>, string> update(shared_ptr<Node>& inputNode, const string& stringPart, const string& rest, const int value);
	string safeCutLastChar(string seq);

    void depthFirstSearch();
    void setCv();
    void setStringDepth();
    void setVk();
    void setLi();

public:
shared_ptr<Node> root;
	shared_ptr<Node> activeLeaf;
    GeneralizedSuffixTree();
     GeneralizedSuffixTree(const vector<string>& vec_str);
    GeneralizedSuffixTree(const string& str0, const string& str1);
    ~GeneralizedSuffixTree();
	list<int> search(string word);
	list<int> search(string word, int results);

	//ResultInfo searchWithCount(string word, int to);   Do we need this function?
	shared_ptr<Node> searchNode(string word);
	void put(const string& word, int index);
    shared_ptr<Node> getRoot();
    int getLast();
	int computeCount();
    class ResultInfo{};
    void printNodes();
    void printRootEdges();
    void printOneDeepEdges();

    void printVkTable();
    void setVkTable();
    string getBestLCS(string refStr, double thresh, int min);
    vector<VkObject> getVk();
    string getLCS_k_strs(int k);
};


#endif /* defined(__CompressBases__GeneralizedSuffixTree__) */
