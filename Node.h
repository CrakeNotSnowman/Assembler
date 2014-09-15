//
//  Node.h
//  CompressBases
//
//  Created by Austin Riffle on 1/9/14.
//  Copyright (c) 2014 Austin Riffle. All rights reserved.
//

#ifndef __CompressBases__Node__
#define __CompressBases__Node__

#include <iostream>
#include <array>
#include <vector>
#include <set>
#include <list>
#include "NodeLabel.h"
#include "EdgeBag.h"
#include "Edge.h"
#include <memory>

#define START_SIZE 0
#define INCREMENT  1

using std::vector;
using std::set;
using std::list;
using std::cout;
using std::endl;
using std::shared_ptr;
using std::weak_ptr;


class Node {
	vector<int> data;

	shared_ptr<EdgeBag> edges;
	weak_ptr<Node>    suffix;
	int  resultCount;//$1/22 figured we don't want this always = -1, so added that to constructor

	bool contains(int index);
	set<int> computeAndCacheCountRecursive();
	void addidx(int index);

protected:
	int computeAndCacheCount();

public://$1/22 changed all of these to be initialized in the Constructor
vector<NodeLabel> labels;
vector<int> nodeSuff;
    vector<int> nodePref;
	int C_v;
    int lastIdx;
	string prefix;
	int leaf_number;
	bool visited;
	bool visited1;
	bool visited2;
	bool visited3;
	bool visited4;
	bool visited5;
    weak_ptr<Node> parent;
	int   nodeDepth;




	Node();
    ~Node();

	void  printIncomingEdges();
	void  printOutgoingEdges();
	void  printNodeLabels();
    shared_ptr<Edge> getIncomingEdge();
	void     addRef(int index);
	int      getResultCount();
	void     addEdge(char ch, shared_ptr<Edge> e);
	shared_ptr<Edge>    getEdge(char ch);
	shared_ptr<EdgeBag> getEdges();
	shared_ptr<Node>    getSuffix();
	void     setSuffix(shared_ptr<Node> suffix);
	void getData(int numElements, bool b, vector<int> *l);
	vector<int> getData(int numElements, bool b);
    vector<int> getData(int numElements);
	vector<int> getData();
};

#endif /* defined(__CompressBases__Node__) */


































