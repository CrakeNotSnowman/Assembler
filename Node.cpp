//
//  Node.cpp
//  CompressBases
//
//  Created by Austin Riffle on 1/9/14.
//  Copyright (c) 2014 Austin Riffle. All rights reserved.
//

#include "Node.h"
#include "NodeLabel.h"
using namespace std;
Node::Node() {
    //private:
	edges  = shared_ptr<EdgeBag> (new EdgeBag() );
	suffix = shared_ptr<Node> (nullptr);
    data   = vector<int>(0);//$1/22 add this, unsure if it is inititialized to 0 or a different size, come back to this
    nodeDepth = 0;

    //public:
    C_v = -1;
    prefix  = "";
    lastIdx = 0;
	leaf_number = -1;

    //$1/22 change this to vectors eventually? vector<bool> visited(6, false);?
    visited  = false;
	visited1 = false;
	visited2 = false;
	visited3 = false;
	visited4 = false;
	visited5 = false;

    parent = weak_ptr<Node>();
    resultCount = -1;
}

Node::~Node() {
    parent.reset();
    edges.reset();
    suffix.reset();
}

bool Node::contains(int index) {
	int low = 0;
	int high = lastIdx - 1;

	while (low <= high) {
		int mid = (low + high) / 2;
		int midVal = data[mid];
		if (midVal < index)
			low = mid + 1;
		else if (midVal > index)
			high = mid - 1;
		else
			return true;
	}
	return false;
}

set<int> Node::computeAndCacheCountRecursive() {
    set<int> ret;
    for (int num : data)
        ret.insert(num);
    for (shared_ptr<Edge> e : edges->getValues()) {
		for (int num : e->getDest()->computeAndCacheCountRecursive())
			ret.insert(num);
    }

    resultCount = (int)ret.size();
    return ret;
}

void Node::addidx(int index) {
	if ((unsigned int)lastIdx == data.size()) {
		vector<int> copy;
		copy.resize(data.size() + INCREMENT);
		for (int i : copy)
			data.push_back(i);
	}
	data[lastIdx++] = index;
}

int Node::computeAndCacheCount() {
	computeAndCacheCountRecursive();
	return resultCount;
}

void Node::printIncomingEdges() {
	shared_ptr<EdgeBag> eb (this->parent.lock()->edges);
	vector<shared_ptr<Edge>> e = eb->getValues();

	for (unsigned int i = 0; i<e.size(); i++)
		cout << "\tincoming: " << e[i]->getLabel() << endl;
}

shared_ptr<Edge> Node::getIncomingEdge() {
	shared_ptr<EdgeBag> eb (this->parent.lock()->edges);
	vector<shared_ptr<Edge> > e = eb->getValues();

	return e[0];
}

void Node::printOutgoingEdges() {
    vector<shared_ptr<Edge>> e = this->edges->getValues();
    for (unsigned int i=0; i<e.size(); i++)
        cout << "\toutgoing: " << e[i]->getLabel() << endl;
}

void Node::printNodeLabels() {
	for (NodeLabel k : labels)
		cout << "(" << k.index << "," << k.position << ")" << endl;
}



vector<int> Node::getData(int numElements) {

	vector<int> l;
	for (int num : this->data) {
        l.push_back(num);
	}
    vector<shared_ptr<Edge> > edgebag = edges->getValues();
	for (shared_ptr<Edge> e : edgebag) {
	    vector<int> tmp_data = e->getDest()->getData(-1);

        for (int num : tmp_data) {
            l.push_back(num);
        }
	}
	for(shared_ptr<Edge> e : edgebag){
        e.reset();
	}
	return l;
}


vector<int> Node::getData(){
	return getData(-1);
}

void Node::addRef(int index) {
	if (contains(index))
		return;

	addidx(index);
	// add this reference to all the suffixes as well
	shared_ptr<Node> iter (this->suffix.lock());
	while (iter != nullptr) {
		if (iter->contains(index))
			break;
		iter->addRef(index);
		iter = iter->suffix.lock();
	}
}

void Node::addEdge(char ch, shared_ptr<Edge> e) {
	edges->put(ch, e);
}

shared_ptr<Edge> Node::getEdge(char ch) {
	if (edges->get(ch) == nullptr)
		return nullptr;
	return edges->get(ch);
}

shared_ptr<EdgeBag> Node::getEdges() {
	return edges;
}

shared_ptr<Node> Node::getSuffix() {
	return suffix.lock();
}

void Node::setSuffix(shared_ptr<Node> suffix) {
	this->suffix = weak_ptr<Node>(suffix);
}





