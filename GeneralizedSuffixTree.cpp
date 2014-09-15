//
//  GeneralizedSuffixTree.cpp
//  CompressBases
//
//  Created by Austin Riffle on 1/9/14.
//  Copyright (c) 2014 Austin Riffle. All rights reserved.
//

#include "GeneralizedSuffixTree.h"
#include <string>
#include <list>
#include <utility>
#include <algorithm>
#include "Node.h"

using namespace std;


GeneralizedSuffixTree::GeneralizedSuffixTree() {
    last = -1;
    root = shared_ptr<Node>(new Node() );
    activeLeaf   = root;
    root->parent = weak_ptr<Node>();
}

shared_ptr<Node> GeneralizedSuffixTree::searchNode(string word){


    shared_ptr<Node> currentNode = root;
    shared_ptr<Edge> currentEdge = nullptr;

    for(unsigned int i=0; i<word.length(); i++){
        char ch = word[i];
        /// follow edge corresponding to ch
        currentEdge = currentNode->getEdge(ch);
        if(nullptr == currentEdge){
            return nullptr;
        } else {
            string label = currentEdge->getLabel();
            int lenToMatch = min(word.length() - i, label.length());
            if(!(word.substr(i, lenToMatch) == (label.substr(0,lenToMatch)))){
                return nullptr;
            }
            if(label.length() >= word.length() - i){
                return currentEdge->getDest();
            } else{
            /// advance to next node
                currentNode = currentEdge->getDest();
                i += lenToMatch - 1;
            }
        }
    }

    return nullptr;
}
GeneralizedSuffixTree::GeneralizedSuffixTree(const vector<string>& vec_str) {
    last = -1;
    root = shared_ptr<Node>(new Node() );
    activeLeaf   = root;
    root->parent = weak_ptr<Node>();

    for(unsigned int i=0; i<vec_str.size(); i++)
        put(vec_str[i] + SP_CHAR[i], i);
}
GeneralizedSuffixTree::GeneralizedSuffixTree(const string& str0, const string& str1) {
    last = -1;
    root = shared_ptr<Node>(new Node() );
    activeLeaf   = root;
    root->parent = weak_ptr<Node>();

    put(str0 + SP_CHAR[0], 0);
    put(str1 + SP_CHAR[1], 1);
}
GeneralizedSuffixTree::~GeneralizedSuffixTree() {
    root.reset();
    activeLeaf.reset();
}

void GeneralizedSuffixTree::printRootEdges() {
    root->printOutgoingEdges();
}

void GeneralizedSuffixTree::printOneDeepEdges() {
    for(shared_ptr<Edge> e : root->getEdges()->getValues())
        e->getDest()->printOutgoingEdges();
}

//private:
void GeneralizedSuffixTree::put(const string& word, int index) {
    last = index;
	activeLeaf  = root;
    string text = "";
	shared_ptr<Node> s (root);

	// iterate over the string, one char at a time
	for (unsigned int i = 0; i < word.length(); i++){
		text += word.at(i);

        pair<shared_ptr<Node>, string> active = update(s, text, word.substr(i, word.size()-i), index);// line 7: update the tree w/ the new transitions due to this new char
        active = canonize(active.first, active.second);//line 8: make sure the active pair is canonical

		s    = active.first;
		text = active.second;
	}

	if ( (activeLeaf->getSuffix() == nullptr) && (activeLeaf != root) && (activeLeaf != s) )
		activeLeaf->setSuffix(s);
}

pair<shared_ptr<Node>, string> GeneralizedSuffixTree::update(shared_ptr<Node>& inputNode, const string& stringPart, const string& rest, int value) {
    shared_ptr<Node>  s (inputNode);
	string tempstr = stringPart;
	char   newChar = stringPart.at(stringPart.length() - 1);

	shared_ptr<Node> oldroot (root);
    pair<bool, shared_ptr<Node> > ret = testAndSplit(s, tempstr.substr(0, tempstr.length() - 1), newChar, rest, value);

    //bool endpoint = ret.first;
    shared_ptr<Node> r (ret.second);
    shared_ptr<Node> leaf(nullptr);

	while (ret.first == false) {
		shared_ptr<Edge> tempEdge (r->getEdge(newChar));
		if (tempEdge != nullptr){
			// such a node is already present. This is one of the main differences fro Ukkonen's case:
			// the tree can contain deeper nodes at this stage becuase different string were added by previous iterations
			leaf = tempEdge->getDest();
		}
		else{
			// must build a new leaf
			leaf = shared_ptr<Node>(new Node() );
			leaf->addRef(value);
			shared_ptr<Edge> newedge = shared_ptr<Edge>(new Edge(rest, leaf));
			r->addEdge(newChar, newedge);
		}

		// update suffix link for newly created leaf
		if (activeLeaf != root)
			this->activeLeaf->setSuffix(leaf);

		this->activeLeaf = leaf;

		if (oldroot != root)
			oldroot->setSuffix(r);

		oldroot = r;
        shared_ptr<Node> tmp (s->getSuffix() );

		if (tmp == nullptr)
            tempstr = tempstr.substr(1, tempstr.size() - 1);
		else{
			pair<shared_ptr<Node>, string> canret = canonize(s->getSuffix(), safeCutLastChar(tempstr));
			s = canret.first;
			tempstr = (canret.second) + tempstr.at(tempstr.length() - 1);
		}
		ret = testAndSplit(s, safeCutLastChar(tempstr), newChar, rest, value);
		//endpoint = ret.first;
        r = ret.second;
	}

	if (oldroot != root)
		oldroot->setSuffix(r);
	oldroot = root;

	return pair<shared_ptr<Node>, string>(s, tempstr);
}

// Tests whether the string stringPart + t is contained in the subtree that has inputs as root.
// If that's not the case, and there exists a path of edges e1, e2, ... such that
// e1.label + e2.label + ... + $end = stringPart
// and there is an edge g such that
// g.label = stringPart + rest
//
// Then g will be split in two different edges, one having $end as label, and the other one
// having rest as label

pair<bool, shared_ptr<Node> > GeneralizedSuffixTree::testAndSplit(shared_ptr<Node>& inputNode, const string& stringPart, const char t, const string& remainder, const int value) {
    // descend the tree as far as possible
	pair<shared_ptr<Node>, string> ret = canonize(inputNode, stringPart);
	shared_ptr<Node>  s (ret.first);
	string str = ret.second;

	if (!str.empty()) {
		shared_ptr<Edge> g (s->getEdge(str.at(0)));
        string label = g->getLabel();

		if ((label.length() > str.length() ) && (label.at(str.length()) == t) ) {//must see whether "str" is substr of the label of an edge
			return pair<bool, shared_ptr<Node>>(true, s);
		} else {
            int n = (int)str.length();
			string newlabel = label.substr(n, label.length()-n);
			// build a new node, edge
			shared_ptr<Node> r = shared_ptr<Node>(new Node());
			shared_ptr<Edge> newedge = shared_ptr<Edge>(new Edge(str, r));

			g->setLabel(newlabel);

			// link s -> r
			r->addEdge(newlabel.at(0), g);
			s->addEdge(str.at(0), newedge);

			return pair<bool, shared_ptr<Node>>(false, r);
		}
	} else {
		shared_ptr<Edge> e (s->getEdge(t));
		if (e == nullptr) {
			return pair<bool, shared_ptr<Node>>(false, s);
        } else {
			if (!remainder.compare(e->getLabel()) ) {
				e->getDest()->addRef(value);
				return pair<bool, shared_ptr<Node>>(true, s);
			} else if (!remainder.compare(0, e->getLabel().length(), e->getLabel()) ) {//here
				return pair<bool, shared_ptr<Node>>(true, s);
            } else if (!e->getLabel().compare(0, remainder.length(), remainder)) {
				shared_ptr<Node> newNode = shared_ptr<Node>(new Node());
				newNode->addRef(value);

				shared_ptr<Edge> newEdge = shared_ptr<Edge>(new Edge(remainder, newNode) );
                int n = (int)remainder.length();
				e->setLabel(e->getLabel().substr(n, e->getLabel().length() - n));

                newNode->addEdge(e->getLabel().at(0), e);
				s->addEdge(t, newEdge);

				return pair<bool, shared_ptr<Node>>(false, s);
			} else
				return pair<bool, shared_ptr<Node>>(true, s);
		}
	}
}

pair<shared_ptr<Node>, string> GeneralizedSuffixTree::canonize(shared_ptr<Node> inputNode, const string& inputstr) {
	if (inputstr == "")
		return pair<shared_ptr<Node>, string>(inputNode, inputstr);
	else {
		shared_ptr<Node> currentNode (inputNode);
		string str = inputstr;
		shared_ptr<Edge>  g (inputNode->getEdge(str.at(0)));

		// descend the tree as long as proper label is found
		while (g != nullptr && !str.compare(0, g->getLabel().length(), g->getLabel()) ){
            int n = (int)g->getLabel().length();
			str = str.substr(n, str.length() - n);
			currentNode = g->getDest();
			if (str.length() > 0)
				g = currentNode->getEdge(str.at(0));
		}
		return pair<shared_ptr<Node>, string>(currentNode, str);
	}
}

shared_ptr<Node> GeneralizedSuffixTree::getRoot() {
	return root;
}

string GeneralizedSuffixTree::safeCutLastChar(string seq) {
    if (seq.length() == 0)
        return "";
    return seq.substr(0, seq.length() - 1);
}

int GeneralizedSuffixTree::getLast() {
    return last;
}

void GeneralizedSuffixTree::depthFirstSearch() {
    shared_ptr<Node> v = this->root;
    stack<shared_ptr<Node>> NodeStack;
    stack<shared_ptr<Node>> ResultStack;
    NodeStack.push(v);

    while(!NodeStack.empty()){
        shared_ptr<Node> u = NodeStack.top();
        NodeStack.pop();
        ResultStack.push(u);

        if(!u->visited){
            u->visited = true;

            /// for each neighbour w of u, push u->NodeStack
            shared_ptr<EdgeBag> eb = u->getEdges();
            vector<shared_ptr<Edge>> e = eb->getValues();
            sort(e.begin(), e.end());

            for(shared_ptr<Edge> edge : e)
                NodeStack.push(edge->getDest());

        }
    }

    // fills in the leaf numbers
    int j=1;
    while(!ResultStack.empty()) {
        shared_ptr<Node> n = ResultStack.top();
        ResultStack.pop();
        shared_ptr<EdgeBag> eb = n->getEdges();
        vector<shared_ptr<Edge>> e = eb->getValues();
        if(e.empty()) {
            n->leaf_number = j;
            j++;
        }
    }
}

void GeneralizedSuffixTree::setCv() {
    shared_ptr<Node> v = this->root;
    stack<shared_ptr<Node>> NodeStack;
    NodeStack.push(v);

    while (!NodeStack.empty()) {
        shared_ptr<Node> u = NodeStack.top();
        NodeStack.pop();

        if (!u->visited1) {
            u->visited1 = true;

            // set the C(v) here
            set<int> coll;
            shared_ptr<EdgeBag> eb = u->getEdges(); //
            vector<shared_ptr<Edge>> e = eb->getValues();
            for (shared_ptr<Edge> edge : e)
                for (int i : edge->getDest()->getData() )
                    coll.insert(i);

            u->C_v = (int)coll.size();

            for(shared_ptr<Edge> edge : e)
                NodeStack.push(edge->getDest());
        }
    }
}

void GeneralizedSuffixTree::setStringDepth() {
    shared_ptr<Node> v = this->root;
    v->nodeDepth = 0;
    stack<shared_ptr<Node>> NodeStack;
    NodeStack.push(v);

    while(!NodeStack.empty()){
        shared_ptr<Node> u = NodeStack.top();
        NodeStack.pop();

        if(!u->visited4){
            u->visited4 = true;
            /// for each neighbour w of u, push u->NodeStack
            shared_ptr<EdgeBag> eb = u->getEdges();
            vector<shared_ptr<Edge>> e = eb->getValues();

            for(shared_ptr<Edge> edge : e){
                NodeStack.push(edge->getDest());
                edge->getDest()->nodeDepth = u->nodeDepth + (int)edge->getLabel().length();
            }

        }
    }
}

void GeneralizedSuffixTree::setLi() {
    for (int i=0; i<last+1; i++)
        Li.push_back(vector<Node>());
    /// DFS to hit every node

    shared_ptr<Node> v = this->root;
    stack<shared_ptr<Node>> NodeStack;
    stack<shared_ptr<Node>> ResultStack;
    NodeStack.push(v);

    while(!NodeStack.empty()){
        shared_ptr<Node> u = NodeStack.top();
        NodeStack.pop();
        ResultStack.push(u);

        if(!u->visited3){
            u->visited3 = true;

            if(u->getEdges()->getValues().size() == 0){// is a leaf  *** trying coll as a vector, previously was a collection.
                vector<int> coll = u->getData();
                for(int i : coll)
                    this->Li[i].push_back(*u);
            }

            shared_ptr<EdgeBag> eb = u->getEdges();
            vector<shared_ptr<Edge>> e = eb->getValues();

            for(shared_ptr<Edge> edge : e)
                NodeStack.push(edge->getDest());
        }
    }

}

void GeneralizedSuffixTree::printVkTable() {
    if (Vk.empty()) {
        cout << "WARNING: VK HAS NOT YET BEEN SET OR THE TREE IS EMPTY\n";
        cout << "Would you like to set Vk? If the tree is empty an\n";
        cout << "empty vector will still be returned and nothing printed.\n";
        cout << "If N is selcted an empty vector will also be printed.\n";
        cout << "(Y/N)?\n";
        string ans;
        cin >> ans;
        if (ans[0] == 'Y' || ans[0] == 'y')
            setVkTable();
    }
//    for(int i=0; i<Vk.size()-1; i++) {
//        if(Vk[i].getNodeDepth() < Vk[i+1].getNodeDepth()) {
//            Vk[i] = Vk[i+1];
//            i-=2;
//        }
//    }

    for(unsigned int i=2; i<Vk.size(); i++)
        cout << "V(" << i << ") = " << Vk[i].getNodeDepth() << "\t" << Vk[i].getSubstr() << endl;
}

void GeneralizedSuffixTree::setVk() {
    this->Vk = vector<VkObject>(last+2);//number of strings + 1 (+1 bc it is zero indexed)

    shared_ptr<Node> v = this->root;
    v->nodeDepth = 0;
    stack<shared_ptr<Node>> NodeStack;
    NodeStack.push(v);

    while(!NodeStack.empty()){
        shared_ptr<Node> u = NodeStack.top();
        NodeStack.pop();

        if(!u->visited5){
            u->visited5 = true;
            if(u->nodeDepth > Vk[u->C_v].getNodeDepth()) {
                VkObject tmp_vk(u->nodeDepth, u->prefix);
                Vk[u->C_v] = tmp_vk;
            }
            /// for each neighbour w of u, push u->NodeStack

            shared_ptr<EdgeBag> eb = u->getEdges();
            vector<shared_ptr<Edge>> e = eb->getValues();

            for(shared_ptr<Edge> edge : e){
                NodeStack.push(edge->getDest());
                edge->getDest()->nodeDepth = u->nodeDepth + (int)edge->getLabel().length();
                edge->getDest()->prefix    = u->prefix + edge->getLabel();
            }
            u->prefix = "";
        }
    }
    for(unsigned int i=0; i<Vk.size()-1; i++) {
        if(Vk[i].getNodeDepth() < Vk[i+1].getNodeDepth()) {
            Vk[i] = Vk[i+1];
            i-=2;
        }
    }
}

vector<VkObject> GeneralizedSuffixTree::getVk() {
    if (Vk.empty()) {
        cout << "WARNING: VK HAS NOT YET BEEN SET OR THE TREE IS EMPTY\n";
        cout << "Would you like to set Vk? If the tree is empty an\n";
        cout << "empty vector will still be returned. If N is selcted\n";
        cout << "an empty vector will also be returned. (Y/N)?\n";
        string ans;
        cin >> ans;
        if (ans[0] == 'Y' || ans[0] == 'y')
            setVkTable();
    }
    return this->Vk;
}

void GeneralizedSuffixTree::printNodes() {
    shared_ptr<Node> v = this->root;
    stack<shared_ptr<Node>> NodeStack;
    NodeStack.push(v);

    while(!NodeStack.empty()){
        shared_ptr<Node> u = NodeStack.top();
        NodeStack.pop();

        if(!u->visited2){
            u->visited2 = true;

            cout << "C(v) = " << u->C_v << endl;
            u->printOutgoingEdges();
            cout << "node depth = " << u->nodeDepth << endl;
            set<int> coll;
            for (int i : u->getData())
                coll.insert(i);
            for(int i : coll)
                cout << i << ",";
            cout << endl;

            shared_ptr<EdgeBag> eb = u->getEdges();
            vector<shared_ptr<Edge>> e = eb->getValues();
            for(shared_ptr<Edge> edge : e)
                NodeStack.push(edge->getDest());
        }
    }
}

void GeneralizedSuffixTree::setVkTable() {
    depthFirstSearch();
    setLi();
    setCv();
    setStringDepth();
    setVk();
}

string GeneralizedSuffixTree::getBestLCS(string refStr, double thresh, int min) {
    setVkTable();//need to put correction functionality in setVk not printVk
    vector<string> vkStr(Vk.size()-2);
    for (unsigned int i=2; i<Vk.size(); i++)
        vkStr[i-2] = Vk[i].getSubstr();
    for (unsigned int i=(int)vkStr.size() - 1; i>=0; i--) {
        if ((int)(vkStr[i].size()) >= min) {
            GeneralizedSuffixTree gst(refStr, vkStr[i]);
            string lcs = gst.getLCS_k_strs(2);
            if ((double)lcs.size()/vkStr[i].size() < thresh)//if the percent amount of new info is lower then the allowed threshold, return curr vkStr
                return vkStr[i];
        }
    }

    return "";
}

string GeneralizedSuffixTree::getLCS_k_strs(int k) {//k must be <= then last + 1
    depthFirstSearch();
    setLi();
    setCv();
    setStringDepth();
    setVk();

    return this->getVk()[k].getSubstr();
}












