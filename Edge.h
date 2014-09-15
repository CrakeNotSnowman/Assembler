//
//  Edge.h
//  CompressBases
//
//  Created by Austin Riffle on 1/9/14.
//  Copyright (c) 2014 Austin Riffle. All rights reserved.
//

#ifndef __CompressBases__Edge__
#define __CompressBases__Edge__

#include <iostream>
#include <string>
#include <memory>

class Node;

using std::string;
using std::shared_ptr;
using std::weak_ptr;

class Edge {
	string label;
	shared_ptr<Node> dest;

public:
    Edge();
    ~Edge();
    Edge(string label, shared_ptr<Node> dest);

	void   setLabel(string label);
    void   setDest(shared_ptr<Node> dest);
	shared_ptr<Node> getDest();
    string getLabel();

    bool operator < (const Edge& tmp) const {return (label < tmp.label);}
};

#endif /* defined(__CompressBases__Edge__) */
