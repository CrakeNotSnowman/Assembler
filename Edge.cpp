//
//  Edge.cpp
//  CompressBases
//
//  Created by Austin Riffle on 1/9/14.
//  Copyright (c) 2014 Austin Riffle. All rights reserved.
//

#include "Edge.h"

Edge::Edge() {
    this->label = "";
    this->dest  = shared_ptr<Node>();
}

Edge::~Edge() {
    dest.reset();
}

Edge::Edge(string label, shared_ptr<Node> dest)  {
	this->label = label;
	this->dest  = dest;
}

void Edge::setLabel(string label) {
	this->label = label;
}

void Edge::setDest(shared_ptr<Node> dest)  {
	this->dest = dest;
}

string Edge::getLabel() {
	return label;
}

shared_ptr<Node> Edge::getDest() {
	return dest;
}
