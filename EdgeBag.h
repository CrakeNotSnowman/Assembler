//
//  EdgeBag.h
//  CompressBases
//
//  Created by Austin Riffle on 1/10/14.
//  Copyright (c) 2014 Austin Riffle. All rights reserved.
//
//  EdgeBag.h
//  CompressBases
//
//  Created by Austin Riffle on 1/10/14.
//  Copyright (c) 2014 Austin Riffle. All rights reserved.
//

#ifndef __CompressBases__EdgeBag__
#define __CompressBases__EdgeBag__

#include <iostream>
#include <vector>
#include <memory>
#include "Edge.h"

#define BSEARCH_THRESHOLD 6

using std::vector;
using std::shared_ptr;

class EdgeBag {
	vector<char> chars;
	vector<shared_ptr<Edge> > values;

	int search(char c);
	void sortArrays();

public:
    EdgeBag();
    ~EdgeBag();
	vector<shared_ptr<Edge> > getValues();
	void put(char c, shared_ptr<Edge>& e);
	shared_ptr<Edge> get(char c);
};

#endif /* defined(__CompressBases__EdgeBag__) */







