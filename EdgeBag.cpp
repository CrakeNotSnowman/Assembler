//
//  EdgeBag.cpp
//  CompressBases
//
//  Created by Austin Riffle on 1/10/14.
//  Copyright (c) 2014 Austin Riffle. All rights reserved.
//

#include "EdgeBag.h"

EdgeBag::EdgeBag() {
    chars  = vector<char>(0);
    values = vector<shared_ptr<Edge> >(0);
}

EdgeBag::~EdgeBag() {
    for (shared_ptr<Edge> i : values)
        i.reset();
}

int EdgeBag::search(char c) {
    if (chars.empty() )
        return -1;

    for (unsigned int i = 0; i < chars.size(); i++)
        if (c == chars[i])
            return i;

    return -1;
}

void EdgeBag::sortArrays() {
    for (unsigned int i = 0; i < chars.size(); i++)
        for (int j = i; j > 0; j--)
            if (chars[j-1] > chars[j]) {
                char swap = chars[j];
                chars[j] = chars[j-1];
                chars[j-1] = swap;

                shared_ptr<Edge> swapEdge = values[j];
                values[j] = values[j-1];
                values[j-1] = swapEdge;
            }
}


//public stuff:

vector<shared_ptr<Edge> > EdgeBag::getValues() {
    return values;
}

void EdgeBag::put(char c, shared_ptr<Edge>& e) {
    int idx = search(c);

    if (idx < 0) {//change chars
        chars.push_back(c);
        values.push_back(e);

        if ((int)chars.size() > BSEARCH_THRESHOLD)
            sortArrays();
    }

    else
        if(!values.empty())
            values[idx] = e;
}

shared_ptr<Edge> EdgeBag::get(char c) {
    int idx = search(c);
    if (idx == -1 || values.size() == 0)
        return nullptr;
    return values[idx];
}
