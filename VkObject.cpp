//
//  Vk_Object.cpp
//  CompressBases
//
//  Created by Austin Riffle on 1/9/14.
//  Copyright (c) 2014 Austin Riffle. All rights reserved.
//  

#include "VkObject.h"

VkObject::VkObject() {
    nodeDepth = 0;
    substring = "";
}

VkObject::VkObject(int a, string b) {
	nodeDepth = a;
	substring = b;
}

string VkObject::getSubstr() {
    return substring;
}
int VkObject::getNodeDepth() {
    return nodeDepth;
}

void VkObject::setSubstr(string substr) {
    substring = substr;
}

void VkObject::setNodeDepth(int nd) {
    nodeDepth = nd;
}