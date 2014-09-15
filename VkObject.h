//
//  Vk_Object.h
//  CompressBases
//
//  Created by Austin Riffle on 1/9/14.
//  Copyright (c) 2014 Austin Riffle. All rights reserved.
//  

#ifndef __CompressBases__Vk_Object__
#define __CompressBases__Vk_Object__

#include <iostream>

using std::string;

class VkObject{
    string substring;
    int nodeDepth;
public:
    VkObject();
	VkObject(int a, string b);
    
    string getSubstr();
    int getNodeDepth();
    void setSubstr(string substr);
    void setNodeDepth(int nd);
};

#endif /* defined(__CompressBases__Vk_Object__) */
