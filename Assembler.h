#include <iostream>
#include <fstream>


#include "GeneralizedSuffixTree.h"

class Assembler{


public:
    //GeneralizedSuffixTree assmTree;
    void populateIndices(GeneralizedSuffixTree* tree, vector<string> set);
    vector<int> chain;
    vector<int> chain_depth;
    vector<vector<shared_ptr<Node> > > getNodeDepth(GeneralizedSuffixTree* tree, int maxDepth);
    //void assemble(vector<string> testSet);
    int assemble(vector<string> testSet, string out_filename);
    void printMatches(vector<string> reads);
    vector<string> printAssembled(vector<string> reads);

    void printMatches(vector<string> reads, int a);
    Assembler();
    ~Assembler();
};
