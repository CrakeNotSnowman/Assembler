
#include "Assembler.h"
#include "fns.h"
#include <algorithm>
#include <time.h>

using namespace std;



void Assembler::printMatches(vector<string> reads, int a){
     for(unsigned int i=0; i<this->chain.size(); i++){

        if(this->chain[i] == 0){

            cout << reads[i]<< endl << "NO MATCH" <<endl;

        }else{

        cout << reads[this->chain[i]] + reads[i].substr(this->chain_depth[i] - 1, reads[i].size() - this->chain_depth[i]) << endl;

        }

        cout << endl;
    }

    return;
}

vector<string> Assembler::printAssembled(vector<string> reads){
    vector<string> comb_reads;
    bool isPrint = false;
    if(isPrint) cout << "---------------------------------------------------------" << endl;
    for(unsigned int i=0; i<this->chain.size(); i++){
        if(isPrint)cout << "i : " << i << endl;
        if(this->chain[i] == 0){



        }else{


        if(isPrint)cout << "reads[" << i << "]" << reads[this->chain[i]] << " + " << reads[i].substr(this->chain_depth[i], reads[i].size() - this->chain_depth[i]) << endl;
        reads[i] = reads[this->chain[i]] + reads[i].substr(this->chain_depth[i], reads[i].size() - this->chain_depth[i]);



        if(isPrint)cout << "chain[" << i << "] = " << chain[chain[i]]<<endl;

        if(chain[i] >= 0) reads[chain[i]] = "";

        int tmp_chain = chain[i];

        chain[i] = chain[chain[i]];
        chain_depth[i] = chain_depth[tmp_chain];
        chain[tmp_chain] = 0;
        i--;
        }

    }

    if(isPrint){for(unsigned int i = 0; i<this->chain.size(); i++){
        cout << i <<" "<< this->chain[i] <<" "<< this->chain_depth[i] << endl;
    }}
    for(string s : reads){
        if(s != "") comb_reads.push_back(s);
    if(isPrint)cout << s << endl;
    }
    return comb_reads;
    //printMatches(reads);
}





int Assembler::assemble(vector<string> testSet, string out_filename){

clock_t t1,t2;
    t1=clock();
    //code goes here


    GeneralizedSuffixTree* assmTree = new GeneralizedSuffixTree();
    unsigned int maxDepth = 0;

    for(unsigned int i=0; i<testSet.size(); i++){
        if(testSet[i].length() > maxDepth) maxDepth = testSet[i].length();
        assmTree->put(testSet[i], i);
    }

    //cout << "setting node depths" << endl;
    setNodeDepths(assmTree->root, 0);
    int numStrings = testSet.size();
    //cout << "populating indices" << endl;
    populateIndices(assmTree, testSet);

    int curr_depth = maxDepth;
    vector<vector<shared_ptr<Node> > > nDepths;
    nDepths = getNodeDepth(assmTree, maxDepth);


    vector<int> *chain_depth_tmp = new vector<int>(testSet.size(),0);
    vector<int> *chain_tmp = new vector<int>(testSet.size(), 0);
    vector<int> *wrap_tmp = new vector<int>(testSet.size(), 0);
    vector<int> u_data;
    //cout << "main loop" << endl;
    int count = 0;

    shared_ptr<Node> u;
    while(numStrings > 1 && curr_depth >= 50){

        while(nDepths[curr_depth].size() == 0){
            curr_depth--;
            }


        u = nDepths[curr_depth][0];

        u_data = u->getData();


        for(NodeLabel nl : u->labels){
            if(nl.position == 1) u->nodePref.push_back(nl.index);  /// need to be 1??
            else u->nodeSuff.push_back(nl.index);
        }


        int i = -1;
        int j = -1;
        vector<int> bad_x;

        for(int x : u->nodePref){
            for(unsigned int m=0; m<u->nodeSuff.size(); m++){
                int y = u->nodeSuff[m];
                if(chain_tmp->at(x) != 0 && x != y){
                    bad_x.push_back(x);
                } else if(chain_tmp->at(x) == 0 && wrap_tmp->at(x) != y && x != y){
                    i = x;
                    j = y;
                }
            }
        }


        for(unsigned int k=0; k<u->nodeSuff.size(); k++){
            for(unsigned int m=0; m<bad_x.size(); m++){
                if(u->nodeSuff[k] == bad_x[m]){ u->nodeSuff.erase(u->nodeSuff.begin() + k-1); /// k not m ??
                k--;}
            }
        }

        if(i!= -1 && j!= -1 && chain_tmp->at(i) != j){
            int uNodeSize = u->nodeSuff.size();
            vector<int> tmpSuff;
            vector<int> tmpPref;
            for(unsigned int k=0; k < uNodeSize; k++){
               // if(u->nodeSuff[k] == i) u->nodeSuff.erase(u->nodeSuff.begin()+i-1);  /// This doesn't matter appa
               if(!(u->nodeSuff[k] = i)) tmpSuff.push_back(u->nodeSuff[k]);

            }


            for(unsigned int k=0; k < uNodeSize; k++){
                //if(u->nodePref[k] == i){
                   // u->nodePref.erase(u->nodePref.begin()+k-1);
                   if(!(u->nodePref[k] = i)) tmpPref.push_back(u->nodePref[k]);

               // }
            }
            u->nodeSuff = tmpSuff;
            u->nodePref = tmpPref;

             chain_tmp->at(i) = j;


             chain_depth_tmp->at(i) = curr_depth;



             if(wrap_tmp->at(i) != 0) wrap_tmp->at(wrap_tmp->at(i)) = wrap_tmp->at(j);
             if(wrap_tmp->at(j)!=0) wrap_tmp->at(wrap_tmp->at(j)) = wrap_tmp->at(i);


        } else {

           shared_ptr<Node> v = u->parent.lock();


            for(int k: u->nodePref){
               v->nodePref.push_back(k);
            }

            vector<shared_ptr<Node> >::iterator position = find(nDepths[curr_depth].begin(), nDepths[curr_depth].end(), u);
            if (position != nDepths[curr_depth].end())
                nDepths[curr_depth].erase(position);

        }

    }
    t2=clock();
    this->chain = *chain_tmp;
    this->chain_depth = *chain_depth_tmp;


//    delete chain_tmp;
//    delete chain_depth_tmp;
//    delete wrap_tmp;
//    delete assmTree;
float diff ((float)t2-(float)t1);



    //cout << "Assembly done - " << diff/CLOCKS_PER_SEC << " seconds "<<endl;

    vector<string> result = printAssembled(testSet);

    ofstream out_stream(out_filename.c_str(), ofstream::binary);
//    if (out_stream.fail())
//    {
//        cerr << "Error creating" << out_filename << endl;
//        cerr << "Please enter a valid output filename\n";
//        exit(-2);
//    }
    int j=0;
    int line_num=0;
    int smax=0;
    for(string s : result){
        if(s.size() > smax){
         line_num = j;
         smax = s.size();
        }
        j++;
        out_stream << s << endl;
    }
//cout << "Results printed!" << endl;
//cout << "longest line: " << line_num << " - " << smax << endl;
 //exit(0);

    return 0;
}

void Assembler::populateIndices(GeneralizedSuffixTree* tree, vector<string> set){


    for(unsigned int i=0; i<set.size(); i++){

        string temp = set[i];

        for(unsigned int j=0; j<temp.length(); j++){

            NodeLabel nlabel;
            nlabel.index = i;
            nlabel.position = j+1;
            shared_ptr<Node> n = tree->searchNode(temp.substr(j,temp.length()));
            n->labels.push_back(nlabel);

        }
    }
    return;
}

vector<vector<shared_ptr<Node> > > Assembler::getNodeDepth(GeneralizedSuffixTree* tree, int maxDepth) {

    shared_ptr<Node> v = tree->root;
    stack<shared_ptr<Node> > NodeStack;
    stack<shared_ptr<Node> > ResultStack;
    NodeStack.push(v);
    vector<vector<shared_ptr<Node> > > nvect;
    for(int i=0; i<=maxDepth; i++){
        nvect.push_back(vector<shared_ptr<Node> >());
    }

    while(!NodeStack.empty()){
        shared_ptr<Node> u = NodeStack.top();
        NodeStack.pop();
        ResultStack.push(u);

        if(!u->visited){
            u->visited = true;

            nvect[u->nodeDepth].push_back(u);
            /// for each neighbour w of u, push u->NodeStack
            shared_ptr<EdgeBag> eb = u->getEdges();
            vector<shared_ptr<Edge> > e = eb->getValues();
            //sort(e.begin(), e.end());

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
        vector<shared_ptr<Edge> > e = eb->getValues();
        if(e.empty()) {
            n->leaf_number = j;
            j++;
        }
    }

    return nvect;
}

void Assembler::printMatches(vector<string> reads){
     for(unsigned int i=0; i<this->chain.size(); i++){

        if(this->chain[i] == 0){

            cout << reads[i]<< endl << "NO MATCH" <<endl;

        }else{

        cout << reads[ this->chain[i] ];

        for(unsigned int j=0; j<reads[i].size() - this->chain_depth[i]; j++) cout << '-';

        cout << endl;

        for(unsigned int j=0; j<reads[i].size() - this->chain_depth[i]; j++) cout << '-';

        cout << reads[i] << endl;

        cout << reads[this->chain[i]] + reads[i].substr(this->chain_depth[i] - 1, reads[i].size() - this->chain_depth[i]) << endl;
        }

        cout << endl;
    }
}

Assembler::Assembler(){}
Assembler::~Assembler(){}


















