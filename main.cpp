
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdint>
#include <utility>
#include <algorithm>
#include "GeneralizedSuffixTree.h"
#include "Assembler.h"
#include "fns.h"
#include <Python.h>

using namespace std;

// Modefied code by Keith Murray, 
// Modification maded base on http://intermediate-and-advanced-software-carpentry.readthedocs.org/en/latest/c++-wrapping.html


int main(string in_filename, string out_filename)
{
    //string in_filename;

    //string out_filename;
    int n_reads = 5000;

    /*
    if(argc > 1){
        in_filename = argv[1];
        out_filename = argv[2];

    }
    else{
        in_filename = "/home/kmurray/Desktop/frags.txt";//change to cmd line eventually
        out_filename = "/home/kmurray/Desktop/AssembledContigs";
        }
    */
    //cout << "Number of reads: " << endl;
    //cin >> n_reads;
    n_reads = 10000000;

    ifstream in_stream(in_filename.c_str());
    if (in_stream.fail())
    {
        cerr << "Error opening " << in_filename << endl;
        cerr << "Please enter a valid input file\n";
        //exit(-1);
    }






    vector<string> vec_reads = getVecReads(in_stream, n_reads);
    //cout << "Number of contigs read: " << vec_reads.size() << endl;


    Assembler assm;// = new Assembler();
    assm.assemble(vec_reads, out_filename);

    return 0;


}

static PyObject * JBGenSuffTree(PyObject * self, PyObject * args)
{
  char * in_filename;
  char * out_filename;
  int result;
  PyObject * ret;

  // parse arguments
  if (!PyArg_ParseTuple(args, "ss", &in_filename, &out_filename)) {
    return NULL;
  }

  // run the actual function
  main(in_filename, out_filename);
  
  ret = PyInt_FromLong(0);
  return ret;
  //Py_RETURN_NONE;// Py_BuildValue("i", result);
}

static PyMethodDef JBGenTreeAlignmentMethods[] = {
	{"GenTreeAlign", JBGenSuffTree, METH_VARARGS,"Generalized Suffix tree based alignment"},
	{NULL,NULL,0,NULL}
};

PyMODINIT_FUNC
initJBGenSuffTree(void)
{
	Py_InitModule("JBGenSuffTree", JBGenTreeAlignmentMethods);
}

