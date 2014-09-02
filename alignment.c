#include "suffix_tree.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <Python.h>
/* Keith Murray 082013
 * Suffix Tree Alignment Beta 
 */


/****************************************************************************
			# of Substrings
****************************************************************************/


/* Added function prototype
 * (Not good to call functions from above) -jb
 */
inline int * ST_FindSubstring_Overlap(SUFFIX_TREE* tree, char  W[], DBL_WORD P);

/// Changed function type from (int) to (int *) -jb ///
inline int * ACS(char String_A[], char String_B[], int size_a, int size_b){
	int* i;
	SUFFIX_TREE* tree = ST_CreateTree(String_A, size_a);
	i = ST_FindSubstring_Overlap(tree, String_B, size_b);
	ST_DeleteTree(tree);
	
return i;
}
/******************************************************************************/
/*
   ST_FindSubstring :
   See suffix_tree.h for description.

	Modificaton:
   This function now returns an array of length two which shows the location
   in string A of the overlap (between string A and B) and the depth in String
   B of said overlap. (in alignment[0] and alignment[1] respectively)
   
*/

/* changed function type from (int) to (int*)
	also changed retrun types	-jb
*/
inline int * ST_FindSubstring_Overlap(
                      /* The suffix array */
                      SUFFIX_TREE*    tree,      
                      /* The substring to find */
                      char  W[],         
                      /* The length of W */
                      DBL_WORD        P)         
{
   /* Starts with the root's son that has the first character of W as its
      incoming edge first character */
   NODE* node   = find_son(tree, tree->root, W[0]);
   DBL_WORD k,j = 0, node_label_end;
   static int alignment[2];


   /* Scan nodes down from the root untill a leaf is reached or the substring is
      found */
   while(node!=0)
   {
      k=node->edge_label_start;
      node_label_end = get_node_label_end(tree,node);
      
      /* Scan a single edge - compare each character with the searched string */
      while(j<P && k<=node_label_end && tree->tree_string[k] == W[j])
      {
         j++;
         k++;
#ifdef STATISTICS
         counter++;
#endif
      }
      
      /* Checking which of the stopping conditions are true */
      if(j == P)
      {
         /* W was found - it is a substring. Return its path starting index */
	     alignment[0] = node->path_position-1;
	 alignment[1] = j-1;
         return alignment; //< EXIT POINT
      }
      else if(k > node_label_end){
         /* Current edge is found to match, continue to next edge */
	 node = find_son(tree, node, W[j]);
	 }
      else
      {
         /* One non-matching symbols is found - W is not a substring */
	 /*
	 IF string B is not completely within string A, (which is most likly)
	 */
	 alignment[0] = node->path_position-1;
	 alignment[1] = j-1;
         return alignment;//< EXIT POINT
      }
   }
   //This is thrown if a new character is seen in the very begining of string B that does not exist in the tree
   alignment[0] = -1;
   alignment[1] = -1;
   return alignment;//< EXIT POINT
}

/******************************************************************************/

/****************************************************************************
			MAIN
****************************************************************************/
static PyObject *
KM_TreeAlign(PyObject *self, PyObject *args)
{	
	unsigned long size_a;
	unsigned long size_b;
	int size_a_int;
	int size_b_int;
	const char* String_A;
	const char* String_B;
	int *alignmentAB, *alignmentBA;
	int i;
	int startAB, endAB, startBA, endBA;

	
	//Assigning from the py file to the c file
	if(!(PyArg_ParseTuple(args,"sksk",&String_A,&size_a,&String_B,&size_b)))
      	{
		  fprintf(stderr,"Error: Incorrect argument types\n");
        	  return NULL;
      	}
	size_a_int = ((int)size_a) ;
	size_b_int = ((int)size_b);
	alignmentAB = ACS(String_A, String_B, size_a_int, size_b_int);
	startAB = alignmentAB[0];
	endAB = alignmentAB[1];
	alignmentBA = ACS(String_B, String_A, size_b_int, size_a_int);
	startBA = alignmentBA[0];
	endBA = alignmentBA[1];
	String_A = NULL;
	String_A = NULL;
	
	
	//c to py format
    return Py_BuildValue("iiii", startAB, endAB, startBA, endBA);
}



static PyMethodDef TreeAlignmentMethods[] = {
	{"Align", KM_TreeAlign, METH_VARARGS,"Suffix tree based alignment"},
	{NULL,NULL,0,NULL}
};

PyMODINIT_FUNC
initKM_TreeAlign(void)
{
	(void) Py_InitModule("KM_TreeAlign", TreeAlignmentMethods);
}














