#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <strings.h>
#include <stdlib.h>
#include <Python.h>
//#include "idc.h"



/**********************************************************************
*                                                                      *
*  File: kamivq.c                                                      *
*  modified from trvqsplit.c on 2/2/2004                               *
*  modified from amivq for Keith 2/17/2014                             *
*  Function:  trains a codebook using the splitting approach of Linde, *
*             Buzo, and Gray                                           *
*  Author (aka person to blame) : K. Sayood                            *
*  Last mod: 5/12/95                                                   *
*  Usage:  see usage().  For more details see trvqsplit.doc or the man *
*          page                                                        *
***********************************************************************/



#define maxit 100

void usage(void);

cluster(char *tfpT, char *ofpT, int final_codebook_size, int  dimension, int ts_size)
{
  FILE *tfp, *ofp;
  char codefile[50];
  float **codebook;
  float **training_set;
  //int  dimension, index, codebook_size, ts_size, c;
  int index, codebook_size, c;
  //int i, j,  iteration, tsvec, final_codebook_size, *count, **membership;
  int i, j,  iteration, tsvec, *count, **membership;
  float measure, total_distortion, total_distortion_old;
  float **new_codebook, *eps, *input, *cluster_distortion;
  float threshold, distortion;
  extern int optind;
  extern char *optarg;

  // final_codebook_size = -1;
  // dimension = -1;
  // ts_size = -1;
   
  if((tfp = fopen(tfpT,"r"))==NULL)
  {
   fprintf(stderr,"File open failed for the training set file %s\n",tfpT);
   usage();
   exit(1);
  }
ofp = fopen(ofpT,"w");
tfp = fopen(tfpT,"r");
/* Obtain the file containing the training set */
/*
  if((tfp = fopen(argv[1],"r"))==NULL)
  {
   fprintf(stderr,"File open failed for the training set file %s\n",argv[1]);
   usage();
   exit(1);
  }
  optind++;

  ofp = stdout;
  strcpy(codefile,"standard out");
 
  while((c=getopt(argc,argv,"o:b:d:t:h")) != EOF)
  {
   switch (c){
   case 'o':
         strcpy(codefile,optarg);
         ofp = fopen(optarg,"w");
         break;
    case 'b':
     sscanf(optarg,"%d", &final_codebook_size);
     break;
    case 'd':
     sscanf(optarg,"%d", &dimension);
     break;
    case 't':
     sscanf(optarg,"%d", &ts_size);
     break;
   case 'h':
         usage();
         exit(1);
             }
   }

   if(final_codebook_size < 0)
   {
    fprintf(stderr,"Enter codebook size:\n");
    scanf("%d",&final_codebook_size);
   }

   if(dimension < 0)
   {
    fprintf(stderr,"Enter vector dimension:\n");
    scanf("%d",&dimension);
   }


   if(ts_size < 0)
   {
    fprintf(stderr,"Enter size of training set:\n");
    scanf("%d",&ts_size);
   }

   fprintf(stderr,"\n\n\t\t Vector Quantizer Design Using the LBG algorithm \n\n");

   fprintf(stderr,"\t This program can be used to train a vector quantizer using the \n");
   fprintf(stderr,"\t LBG (GLA) algorithm.  The codebook is initialized using \n");
   fprintf(stderr,"\t the splitting method. If the user wishes to use an initial\n");
   fprintf(stderr,"\t codebook for initialization s/he should use the program trainvq.\n");  
   fprintf(stderr,"\t This program reads in the training set as ints.  The user\n");
   fprintf(stderr,"\t  can modify this by changing the read statement. \n");
*/

   //printf("Codebook size = %d, Training set = %d, dimension = %d\n", final_codebook_size, ts_size, dimension);
   codebook_size = final_codebook_size;
  
/*  Allocate space for codebook, training set, input and associated work
    arrays */

    codebook = (float **) calloc(codebook_size,sizeof(float *));
    membership = (int **) calloc(codebook_size,sizeof(int *));
    cluster_distortion = (float *) calloc(codebook_size, sizeof(float));

    //fprintf(stderr,"First Allocation for codebook\n");
    for(i=0; i< codebook_size; i++)
    {
      codebook[i] = (float *) calloc(dimension,sizeof(float));
      membership[i] = (int *) calloc(ts_size,sizeof(int));
    }
    //fprintf(stderr,"Allocation for codebook\n");

    training_set = (float **) calloc(ts_size,sizeof(float *));
    for(i=0; i< ts_size; i++)
      training_set[i] = (float *) calloc(dimension,sizeof(float));
    //fprintf(stderr,"Allocation for training set\n");


    input = (float *) calloc(dimension,sizeof(float));

    new_codebook = (float **) calloc(codebook_size,sizeof(float *));
    for(i=0; i< codebook_size; i++)
      new_codebook[i] = (float *) calloc(dimension,sizeof(float));

    count = (int *) calloc(codebook_size,sizeof(int));

    eps = (float *) calloc(dimension,sizeof(float));
//Calloc'ed: count, eps, new_codebook, training_set, codebook, membership, cluster_distortion, input

    //fprintf(stderr,"Allocations complete\n");
 
   for(i=0; i< ts_size; i++)
   {
    for(j=0; j < dimension; j++)
    {
     fscanf(tfp,"%f",&training_set[i][j]);
//     fprintf(stderr,"%f ",training_set[i][j]);
    }
//    fprintf(stderr,"\n");
   }

/* Set up perturbation vector */
/* For now use fixed perturbation - this can be changed to a random perturbation
   vector    */

   for(i=0; i<dimension; i++)
     eps[i] = .0001;

/* Obtain one level VQ */

   for(j=0; j<dimension; j++)
    codebook[0][j] = 0;

   for(i=0; i<ts_size; i++)
    for(j=0; j<dimension; j++)
     codebook[0][j] += training_set[i][j];

   for(j=0; j<dimension; j++)
   {
     codebook[0][j] = ((float) codebook[0][j]/(float) ts_size );
     //fprintf(stderr,"%f  ",codebook[0][j]);
   }
   //fprintf(stderr,"\n");

   codebook_size = 1;
   while(codebook_size < final_codebook_size)
   {
    codebook_size = 2*codebook_size;
    
    for(i=0; i<codebook_size/2; i++)
     for(j=0; j<dimension; j++)
     {
      new_codebook[2*i][j] = codebook[i][j] + eps[j];
      new_codebook[2*i+1][j] = codebook[i][j];
     }
    for(i=0; i<codebook_size; i++)
     for(j=0; j<dimension; j++)
      codebook[i][j] = new_codebook[i][j];

    measure = 1000.;
    threshold = 0.0000000001;
    iteration = 0;
    while(measure > threshold && iteration < maxit)
    {
     if(iteration > 0)
     {
      total_distortion_old = total_distortion;
      for(i=0;i<codebook_size;i++)
      {
       for(j=0; j<dimension; j++)
        if(count[i] > 0)
         codebook[i][j] =((float) new_codebook[i][j]/(float) count[i] );
      }
     }


/* do something about zero count */



     for(i=0;i<codebook_size;i++)
      for(j=0; j<dimension; j++)
       new_codebook[i][j] = 0;
     for(i=0;i<codebook_size;i++)
       count[i] = 0;
     total_distortion = 0;
     for(tsvec=0; tsvec<ts_size; tsvec++)
     {
      for(i=0; i<dimension; i++)
       input[i] = training_set[tsvec][i];
      
      index = fvqe(input,codebook,codebook_size,dimension,&distortion);
 
      total_distortion += distortion;
      count[index]++;
      for(j=0; j<dimension; j++)
       new_codebook[index][j] += input[j];
      }
      total_distortion = total_distortion/(float) ts_size;
      if(iteration > 0)
       measure = (total_distortion_old - total_distortion)/total_distortion_old;
      else
       measure = total_distortion;
      iteration++;
     }
     //fprintf(stderr,"# of iterations: %d\n",iteration);
 
   }

/* Putting training set vectors into final clusters */

     for(i=0; i<final_codebook_size; i++)
     {
       cluster_distortion[i] = 0;
       count[i] = 0;
     }

     for(tsvec=0; tsvec<ts_size; tsvec++)
     {
      for(i=0; i<dimension; i++)
       input[i] = training_set[tsvec][i];

      index = fvqe(input,codebook,codebook_size,dimension,&distortion);
      cluster_distortion[index] += distortion;
      membership[index][count[index]] = tsvec;
      count[index]++;
     }

   for(i=0; i<final_codebook_size; i++)
   {
    if(count[i] > 0)
    cluster_distortion[i] = cluster_distortion[i]/count[i];
   }

   for(i=0; i<final_codebook_size; i++)
   {
    fprintf(ofp,"Codebook centroid: \n");
    for(j=0; j<dimension; j++)
     fprintf(ofp,"%f  ",codebook[i][j]);
    fprintf(ofp,"\n");
    fprintf(ofp,"count[%d]=%d\n",i,count[i]);
    fprintf(ofp,"Cluster membership: \n");
    for(j=0; j<count[i]; j++)
      fprintf(ofp,"%d ",membership[i][j]);
    fprintf(ofp,"\n");
    fprintf(ofp,"Cluster distortion: %f\n",cluster_distortion[i]);
    fprintf(ofp,"\n"); 
   }
   //printf(" ********  Final distortion :  %f  ********\n",total_distortion);
   fclose(ofp);
   fclose(tfp);
   
//Calloc'ed: count, eps, new_codebook, training_set, codebook, membership, cluster_distortion, input);
   free(count);
   free(eps);
   free(new_codebook);
   free(training_set);
   free(codebook);
   free(membership);
   free(cluster_distortion);
   free(input);
   return;
 }
 
 void usage(void)
{
  fprintf(stderr,"Usage:\n");
 fprintf(stderr,"amivq file1 [-o outfile] [-b codebook size] [-d dimension] [-t training set size] [-h]\n\n");
  fprintf(stderr,"file1  : File containing training vectors \n");

  fprintf(stderr,"outfile: File containing codebook vectors.  If a filename is\n");
  fprintf(stderr,"not specified, the values are written to standard out\n \n");
  fprintf(stderr," The elements of the vectors are read in as int.  If this is \n");
  fprintf(stderr," not desirable the user should modify the read statement.\n\n");
  fprintf(stderr,"\t codebook size    :  The number of vectors in the codebook.\n");
  fprintf(stderr,"\t dimension        :  The dimension of the vectors (total number of \n");
  fprintf(stderr,"\t\t\t\t elements in each vector).  \n");
  fprintf(stderr,"\t training set size: The number of vectors in the training set. \n");
}  

static PyObject *
kvqsplitForClust(PyObject *self, PyObject *args)
{
	int i;
	char *tfpP, *ofpP;
	int *final_codebook_sizeP, *dimensionP, *ts_sizeP;
	char *tfp, *ofp;
	int final_codebook_size, dimension, ts_size;
	float total_distortion;

	
	//Assigning from the py file to the c file
	if(!(PyArg_ParseTuple(args,"ssiii",&tfpP,&ofpP,&final_codebook_sizeP,&dimensionP,&ts_sizeP)))
      	{
		  fprintf(stderr,"Error: Incorrect argument types\n");
        	  return NULL;
      	}

	tfp = tfpP;
	ofp = ofpP;
	final_codebook_size = (int)final_codebook_sizeP;
	dimension = (int)dimensionP;
	ts_size = (int)ts_sizeP;
	
	cluster(tfp, ofp, final_codebook_size, dimension, ts_size);
	//printf("MADE IT\n");
	//c to py format
    return Py_BuildValue("f", total_distortion);
}



static PyMethodDef ClusterMethods[] = {
	{"cluster", kvqsplitForClust, METH_VARARGS,"Linde, Buzo, and Gray Clustering Algorithm by K. Sayood"},
	{NULL,NULL,0,NULL}
};

PyMODINIT_FUNC
initkvqsplitForClust(void)
{
	(void) Py_InitModule("kvqsplitForClust", ClusterMethods);
}




























