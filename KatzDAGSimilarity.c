#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "aux_functions.h"
#include "similarity_measures.h"

void main(int argc, char **argv){
	//------------------------- Read input and create variables------------------------
	if(argc <9){
	  printf("Not enough input arguments\\");
	  return;
	}
	char * file_name1 = argv[1];
	char * file_name2 = argv[2];
	double alpha = atof(argv[3]);
	double gamma = atof(argv[4]);
	int num_vertices = atoi(argv[5]);
	int method_id = atoi(argv[6]); /* set to
	1 for Katz Similarity for Graphs,
	2 for Grouped Katz Similarity for Graphs,
	3 for Katz Index Similarity for Graphs */
	char * output_file_name = argv[7];
	
	// declare other variables
  int num_groups;
	KSV* K1, *K2;
	GKSV* GKS1, *GKS2;
	double* KIV1mat, *KIV2mat;
	int * GroupIdx;
	
	if(method_id ==2){
	  // Grouped Katz Similarity
	  if(argc <10){
	    printf("Not enough input arguments\\");
	    return;
	  }
	  num_groups = atoi(argv[8]);
	  srand(time(NULL));   // should only be called once
	  GroupIdx = (int*)malloc(num_vertices*sizeof(int));
    int i, r;
    for(i=1;i<=num_vertices;i++) GroupIdx[i-1] = rand()%num_groups + 1; ; // randomly assign a positive integer from 1 to num_groups
	}
	
	// Read in first taxonomy and compute topological ordering
	int*** Tx1 = readEdgeListFile(file_name1, num_vertices);
  int * Ord1 = getTopologicalOrdering(Tx1, num_vertices);
  
  // Read in Second taxonomy and compute topological ordering
  int*** Tx2 = readEdgeListFile(file_name2, num_vertices);
  int * Ord2 = getTopologicalOrdering(Tx2, num_vertices);
  
  // Create array for vertex imporatnce
  double * vertexImportance = (double*)malloc(num_vertices*sizeof(double));
  
  
  //-------------------- compute the similarity ---------------------------------
  double graph_similarity = 0;
  double vector_difference_norm = 0;
  if(method_id==1){
    // compute Katz Similarity Vectors
    K1 = computeKSV( Tx1,  Ord1,  num_vertices, alpha);
    K2 = computeKSV( Tx2,  Ord2,  num_vertices, alpha);
    
    // Compute the one-norm of the difference between the two vectors
    vector_difference_norm = computeKSG( K1->KSVmat, K1->KSVnz, K2->KSVmat, K2->KSVnz,num_vertices, vertexImportance);
  }
  
  if(method_id ==2){
    // Compute the Group Katz Similarity Vectors
    GKS1 = computeGKSV( Tx1,  Ord1,  num_vertices, alpha, GroupIdx, num_groups);
    GKS2 = computeGKSV( Tx2,  Ord2,  num_vertices, alpha, GroupIdx, num_groups);
    
    // Compute the one-norm of the difference between the two vectors
    vector_difference_norm = computeGKSG( GKS1, GKS2, vertexImportance);
  }
  
  if(method_id ==3){
    KIV1mat = computeKIV( Tx1,  Ord1,  num_vertices, alpha);
    KIV2mat = computeKIV( Tx2,  Ord2,  num_vertices, alpha);
    vector_difference_norm = computeKIG( KIV1mat, KIV2mat, num_vertices, vertexImportance);
  }
  
  graph_similarity = (2.0*exp(-gamma*vector_difference_norm))/(1.0+exp(-gamma*vector_difference_norm));
  
  //---------------------- Write output to files--------------------------------
  // write the similarity value and the one norm of the vector difference to a file
  FILE * fp;
  fp = fopen(output_file_name,"w");
  fprintf(fp,"%lf %lf", graph_similarity, vector_difference_norm);
  fclose(fp);
  
  // write vertex importance array to a file
  char * output_file_name_importance = argv[argc-1]; //last argument is the name of the file to write vertex importance to
	fp = fopen(output_file_name_importance,"w");
	int i =0;
	for(i=0; i<num_vertices;i++) fprintf(fp, "%30.4lf\n", vertexImportance[i]);
	fclose(fp);
	
	
	//------------------------ free variables ------------------------------------
	freeTx(Tx1,num_vertices);
  freeTx(Tx2,num_vertices);
	free(Ord1);
  free(Ord2);
  free(vertexImportance);
  
  if(method_id==1){
    for(i=1;i<=num_vertices;i++) free(K1->KSVmat[i-1]);
    free(K1->KSVmat);
    for(i=1;i<=num_vertices;i++) free(K1->KSVnz[i-1]);
    free(K1->KSVnz);
    free(K1);
  
    for(i=1;i<=num_vertices;i++) free(K2->KSVmat[i-1]);
    free(K2->KSVmat);
    for(i=1;i<=num_vertices;i++) free(K2->KSVnz[i-1]);
    free(K2->KSVnz);
    free(K2);
  }
  
  if(method_id==2){
    for(i=1;i<=num_groups;i++) free(GKS1->GKSVmat[i-1]);
    free(GKS1->GKSVmat);
    free(GKS1);
    for(i=1;i<=num_groups;i++) free(GKS2->GKSVmat[i-1]);
    free(GKS2->GKSVmat);
    free(GKS2);
    free(GroupIdx);
  }
  
  if(method_id==3){
    free(KIV1mat);
    free(KIV2mat);
  }
  
  return;
}
  
  