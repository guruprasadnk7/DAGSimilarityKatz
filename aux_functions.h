/* Auxiliary functions -
1. Reading the edge list for a DAG
2. Computing a topological ordering for the DAG
3. Free memory allocared to a Taxonomy
*/


//------------------- Read the adj matrix from the text file ----------------------------
// each line of the file should have the format 'u v'
// This shows that there is an edge from vertex id u to vertex id v in the graph
int*** readEdgeListFile(const char * file_name, int num_vertices){
  FILE *fp;
  
  int i,j,u,v;
  
  
  // Count the number of lines in the file
  fp = fopen(file_name,"r");
  if(fp == NULL) printf("Error opening file\n");
  char str[100];
  int line_count = 0;
  while(fgets (str, 100, fp)!=NULL) line_count++;
  fclose(fp);
  //printf("Line count = %d\n", line_count);
  
  // count parents for each vertex
  int * parent_count = (int*)malloc(num_vertices*sizeof(int));
  for(i=1;i<=num_vertices;i++) parent_count[i-1] =0;
  
  // count children for each vertex
  int * children_count = (int*)malloc(num_vertices*sizeof(int));
  for(i=1;i<=num_vertices;i++) children_count[i-1] =0;
  
  // read into adjascency mat
  fp = fopen(file_name,"r");
  for(i=0;i<line_count;i++){
    fscanf(fp, "%d %d", &u, &v);
    parent_count[v-1] = parent_count[v-1] + 1;
    children_count[u-1] = children_count[u-1] + 1;
  }
  fclose(fp);
  
  // store parents and children for each vertex
  int ** parents = (int**)malloc(num_vertices*sizeof(int*));
  for(i=1;i<=num_vertices;i++) {
    parents[i-1] = (int*)malloc((parent_count[i-1]+1)*sizeof(int));
    parents[i-1][0] = parent_count[i-1];
  }
  
  int ** children = (int**)malloc(num_vertices*sizeof(int*));
  for(i=1;i<=num_vertices;i++) {
    children[i-1] = (int*)malloc((children_count[i-1]+1)*sizeof(int));
    children[i-1][0] = children_count[i-1];
  }
  
  int *curCountParents = (int*)malloc(num_vertices*sizeof(int));
  int *curCountChildren = (int*)malloc(num_vertices*sizeof(int));
  for(i=1;i<=num_vertices;i++){
    curCountParents[i-1] =0;
    curCountChildren[i-1] =0;
  }
  
  // read into parents, children mat
  fp = fopen(file_name,"r");
  int ind;
  for(i=0;i<line_count;i++){
    fscanf(fp, "%d %d", &u, &v);
    ind = curCountParents[v-1] + 1;
    parents[v-1][ind] = u;
    curCountParents[v-1] = curCountParents[v-1] + 1;
    if(curCountParents[v-1] > parent_count[v-1]) printf("Error\n");
    
    
    ind = curCountChildren[u-1] + 1;
    children[u-1][ind] = v;
    curCountChildren[u-1] = curCountChildren[u-1] + 1;
    if(curCountChildren[u-1] > children_count[u-1]) printf("Error\n");
    
  }
  fclose(fp);
  
  // clear temporary variables
  free(curCountParents);
  free(curCountChildren);
  free(parent_count);
  free(children_count);
  
  int *** Tx = (int***)malloc(3*sizeof(int**));
  Tx[0] = NULL;
  Tx[1] = parents;
  Tx[2] = children;
  return Tx;
}

//------------------------- Free the adj matrix, parents, children arrays --------------------
void freeTx(int*** Tx, int num_vertices){
  int i;
  //for(i=1;i<=num_vertices;i++) free(Tx[0][i-1]);
  //free(Tx[0]);
  
  for(i=1;i<=num_vertices;i++) free(Tx[1][i-1]);
  free(Tx[1]);
  
  for(i=1;i<=num_vertices;i++) free(Tx[2][i-1]);
  free(Tx[2]);
  
  free(Tx);
  
  return;
}


//------------------------------- Create topological ordering with adj matrix as input -------------------------
int * getTopologicalOrdering(int*** Tx, int num_vertices){
	//Count the number of edges incident on each vertex
	int i,j,k, rem_set, flag =0;
	int ** parents = Tx[1];
	int ** children = Tx[2];
	int * inDegree = (int*)malloc(num_vertices*sizeof(int));
	for(i = 1;i<=num_vertices;i++) inDegree[i-1] = parents[i-1][0];
	
	int * Ord = (int*)malloc(num_vertices*sizeof(int));
	for(i=1;i<=num_vertices;i++) Ord[i-1] = 0;
	int cur_level = 1; // level starts from 1, root node is at level 1
	rem_set = num_vertices; //counts the number of vertices left after each iteration
	while(rem_set>0){
		for(i=1;i<=num_vertices;i++){
			if(inDegree[i-1] == 0){
				//add i to cur level
				Ord[i-1] = cur_level;
				inDegree[i-1] = -1; //So that it is not picked on next iteration
				rem_set = rem_set -1;
			}
		}
		// Reduce indegrees separately otherwise next iter vertices may get into this level
		for(i=1;i<=num_vertices;i++){
			if(Ord[i-1] == cur_level){
				//this vertex was added in this iteration
				//remove all edges starting from i, i.e reduce destination vertex indegrees by 1
				for(j=1;j<=children[i-1][0];j++) inDegree[children[i-1][j]-1]  = inDegree[children[i-1][j]-1] -1;
			}
		}
		cur_level = cur_level + 1;
	}
	
	// free vars
	free(inDegree);
	return Ord;
}
