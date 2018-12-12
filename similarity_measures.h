/*
Includes the functions to compute the 3 kinds of similarity measures proposed in the paper -

1 Katz Similarity for Graphs,
2 Grouped Katz Similarity for Graphs,
3 Katz Index Similarity for Graphs

*/

//----------------------------- Katz similarity for Graphs -------------------------------------

typedef struct KSV {
    double** KSVmat;
    int** KSVnz;
}KSV;

// Compute Katz similarity vector using adj mat and Ord as inputs
KSV* computeKSV(int*** Tx, int* Ord, int num_vertices, double alpha){
  // allocate memory
  int i,j;
  double ** KSVmat = (double**) malloc(num_vertices*sizeof(double*));
  for(i=1;i<=num_vertices;i++) KSVmat[i-1] = (double*)malloc(2*sizeof(double));
  
  //initialize to 0
  for(i=1;i<=num_vertices;i++) {KSVmat[i-1][0] = 0;KSVmat[i-1][1] = 0;}
  
  //get maximum number of levels
  int maxLevel = 0;
  for(i=1;i<=num_vertices;i++) if(Ord[i-1]>maxLevel) maxLevel = Ord[i-1];
  
  // precompute the ordering of vertices such that they come in order of top levels - saves us a loop
  // save in the array OrdVertices
  int levelCount[maxLevel];
  for(i=1;i<=maxLevel;i++) levelCount[i-1] =0;
  for(i=1;i<=num_vertices;i++) levelCount[Ord[i-1]-1] = levelCount[Ord[i-1]-1] + 1;
  int cumLevelCount[maxLevel];
  for(i=1;i<=maxLevel;i++) cumLevelCount[i-1] =0;
  cumLevelCount[0] = 0;
  for(i=2;i<=maxLevel;i++) cumLevelCount[i-1] = cumLevelCount[i-2] + levelCount[i-2];
  
  //order vertices in the level order
  int * OrdVertices = (int*)malloc(num_vertices*sizeof(int));
  int curLevelCount[maxLevel];
  for(i=1;i<=maxLevel;i++) curLevelCount[i-1] =0;
  int o, ind;
  for(i=1;i<=num_vertices;i++){
    o = Ord[i-1];
    ind = curLevelCount[o-1] + cumLevelCount[o-1];
    OrdVertices[ind] = i;
    curLevelCount[o-1] = curLevelCount[o-1] + 1;
  }
  
  //loop over ordered vertices
  int cur_level = 1,u,v, num_parents =0, p, k, i1,i3;
  int ** parents = Tx[1];
  // also store non-zero KSV entries for each vertex to avoid the double loop
  int ** KSVnonZero = (int**)malloc(num_vertices*sizeof(int*));
  for(i=1;i<=num_vertices;i++){
    KSVnonZero[i-1] = (int*)malloc(2*sizeof(int));
    KSVnonZero[i-1][0] =0;
    KSVnonZero[i-1][1] =0; //dummy entry
  }
  
  int * tmp1;
  double * tmp2;
  int isPresent;
  for(i=1;i<=num_vertices;i++){
    if(i>cumLevelCount[cur_level]) cur_level++;
    u = OrdVertices[i-1];
    num_parents = parents[u-1][0];
    for(j=1;j<=num_parents;j++){
      p = parents[u-1][j];
      for(k=1;k<=KSVnonZero[p-1][0];k++){
        v = KSVnonZero[p-1][k];
        isPresent = 0;
        for(i1=1;i1<=KSVnonZero[u-1][0];i1++){
          if(KSVnonZero[u-1][i1] == v){
            isPresent = 1;
            // the i1th entry of KSVmat corresponds to KSV(v,u) since v is the ith node from which u is reachable
            KSVmat[u-1][i1] = KSVmat[u-1][i1] + alpha*KSVmat[p-1][k];
            break;
          }
        }
        if(isPresent ==0){
          tmp1 = (int*)malloc((KSVnonZero[u-1][0]+2)*sizeof(int));
          tmp1[0] = KSVnonZero[u-1][0] + 1;
          for(i3=1;i3<=KSVnonZero[u-1][0];i3++) tmp1[i3] = KSVnonZero[u-1][i3];
          tmp1[KSVnonZero[u-1][0] + 1] = v;
          free(KSVnonZero[u-1]);
          KSVnonZero[u-1] = tmp1;
          
          tmp2 = (double*)malloc((KSVnonZero[u-1][0]+1)*sizeof(double));
          tmp2[0] = KSVmat[u-1][0] + 1;
          for(i3=1;i3<=KSVmat[u-1][0];i3++) tmp2[i3] = KSVmat[u-1][i3];
          tmp2[KSVnonZero[u-1][0]] = alpha*KSVmat[p-1][k];
          free(KSVmat[u-1]);
          KSVmat[u-1] = tmp2;
        }
      }
      
      //add parent to nonzero KSV list
      isPresent =0;
      for(i1=1;i1<=KSVnonZero[u-1][0];i1++){
        if(KSVnonZero[u-1][i1] == p){
          isPresent = 1;
          KSVmat[u-1][i1] = KSVmat[u-1][i1] + alpha;
          break;
        }
      }
      if(isPresent ==0){
        tmp1 = (int*)malloc((KSVnonZero[u-1][0]+2)*sizeof(int));
        tmp1[0] = KSVnonZero[u-1][0] + 1;
        for(i1=1;i1<=KSVnonZero[u-1][0];i1++) tmp1[i1] = KSVnonZero[u-1][i1];
        tmp1[KSVnonZero[u-1][0] + 1] = p;
        free(KSVnonZero[u-1]);
        KSVnonZero[u-1] = tmp1;
        
        tmp2 = (double*)malloc((KSVnonZero[u-1][0]+1)*sizeof(double));
        tmp2[0] = KSVmat[u-1][0] + 1;
        for(i1=1;i1<=KSVnonZero[u-1][0]-1;i1++) tmp2[i1] = KSVmat[u-1][i1];
        tmp2[KSVnonZero[u-1][0]] = alpha;
        free(KSVmat[u-1]);
        KSVmat[u-1] = tmp2;
      }
    }
  }
  
  free(OrdVertices);
  
  KSV * K1 = (KSV*)malloc(sizeof(KSV));
  K1->KSVmat = KSVmat;
  K1->KSVnz = KSVnonZero;
  return K1;
}


// compute Katz similarity for graphs given KSV for each graph
/* Note that this function outputs the 1-norm of the difference between the Katz Index vectors.
    The similarity value is a scaled transform of this
Inputs -
1. KSVmat - Katz similarity vectors for each taxonomy
2. num_vertices - Number of vertices in the union set of two taxonomies, vertices in one taxonomy and not in the other are added as singletons. Note that the indexing of vertices is consistent across the two taxonomies. i.e A particular vertex id represents the same node in both taxonomies.
3. KSV1nonZero, KSV2nonZero - stores the non-zero entries of the KSVmat, so the indices of vertices that are reachable from one other
Output - the one norm of the difference between the Katz Similarity vectors of the two taxonomies
*/
double computeKSG(double** KSV1, int** KSV1nonZero, double **KSV2, int** KSV2nonZero, int num_vertices, double * vertexImportance){
  int i1,j1,i2, j2, i;
  double c =0;
  // initialize vertexImportance to 0
  for(i=0;i<num_vertices;i++) vertexImportance[i] =0;
  
  int* flag_arr;
  int *nz1;
  int * nz2;
  int flag =0;
  double d1=0, d2=0;
  for(i1= 1;i1<=num_vertices;i1++){
    nz1 = KSV1nonZero[i1-1];
    nz2 = KSV2nonZero[i1-1];
    
    if(nz2[0] >0) flag_arr = (int*)malloc(nz2[0]*sizeof(int)); //free later
    for(i=1;i<=nz2[0];i++) flag_arr[i-1] =0;
    for(j1 = 1;j1<=nz1[0]; j1++){
      flag =0; // is 1 if j1 has a match in the second array
      for(j2=1;j2<=nz2[0];j2++){
        if(nz1[j1] == nz2[j2]){
          
          d1 = KSV1[i1-1][j1];
          d2 = KSV2[i1-1][j2];
          
          if(d1>=d2) {c+= (d1-d2); vertexImportance[i1-1] += (d1-d2);}
          else {c+= (d2-d1); vertexImportance[i1-1] += (d2-d1);}
          flag_arr[j2-1] = 1;
          flag = 1;
          break;
        }
      }
      if(flag ==0){
        // j1 is not present in second array
        d1 = KSV1[i1-1][j1];
        c += d1;
        vertexImportance[i1-1] += (d1);
      }
    }
    
    for(j2=1;j2<=nz2[0];j2++){
      if(flag_arr[j2-1]==0){
        //this element had no match in j1
        //d2 = KSV2[nz2[j2]-1][i1-1];
        d2 = KSV2[i1-1][j2];
        c += d2;
        vertexImportance[i1-1] += (d2);
      }
    }
    //printf("Done with vertex %d\n", i1);
    if(nz2[0] >0) free(flag_arr); //free if assigned
  }
  
  return c;
  
}


//-----------------------Group Katz similarity for Graphs--------------------------

typedef struct GKSV {
    double** GKSVmat;
    int num_groups;
    int num_vertices;
}GKSV;


// Compute Grouped Katz similarity vector using adj mat and Ord as inputs
// GKSV is a dense matrix of size G*N so we would assume that the number of groups << N for it to be scalable
GKSV* computeGKSV(int*** Tx, int* Ord, int num_vertices, double alpha, int * GroupIdx, int num_groups){
  // allocate memory
  int i,j;
  double ** GKSVmat = (double**) malloc(num_groups*sizeof(double*));
  for(i=1;i<=num_groups;i++) GKSVmat[i-1] = (double*)malloc(num_vertices*sizeof(double));
  
  //initialize to 0
  for(i=1;i<=num_groups;i++) for(j=1;j<=num_vertices;j++) GKSVmat[i-1][j-1] = 0;
  
  //get maximum number of levels
  int maxLevel = 0;
  for(i=1;i<=num_vertices;i++) if(Ord[i-1]>maxLevel) maxLevel = Ord[i-1];
  
  // precompute the ordering of vertices such that they come in order of top levels - saves us a loop
  // save in the array OrdVertices
  int levelCount[maxLevel];
  for(i=1;i<=maxLevel;i++) levelCount[i-1] =0;
  for(i=1;i<=num_vertices;i++) levelCount[Ord[i-1]-1] = levelCount[Ord[i-1]-1] + 1;
  int cumLevelCount[maxLevel];
  for(i=1;i<=maxLevel;i++) cumLevelCount[i-1] =0;
  cumLevelCount[0] = 0;
  for(i=2;i<=maxLevel;i++) cumLevelCount[i-1] = cumLevelCount[i-2] + levelCount[i-2];
  
  //order vertices in the level order
  int * OrdVertices = (int*)malloc(num_vertices*sizeof(int));
  int curLevelCount[maxLevel];
  for(i=1;i<=maxLevel;i++) curLevelCount[i-1] =0;
  int o, ind;
  for(i=1;i<=num_vertices;i++){
    o = Ord[i-1];
    ind = curLevelCount[o-1] + cumLevelCount[o-1];
    OrdVertices[ind] = i;
    curLevelCount[o-1] = curLevelCount[o-1] + 1;
  }
  
  //loop over ordered vertices
  int cur_level = 1,u,v, num_parents =0, p, k, i1,i3;
  int ** parents = Tx[1];
  
  
  int pgroup=0;
  for(i=1;i<=num_vertices;i++){
    if(i>cumLevelCount[cur_level]) cur_level++;
    u = OrdVertices[i-1];
    num_parents = parents[u-1][0];
    for(j=1;j<=num_parents;j++){
      p = parents[u-1][j];
      pgroup = GroupIdx[p-1];
      if(pgroup> num_groups || pgroup<=0) printf("Error, pgroup = %d\n",pgroup);
      for(k=1;k<=num_groups;k++) GKSVmat[k-1][u-1] = GKSVmat[k-1][u-1] + alpha*GKSVmat[k-1][p-1]; // find a way to update only non-zero entries, else this is v^2
      GKSVmat[pgroup-1][u-1] = GKSVmat[pgroup-1][u-1] + alpha;
    }
  }
  
  
  // free variables
  free(OrdVertices);
  
  GKSV * GK1 = (GKSV*)malloc(sizeof(GKSV));
  GK1->GKSVmat = GKSVmat;
  GK1->num_vertices = num_vertices;
  GK1->num_groups = num_groups;
  return GK1;
}


// compute Grouped Katz similarity for graphs given GKSV for each graph
/* Note that this function outputs the 1-norm of the difference between the Grouped Katz Similarity vectors.
    The similarity value is a scaled transform of this
Inputs -
1. GKSV - Grouped Katz similarity vectors for each taxonomy
2. num_vertices - Number of vertices in the union set of two taxonomies, vertices in one taxonomy and not in the other are added as singletons. Note that the indexing of vertices is consistent across the two taxonomies. i.e A particular vertex id represents the same node in both taxonomies.
Output - the one norm of the difference between the Grouped Katz Similarity vectors of the two taxonomies
*/
double computeGKSG(GKSV* GKSV1, GKSV* GKSV2, double * vertexImportance){
  int i1,j1;
  int num_groups = GKSV1->num_groups;
  int num_vertices = GKSV1->num_vertices;
  double c =0, d1 =0, d2=0;
  
  // initialize vertexImportance to 0
  for(i1=0;i1<num_vertices;i1++) vertexImportance[i1] =0;
  
  for(i1= 1;i1<=num_groups;i1++){
    for(j1=1;j1<=num_vertices;j1++){
      d1 = GKSV1->GKSVmat[i1-1][j1-1];
      d2 = GKSV2->GKSVmat[i1-1][j1-1];
      if(d1>=d2) {c+= (d1-d2);  vertexImportance[j1-1] += (d1-d2);}
      else {c+= (d2-d1); vertexImportance[j1-1] += (d2-d1);}
    }
  }
  
  return c;
}



//------------------------------------  Katz Index Similarity for Graphs ------------------------------------
// Compute Katz index vector for a given taxonomy using adjascency matrix and topological ordering as inputs
double* computeKIV(int*** Tx, int* Ord, int num_vertices, double alpha){
  // allocate memory
  int i,j;
  double * KIVmat = (double*) malloc(num_vertices*sizeof(double));
  
  //initialize to 0
  for(j=1;j<=num_vertices;j++) KIVmat[j-1] = 0;
  
  //get maximum number of levels
  int maxLevel = 0;
  for(i=1;i<=num_vertices;i++) if(Ord[i-1]>maxLevel) maxLevel = Ord[i-1];
  
  printf("Maximum levels = %d\n", maxLevel);
  // precompute the ordering of vertices such that they come in order of top levels - saves us a loop
  // save in the array OrdVertices
  int levelCount[maxLevel];
  for(i=1;i<=maxLevel;i++) levelCount[i-1] =0;
  for(i=1;i<=num_vertices;i++) levelCount[Ord[i-1]-1] = levelCount[Ord[i-1]-1] + 1;
  int cumLevelCount[maxLevel];
  for(i=1;i<=maxLevel;i++) cumLevelCount[i-1] =0;
  cumLevelCount[0] = 0;
  for(i=2;i<=maxLevel;i++) cumLevelCount[i-1] = cumLevelCount[i-2] + levelCount[i-2];
  
  //order vertices in the level order
  int * OrdVertices = (int*)malloc(num_vertices*sizeof(int));
  int curLevelCount[maxLevel];
  for(i=1;i<=maxLevel;i++) curLevelCount[i-1] =0;
  int o, ind;
  for(i=1;i<=num_vertices;i++){
    o = Ord[i-1];
    ind = curLevelCount[o-1] + cumLevelCount[o-1];
    OrdVertices[ind] = i;
    curLevelCount[o-1] = curLevelCount[o-1] + 1;
  }
  
  //loop over ordered vertices
  int cur_level = 1,u,v, num_parents =0, p, k, i1,i3;
  int ** parents = Tx[1];
  
  for(i=1;i<=num_vertices;i++){
    if(i>cumLevelCount[cur_level]) cur_level++;
    u = OrdVertices[i-1];
    num_parents = parents[u-1][0];
    for(j=1;j<=num_parents;j++){
      p = parents[u-1][j];
      KIVmat[u-1] = KIVmat[u-1] + (KIVmat[p-1] + 1)*alpha;
    }
  }
  
  // free variables
  free(OrdVertices);
  
  return KIVmat;
}


// compute Katz index similarity for graphs given KSV for each graph
/* Note that this function outputs the 1-norm of the difference between the Katz Index vectors.
    The similarity value is a scaled transform of this
Inputs -
1. KIVmat - Katz index vectors for each taxonomy
2. num_vertices - Number of vertices in the union set of two taxonomies, vertices in one taxonomy and not in the other are added as singletons. Note that the indexing of vertices is consistent across the two taxonomies. i.e A particular vertex id represents the same node in both taxonomies.

Output - the one norm of the difference between the Katz Index vectors of the two taxonomies
*/
double computeKIG(double* KIV1mat, double* KIV2mat, int num_vertices, double* vertexImportance){
  int i1,j1;
  double c =0, d1 =0, d2=0;
  
  for(j1=1;j1<=num_vertices;j1++){
      d1 = KIV1mat[j1-1];
      d2 = KIV2mat[j1-1];
      if(d1>=d2){
        c= c + (d1-d2);
        vertexImportance[j1-1] += (d1-d2);
      }
      else{
        c = c + (d2-d1);
        vertexImportance[j1-1] += (d2-d1);
      }
  }
  
  return c;
}