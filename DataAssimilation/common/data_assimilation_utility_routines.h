#define PrintLine(LINE_STR)  PetscPrintf(PETSC_COMM_WORLD, "  " LINE_STR "\n");


Mat petsc_create_matrix(int type, int size_1, int size_2, MPI_Comm comm)
{
  Mat new_mat;
  if(type==1){
    MatCreateMPIDense(comm, PETSC_DECIDE, PETSC_DECIDE, size_1,size_2, PETSC_NULL, &new_mat); 
    MatAssemblyBegin(new_mat  ,MAT_FINAL_ASSEMBLY ); 
    MatAssemblyEnd  (new_mat  ,MAT_FINAL_ASSEMBLY ); 
    MatSetFromOptions(new_mat);
  }
  return new_mat;
}
Vec petsc_create_vector(int type, int length, MPI_Comm comm)
{
  Vec new_vec;
  if(type==1){
    VecCreate(comm, &new_vec);
    VecSetSizes(new_vec, PETSC_DECIDE,length);
    VecSetFromOptions(new_vec);
  }
  return new_vec;
}


//  ierr=VecCreateMPI(comm, PETSC_DECIDE, m_A, &yy);CHKERRQ(ierr);
//PetscErrorCode PETSCMAT_DLLEXPORT MatGetColumnVector(Mat A,Vec yy,PetscInt col);
//The vector yy must have the same parallel row layout as the matrix A. 
//Z(rowShift_Z: ,columnNum) = x(low_x_set:high_x_set-1)
//
//Set a column of Z to be x
//Afterwards, MatAssemblyBegin(*Z,MAT_FINAL_ASSEMBLY) and MatAssemblyEnd(*Z,MAT_FINAL_ASSEMBLY) should be called
//
//PetscErrorCode MatSetColumnVec(Mat* Z, Vec& x,  PetscInt* indices_v, PetscInt columnNum)
//attention: call the function below afterwards!!!
//MatAssemblyBegin(*Z,MAT_FINAL_ASSEMBLY );
//MatAssemblyEnd  (*Z,MAT_FINAL_ASSEMBLY );
//rowShift default should be 0
//low_x_set default should be 0
//high_x_set default should be VecGetSize(x, &high_x_set)
//indices_v: has allocated space outside, with size of (high_x_set-low_x_set)
//
//TODO: check the parallel version of this!!
PetscErrorCode MatSetColumnVec(Mat* Z, Vec* x,  PetscInt* indices_v, PetscInt columnNum, PetscInt rowShift_Z, PetscInt low_x_set, PetscInt high_x_set)
{
  PetscInt n_x_local,low_x,high_x, portion_a, portion_b, j, m_Z, n_Z, i;
  PetscErrorCode    ierr;
  PetscScalar* v;
    
  ierr=VecGetLocalSize(*x, &n_x_local);CHKERRQ(ierr); //get the result back 
  ierr=VecGetArray(*x, &v);CHKERRQ(ierr); 
  ierr=VecGetOwnershipRange(*x, &low_x, &high_x);CHKERRQ(ierr); 
  portion_a = low_x;
  portion_b = high_x;
  if(low_x<low_x_set){
    low_x = low_x_set;
  }
  if(high_x>high_x_set){
    high_x = high_x_set;
  }
  MatGetSize(*Z, &m_Z, &n_Z);
//  for (j=low_x; j<high_x; j++) indices_v[j-low_x]=j-low_x_set + rowShift_Z;
  for (j=0; j<high_x-low_x; j++) indices_v[j]=j+rowShift_Z+(low_x-low_x_set); //+shift of vector +shift of row in matrix
  assert(n_x_local>=high_x-low_x);
  n_x_local = high_x-low_x;
  if(n_x_local<0){
    //assert(0); //probably becaue there is no intersection between the elements of the vector this rank holds and the portor we are setting, but need to be confirmed if this happens!
    n_x_local = 0;
  }    
  
  assert(high_x-low_x_set + rowShift_Z<=m_Z);
  assert(columnNum<=n_Z);
#ifdef PETSC_OPERATIONS_DEBUG
  for(i=0;i<n_x_local; i++){
        PetscPrintf(PETSC_COMM_WORLD,"MatSetColumnVec: set MatSetValues[%d,%d]=%f\n", indices_v[i], columnNum, *(v+(low_x-portion_a)+i));
  }
#endif    
  ierr=MatSetValues(*Z, n_x_local, indices_v, 1, &columnNum, v+(low_x-portion_a), INSERT_VALUES );CHKERRQ(ierr); //set Global value
  ierr=VecRestoreArray(*x, &v);CHKERRQ(ierr); //only can get the local values
 
//  MatAssemblyBegin(*Z,MAT_FINAL_ASSEMBLY );
//  MatAssemblyEnd  (*Z,MAT_FINAL_ASSEMBLY );
  
  VecAssemblyBegin(*x);
  VecAssemblyEnd  (*x);
  
  return ierr;
}
// 
//Z(row_Z:,col_Z:) = B(row_B+row_span_B-1:,col_B:col_B+col_span_B-1)
PetscErrorCode MatSetMatBlock(MPI_Comm comm, Mat* Z, Mat B, PetscInt row_Z, PetscInt col_Z, PetscInt row_B, PetscInt col_B, PetscInt row_span_B, PetscInt col_span_B)
{

  PetscInt j, m_B, n_B;
  PetscErrorCode    ierr;
  Vec b;
  PetscInt* indices_v;
#ifdef PETSC_OPERATIONS_DEBUG
  PetscViewer pviewer = PETSC_VIEWER_STDOUT_(comm);
#endif

  //allocate space for later usage
  MatGetSize(B, &m_B, &n_B);
  ierr=VecCreateMPI(comm, PETSC_DECIDE, m_B, &b);CHKERRQ(ierr);
  ierr=PetscMalloc(sizeof(PetscInt)*col_span_B, &indices_v);CHKERRQ(ierr); 
#ifdef PETSC_OPERATIONS_DEBUG
  MatView(B,pviewer);
#endif
  //set each column
  for(j=0; j<col_span_B; j++){
    ierr=MatGetColumnVector(B,b,col_B+j);CHKERRQ(ierr); 
    ierr=MatSetColumnVec(Z, &b, indices_v, col_Z+j,row_Z, row_B, row_B+row_span_B);CHKERRQ(ierr); 
#ifdef PETSC_OPERATIONS_DEBUG
  VecView(b,pviewer);
#endif
  }  
  
  //begin assembly
  MatAssemblyBegin(*Z, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*Z, MAT_FINAL_ASSEMBLY);

#ifdef PETSC_OPERATIONS_DEBUG
  MatView(*Z,pviewer);
#endif
  
  return ierr;
}

// 
//z(row_Z:) = b(row_B+row_span_B-1)
PetscErrorCode VecSetVecBlock(MPI_Comm comm, Vec* z, Vec b, PetscInt row_Z, PetscInt row_B, PetscInt row_span_B)
{

  PetscErrorCode    ierr;
  PetscInt *indices_b, *indices_z; 
  PetscReal* values_v; 

  PetscMalloc(sizeof(PetscInt)  * row_span_B, &indices_b);
  PetscMalloc(sizeof(PetscInt)  * row_span_B, &indices_z);
  PetscMalloc(sizeof(PetscReal) * row_span_B, &values_v );
  
  for(int i = row_B; i<row_B+row_span_B; i++) indices_b[i-row_B]=i;
  for(int i = row_Z; i<row_Z+row_span_B; i++) indices_z[i-row_Z]=i;

#ifdef PETSC_OPERATIONS_DEBUG
  PetscViewer pviewer = PETSC_VIEWER_STDOUT_(comm);
#endif

#ifdef PETSC_OPERATIONS_DEBUG
  PrintLine("input vector used by VecSetVecBlock:");
  VecView(b,pviewer);
#endif
  ierr=VecGetValues( b, row_span_B, indices_b, values_v);CHKERRQ(ierr); 
  ierr=VecSetValues(*z, row_span_B, indices_z, values_v, INSERT_VALUES);CHKERRQ(ierr); 
  
  VecAssemblyBegin(*z);
  VecAssemblyEnd(*z);

#ifdef PETSC_OPERATIONS_DEBUG
  PrintLine("ouput vector modifed by VecSetVecBlock:");
  VecView(*z,pviewer);
#endif
  
  return ierr;
}

//Z = X*Y, comm should be the same as that in X,Y, *Z should not be the same as X
PetscErrorCode MatMultilXY(Mat X, Mat Y, Mat* Z, MPI_Comm comm)
{
  PetscInt m_X, n_X, m_Y,n_Y, i;
  Vec y, x;
  PetscInt* indices_v;
  PetscErrorCode    ierr;

  if(*Z==X || *Z==Y){PetscPrintf(comm, "MatMultilXY: output matrix should not be the same as input matrix\n"); assert(0);}
  
  MatGetSize(X, &m_X, &n_X);
  MatGetSize(Y, &m_Y, &n_Y);

  assert(n_X == m_Y);
  
  ierr=VecCreateMPI(comm, PETSC_DECIDE, m_Y, &y);CHKERRQ(ierr);
  ierr=VecCreateMPI(comm, PETSC_DECIDE, m_X, &x);CHKERRQ(ierr);
  
  PetscMalloc(sizeof(PetscInt)*m_X, (void**) &indices_v);
  
  for (i=0; i<m_X; i++) indices_v[i]=i;
  
  for(i=0; i<n_Y; i++){ //for each column of Y
    ierr=MatGetColumnVector(Y, y, i); CHKERRQ(ierr);
    ierr=MatMult(X,y, x);CHKERRQ(ierr); //X*y->x
    ierr=MatSetColumnVec(Z, &x, indices_v, i, 0, 0, m_X); CHKERRQ(ierr);
  }
  MatAssemblyBegin(*Z, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*Z, MAT_FINAL_ASSEMBLY);
  
//  PetscFree(v);
  PetscFree(indices_v);
  VecDestroy(x);
  VecDestroy(y);
  
  return ierr;
}
//Z = X*Y', comm should be the same as that in X,Y
PetscErrorCode MatMultilXYT(Mat X, Mat Y, Mat* Z, MPI_Comm comm)
{
  Mat Ytrans;  
  PetscErrorCode    ierr;
  ierr=MatTranspose(Y, MAT_INITIAL_MATRIX,&Ytrans);
  //ierr=MatCreateTranspose(Y,&Ytrans);
  ierr=MatMultilXY(X,Ytrans,Z, comm); CHKERRQ(ierr);
  ierr=MatDestroy(Ytrans); CHKERRQ(ierr);
  return ierr;  
}
//Z = X'*Y, comm should be the same as that in X,Y
PetscErrorCode MatMultilXTY(Mat X, Mat Y, Mat* Z, MPI_Comm comm)
{
  Mat Xtrans;  
  PetscErrorCode    ierr;
  ierr=MatTranspose(X, MAT_INITIAL_MATRIX,&Xtrans);
  ierr=MatMultilXY(Xtrans,Y,Z, comm); CHKERRQ(ierr);
  ierr=MatDestroy(Xtrans); CHKERRQ(ierr);
  return ierr;  
}
//Z = X+[y y ... y y], comm should be the same as that in X,y
//Z can not be X itself!
PetscErrorCode MatAddVec(Mat X, Vec y, Mat* Z, MPI_Comm comm)
{
  PetscInt m_X, n_X, m_y, n_Z, m_Z, i;
  Vec x;
  PetscInt* indices_v;
  PetscErrorCode    ierr;
  
  MatGetSize(X, &m_X, &n_X);
  MatGetSize(*Z, &m_Z, &n_Z);
  VecGetSize(y, &m_y);

  assert(m_X == m_y); //row-size equal
  assert(m_X == m_Z); //X,Z
  assert(n_X == n_Z); //X,Z
  
  ierr=VecCreateMPI(comm, PETSC_DECIDE, m_X, &x);CHKERRQ(ierr);
  
  ierr=PetscMalloc(sizeof(PetscInt)*m_X, (void**) &indices_v);CHKERRQ(ierr);
    
  for(i=0; i<n_X; i++){ //for each column of X
    //only difference with MatMultilXY!!!
    ierr=MatGetColumnVector(X, x, i); CHKERRQ(ierr);
    ierr=VecAXPY(x,1.0, y);CHKERRQ(ierr); //x+y->x

    ierr=MatSetColumnVec(Z, &x, indices_v, i, 0, 0, m_X); CHKERRQ(ierr);
  }
  MatAssemblyBegin(*Z, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*Z, MAT_FINAL_ASSEMBLY);
  
  PetscFree(indices_v);
  VecDestroy(x);
  
  return ierr;
}


//Z = X+[y y ... y y], comm should be the same as that in X,y
//Z can not be X itself!
//a_X default is 1
PetscErrorCode MatAddVecToVecs(Mat X, Vec y, Vec* zs, PetscInt n_zs, MPI_Comm comm, PetscScalar a_X)
{
  PetscInt m_X, n_X, m_y, m_zs, i;
  PetscErrorCode    ierr;
    
  assert(n_zs > 0 && zs[0]!=PETSC_NULL);
  
  MatGetSize(X, &m_X, &n_X);
  VecGetSize(zs[0], &m_zs);
  VecGetSize(y, &m_y);

  assert(m_X == m_y); //row-size equal
  assert(m_X == m_zs); //X,Z
  assert(n_X == n_zs); //X,Z
        
  for(i=0; i<n_X; i++){ //for each column of X
    //only difference with MatAddVec!!! store it to an array
    ierr=MatGetColumnVector(X, zs[i], i); CHKERRQ(ierr);
    assert(zs[i]!=PETSC_NULL);
    ierr=VecAYPX(zs[i],a_X, y);CHKERRQ(ierr); //a_X*zs[i]+y->zs[i]
  }
  
  return ierr;
}

//A = [x[i]], i=1...n; comm should be the same as that in x,x[i] is a vec type and A points to Mat type 
PetscErrorCode CreateMatFromVecs(Vec x[], PetscInt n, Mat* A, MPI_Comm comm)
{
  PetscInt i=0, m_A=0, n_A=0, n_x=0;
  PetscInt* indices_v;
  PetscErrorCode    ierr;
    
  assert(n>0);
  VecGetSize(x[0], &n_x);
  
  if(n>1){
    VecGetSize(x[1], &i);
    assert(i==n_x);
  }
  
  m_A = n_x; //row-size of A
  n_A = n;   //column-size of A
  
  ierr=MatCreateMPIDense(comm, PETSC_DECIDE, PETSC_DECIDE, m_A ,n_A, PETSC_NULL, A); CHKERRQ(ierr);
  
  ierr=PetscMalloc(sizeof(PetscInt)*n_x, (void**) &indices_v);CHKERRQ(ierr);

  for(i=0; i<n; i++)
    MatSetColumnVec(A, &x[i], indices_v, i, 0, 0, n_x);
  
  MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);
 
  PetscFree(indices_v); 
  return ierr;
}

//U*SIGMA*V' = A
//put PETSC_NULL for U and V if not interested in retrieve these matrices
//PetscErrorCode MatSVD(Mat& A, MPI_Comm comm, PetscInt& nconv, PetscReal*& sigma, Mat* U, Mat* V)
//multiplySigmaSqrt==1: U = U*sqrt(sigma), V = V*sqrt(sigma)
//default values: nconv=0, sigma=NULL, U=NULL, V=NULL
//*nconv: is the number of eig vector interested in, set *nconv=0 to get all convergence ones
//upon return, nconv is the actual number of eig vectors returned! 
//sigma: do not need to pre-allocated space for sigma, just PetscScalar**
#include "slepcsvd.h"
PetscErrorCode MatSVD(Mat* A, MPI_Comm comm, PetscInt* nconv, PetscScalar** sigma, Mat* U, Mat* V, PetscInt multiplySigmaSqrt)
{
  SVD            svd=PETSC_NULL;     /* singular value problem solver context */
  const SVDType  type=PETSC_NULL;
  PetscReal      error, tol;
  PetscErrorCode ierr;
  PetscInt       nconv_local, nsv, maxit, i, count, its, m_A, n_A, tmp_ia, tmp_ib;
  PetscViewer pviewer = PETSC_VIEWER_STDOUT_(comm);
  Vec u=PETSC_NULL, v=PETSC_NULL;                                   
  PetscInt *indices_v=PETSC_NULL, *indices_u=PETSC_NULL; //for create U and V, MatSetColumnVec
  PetscScalar* sigma_local=PETSC_NULL; //eigenvalue array
  PetscScalar sigma_scalar;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
               Allocate space for u, v and U, V
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = MatGetSize(*A, &m_A, &n_A); CHKERRQ(ierr);
  ierr = VecCreateMPI(comm,PETSC_DECIDE,m_A,&u); CHKERRQ(ierr);
  ierr = VecCreateMPI(comm,PETSC_DECIDE,n_A,&v); CHKERRQ(ierr);
  if(U!=PETSC_NULL){
    ierr = PetscMalloc(sizeof(PetscInt)*m_A, (void**) &indices_u); CHKERRQ(ierr);
  }
  if(V!=PETSC_NULL) {
    ierr = PetscMalloc(sizeof(PetscInt)*n_A, (void**) &indices_v); CHKERRQ(ierr);
  }
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                Create the singular value solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  //  Create singular value solver context
  ierr = SVDCreate(comm,&svd);CHKERRQ(ierr);

  // Set operator
  ierr = SVDSetOperator(svd,*A);CHKERRQ(ierr);

  
  //   Set solver parameters at runtime
  ierr = SVDSetFromOptions(svd);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                      Solve the singular value system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = SVDSolve(svd);CHKERRQ(ierr);
  if(*nconv>0){ //set the number of singular values to compute 
    ierr = SVDGetDimensions(svd,&nsv,&tmp_ia,&tmp_ib);CHKERRQ(ierr);
#ifdef PETSC_OPERATIONS_DEBUG
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested singular values: %d\n",nsv);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," the maximum dimension of the subspace to be used by the solver: %d\n",tmp_ia);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," the maximum dimension allowed for the projected problem : %d\n",tmp_ib);CHKERRQ(ierr);
#endif
    ierr = SVDSetDimensions(svd, *nconv, tmp_ia, tmp_ib);CHKERRQ(ierr);    
  }
  
  //   Optional: Get some information from the solver and display it
#ifdef PETSC_OPERATIONS_DEBUG
  ierr = SVDGetIterationNumber(svd, &its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %d\n",its);CHKERRQ(ierr);
  ierr = SVDGetType(svd,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = SVDGetDimensions(svd,&nsv,&tmp_ia,&tmp_ib);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested singular values: %d\n",nsv);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," the maximum dimension of the subspace to be used by the solver: %d\n",tmp_ia);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," the maximum dimension allowed for the projected problem : %d\n",tmp_ib);CHKERRQ(ierr);
  ierr = SVDGetTolerances(svd,&tol,&maxit);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%d\n",tol,maxit);CHKERRQ(ierr);
#endif

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  //TODO: [ 0 0 0 0; 0 0 0 0; 0 0 4 0; 0 0 0 0] -->nan!
  //sort the egi values and get rid of nans
  
  //   Get number of converged singular triplets
  ierr = SVDGetConverged(svd,&nconv_local);CHKERRQ(ierr);
#ifdef PETSC_OPERATIONS_DEBUG
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged approximate singular triplets: %d\n\n",nconv_local);CHKERRQ(ierr);
#endif
  if(*nconv==0){
    *nconv = nconv_local;
  }else if(*nconv > nconv_local){ //SVDSetDimensions is not necessrily strict! So not necessarily equal
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of Requested eigs (%d) is bigger than the number of converged approximate singular triplets: %d\n\n",*nconv, nconv_local);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Possible reason: matrix is too singular\n");CHKERRQ(ierr);
    *nconv = nconv_local; //trucate it
    assert(0);
  }else{//*nconv <=nconv_local
    //do nothing
  }

  if (*nconv>0) {
    //   Display singular values and relative errors
#ifdef PETSC_OPERATIONS_DEBUG
    ierr = PetscPrintf(PETSC_COMM_WORLD,
         "          sigma           residual norm\n"
         "  --------------------- ------------------\n" );CHKERRQ(ierr);
#endif
    if(sigma!=NULL){
      ierr = PetscMalloc(sizeof(PetscReal)*(*nconv),&sigma_local); CHKERRQ(ierr);
      *sigma = sigma_local; //asign to the output argument
    }
    if(U!=PETSC_NULL){ //return U
      //check if *U is allocated already and it dimension, U: m_A by *nconv
      if(*U!=PETSC_NULL){ //Mat already allocated
        MatGetSize(*U, &tmp_ia,&tmp_ib);
        if(tmp_ia!=m_A){
          ierr = PetscPrintf(PETSC_COMM_WORLD,"U:tmp_ia(%d)!=m_A(%d)\n",tmp_ia, m_A);CHKERRQ(ierr);
          assert(0);
        }
        if(tmp_ib!=*nconv){
          ierr = PetscPrintf(PETSC_COMM_WORLD,"U:tmp_ib(%d)!=nconv(%d)\n",tmp_ib, *nconv);CHKERRQ(ierr);
          assert(0);
        }
        MatScale(*U, 0.0); //set to be zeros
      }else{         
        MatCreateMPIDense(comm,PETSC_DECIDE,PETSC_DECIDE,m_A, (*nconv),PETSC_NULL, U);
      }
    }
    if(V!=PETSC_NULL){ //return V
      //check if *V is allocated already and it dimension, V: n_A by *nconv
      if(*V!=PETSC_NULL){ //Mat already allocated
        MatGetSize(*V, &tmp_ia,&tmp_ib);
        assert(tmp_ia==n_A); //V
        assert(tmp_ib==*nconv); //V
        MatScale(*V, 0.0);//set to be zeros
      }else{ 
        MatCreateMPIDense(comm,PETSC_DECIDE,PETSC_DECIDE,n_A, (*nconv), PETSC_NULL,V);
      }
    }
    
    //nconv_local is the total convergene solution
    
    for( i=0, count = 0; i<(*nconv) && count<nconv_local; i++, count++) {  
      //   Get converged singular triplets: i-th singular value is stored in sigma_scalar
      ierr = SVDGetSingularTriplet(svd,count,&sigma_scalar,u,v);CHKERRQ(ierr);
      //   Compute the error associated to each singular triplet 
      ierr = SVDComputeRelativeError(svd,count,&error);CHKERRQ(ierr);
      if(isnan(sigma_scalar) ||isnan(error)){ //skip this pair
        i--; //indices of U and V
        continue;
      }
      #ifdef PETSC_OPERATIONS_DEBUG
      ierr = PetscPrintf(PETSC_COMM_WORLD,"       % 6f      ",sigma_scalar); CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD," % 12g\n",error);CHKERRQ(ierr);
      #endif

      if(sigma!=NULL){
        (*sigma)[i]=sigma_scalar;
      }
        
      if(U!=PETSC_NULL){
        if(multiplySigmaSqrt==1){
          ierr=VecScale(u,sqrt(sigma_scalar));CHKERRQ(ierr);
          ierr=MatSetColumnVec(U, &u, indices_u, i, 0,0,m_A);CHKERRQ(ierr);
        }else{
          ierr=MatSetColumnVec(U, &u, indices_u, i, 0,0,m_A);CHKERRQ(ierr);
        }
        #ifdef PETSC_OPERATIONS_DEBUG
        PetscPrintf(PETSC_COMM_WORLD," -- SVD, u[%d]:",i); 
        VecView(u, pviewer);
        #endif
      }
      if(V!=PETSC_NULL){
        if(multiplySigmaSqrt==1){
          ierr=VecScale(v,sqrt(sigma_scalar));CHKERRQ(ierr);
          ierr=MatSetColumnVec(V, &v, indices_v, i, 0,0,n_A);CHKERRQ(ierr);
        }else{
          ierr=MatSetColumnVec(V, &v, indices_v, i, 0,0,n_A);CHKERRQ(ierr);
        }
        #ifdef PETSC_OPERATIONS_DEBUG
        PetscPrintf(PETSC_COMM_WORLD," -- SVD, v[%d]:",i); 
        VecView(v, pviewer);
        #endif
      }  
    }
    #ifdef PETSC_OPERATIONS_DEBUG
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n" );CHKERRQ(ierr);
    #endif
  }
  *nconv = i;
  if(U!=PETSC_NULL){
    MatAssemblyBegin(*U, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*U, MAT_FINAL_ASSEMBLY);
#ifdef PETSC_OPERATIONS_DEBUG
    MatView(*U,pviewer);CHKERRQ(ierr);
#endif
  }
  if(V!=PETSC_NULL){
    MatAssemblyBegin(*V, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*V, MAT_FINAL_ASSEMBLY);
#ifdef PETSC_OPERATIONS_DEBUG
    MatView(*V,pviewer);CHKERRQ(ierr);
#endif
  }

  
  //   Free work space
  ierr = SVDDestroy(svd);CHKERRQ(ierr);
  ierr = VecDestroy(u);CHKERRQ(ierr);
  ierr = VecDestroy(v);CHKERRQ(ierr);
  
  return ierr;
}

//Limitation: only application to one-rank environment!
//because of MatMatSolve! Need to install lapack! 
PetscErrorCode  MatInverse(Mat init_covar, Mat* A_inv, MPI_Comm comm)
{
  PetscInt dim , i, m_X, n_X;
  IS idx;
  PetscInt* is;
  Mat fact, eye, init_covar_inverse;
  PetscErrorCode ierr;
//  MPI_Comm comm = PETSC_COMM_WORLD;

  //initialize dim and is
  MatGetSize(init_covar, &m_X, &n_X);
  assert(m_X==n_X);
  dim=m_X;
  ierr = PetscMalloc(sizeof(PetscInt)*dim, &is);CHKERRQ(ierr); 
  for(i=0;i<dim; i++){
    is[i]=i;
  }
  
  //create identity matrix
  ierr = MatCreateMPIDense(comm, PETSC_DECIDE, PETSC_DECIDE, dim, dim, PETSC_NULL, &eye);CHKERRQ(ierr); 
  MatSetFromOptions( eye );
  for( i = 0; i < dim; ++i)
  {
    MatSetValue(eye, i, i, 1.0, INSERT_VALUES);
  }
  MatAssemblyBegin(eye,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(eye,MAT_FINAL_ASSEMBLY);

  //create solution matrix
  ierr = MatCreateMPIDense(comm, PETSC_DECIDE, PETSC_DECIDE, dim, dim, PETSC_NULL, &init_covar_inverse);CHKERRQ(ierr); 
  MatSetFromOptions( init_covar_inverse );
  
  /** cholesky decomposition **/
  ierr = ISCreateGeneral(comm, dim, is, &idx);CHKERRQ(ierr); 
  MatFactorInfo info;
  info.fill = 1.f;
  info.dtcol = 0.f;

  if(0){ //for symmetric matrix only
    ierr = MatGetFactor(init_covar,MAT_SOLVER_PETSC, MAT_FACTOR_CHOLESKY, &fact);CHKERRQ(ierr); 
    ierr = MatCholeskyFactorSymbolic(fact, init_covar, idx, &info);CHKERRQ(ierr);   
    ierr = MatCholeskyFactorNumeric( fact, init_covar, &info);CHKERRQ(ierr); 
  }else{//for all matrix
    ierr = MatGetFactor(init_covar,MAT_SOLVER_PETSC, MAT_FACTOR_LU, &fact);CHKERRQ(ierr); 
    ierr = MatLUFactorSymbolic(fact, init_covar, idx,idx, &info);CHKERRQ(ierr);   
    ierr = MatLUFactorNumeric( fact, init_covar, &info);CHKERRQ(ierr); 
  }  
  //if this fails --> the matrix is not invertible
  
  //MatMatSolve(init_covar,eye,init_covar_inverse);
  MatMatSolve(fact,eye, init_covar_inverse);
  MatAssemblyBegin(init_covar_inverse,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(init_covar_inverse,MAT_FINAL_ASSEMBLY);
  
#ifdef PETSC_OPERATIONS_DEBUG
  PetscPrintf(PETSC_COMM_WORLD, "-------MatInverse: eye ----------\n");
  ierr = MatView( eye, PETSC_VIEWER_STDOUT_WORLD );CHKERRQ(ierr); 
  PetscPrintf(PETSC_COMM_WORLD, "-------MatInverse: fact----------\n");
  ierr = MatView( fact, PETSC_VIEWER_STDOUT_WORLD );CHKERRQ(ierr); 
  PetscPrintf(PETSC_COMM_WORLD, "-------MatInverse: inverse-------\n");
  ierr = MatView( init_covar_inverse, PETSC_VIEWER_STDOUT_WORLD );CHKERRQ(ierr); 
#endif
  *A_inv  = init_covar_inverse;
  return ierr;  
}

//Cholesky decomposition, result overwritten in the same matrix
PetscErrorCode  MatFactorize(Mat& init_covar, MPI_Comm comm, int type=1)
{
  PetscInt dim , i, m_X, n_X, *is;
  IS idx;
  Mat fact;    
  PetscErrorCode ierr;

  //initialize dim and is
  MatGetSize(init_covar, &m_X, &n_X);
  assert(m_X==n_X);
  dim=m_X;
  ierr = PetscMalloc(sizeof(PetscInt)*dim, &is);CHKERRQ(ierr); 
  for(i=0;i<dim; i++){
    is[i]=i;
  }
    
  /** cholesky or LU decomposition **/
  ierr = ISCreateGeneral(comm, dim, is, &idx);CHKERRQ(ierr); 
  MatFactorInfo info;
  //info.fill = 1.f;
  //info.dtcol = 0.f;

  if(type==1){ //for symmetric matrix only
    ierr = MatGetFactor(init_covar,MAT_SOLVER_PETSC, MAT_FACTOR_CHOLESKY, &fact);CHKERRQ(ierr); 
    ierr = MatCholeskyFactorSymbolic(fact, init_covar, idx, &info);CHKERRQ(ierr);   
    ierr = MatCholeskyFactorNumeric( fact, init_covar, &info);CHKERRQ(ierr); 
  }else if(type==2){//for all matrix
    ierr = MatGetFactor(init_covar,MAT_SOLVER_PETSC, MAT_FACTOR_LU, &fact);CHKERRQ(ierr); 
    ierr = MatLUFactorSymbolic(fact, init_covar, idx,idx, &info);CHKERRQ(ierr);   
    ierr = MatLUFactorNumeric( fact, init_covar, &info);CHKERRQ(ierr); 
  }else{
    PetscPrintf(PETSC_COMM_WORLD, "-------MatFactorize: wrong type----------\n");
    assert(0);
  }  
  
  //fact->factor = 0;
  assert(  *((int*) ((void*)fact + 336)) == MAT_FACTOR_CHOLESKY);
  *((int*) ((void*)fact + 336)) = 0;  //hack it
  if(0){
    MatCopy(fact, init_covar, SAME_NONZERO_PATTERN);
  }if(0){
    Mat new_mat;
    MatDuplicate(fact, MAT_COPY_VALUES, &new_mat);
    MatDestroy(init_covar);
    init_covar = new_mat;
  }else{
    Vec tmp = petsc_create_vector(1, n_X, comm);
    Mat new_mat = petsc_create_matrix(1, m_X, n_X, comm);
    PetscInt* index = new int[m_X]; 
    for(int i=0;i<m_X;i++) index[i]=i;
    for(int i=0;i<n_X; i++){
      MatGetColumnVector(fact, tmp, i);
      MatSetColumnVec(&new_mat, &tmp,   index, i, i, i, m_X);
    }
    MatAssemblyBegin(new_mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd  (new_mat, MAT_FINAL_ASSEMBLY);
    MatDestroy(init_covar);
    VecDestroy(tmp);
    free(index);
    init_covar = new_mat;
  }

#ifdef PETSC_OPERATIONS_DEBUG
  PetscPrintf(PETSC_COMM_WORLD, "-------MatInverse: fact----------\n");
  ierr = MatView( fact, PETSC_VIEWER_STDOUT_WORLD );CHKERRQ(ierr); 
  PetscPrintf(PETSC_COMM_WORLD, "-------MatInverse: result returned----------\n");
  ierr = MatView( init_covar, PETSC_VIEWER_STDOUT_WORLD );CHKERRQ(ierr); 
#endif
  return ierr;  
}


//interploate a vector at discrete time points into a continuous time field
//y: nx by nt, x: 1 by nt, yi: 1 by nt
//all of the matrices/vectors should be pre-allocated
//example:
//  interp_v(0, NX, NT, x, (double*)y, -1, NULL);  
//  xi = 0.15;
//  interp_v(1, NX, NT, NULL, NULL, xi, (double*)yi);
//  interp_v(2, NX, NT, NULL, NULL, -1, NULL);
void interp_v(int usage, int nx, int nt, double* x, double *y , double xi, double* yi)
{
  static gsl_interp_accel  ** acc; //[NX]
  static gsl_spline        ** spline; //[NX]
  int j;
  
  if(usage==0){ //initialize
    acc = (gsl_interp_accel**) malloc(nx * sizeof(gsl_interp_accel*));
    spline = (gsl_spline**) malloc(nx* sizeof(gsl_spline*));
    for (j=0;j<nx;j++){
      acc[j]       = gsl_interp_accel_alloc ();
      spline[j]    = gsl_spline_alloc (gsl_interp_cspline, nt);
      gsl_spline_init (spline[j], x, &y[j*nt+0], nt);
    }
  }  

  if(usage==1){ //evaulation
    printf ("#interpolated\n");
    for(j=0;j<nx; j++){
      yi[j] = gsl_spline_eval (spline[j], xi, acc[j]);
      printf ("%04d: %5.5g, %5.5g\n", j, xi, yi[j]);
    }
  }
  
  if(usage==2){ //finalize
    for (j=0;j<nx;j++){    
      gsl_spline_free (spline[j]);
      gsl_interp_accel_free (acc[j]);
    }
    free(spline);
    free(acc);
  }
  
}


//----------------io--------------------------
//load a 2D matrix from an ascii file
//file format:
//comment starts with %
//%dimension: 4 by 4
//[matrix data]
// example:
//  save_matrix_2D_ascii("/users/jiaxi/matrix/tmp.txt", (double*) y, NX, NT, "written from test004.cpp");
//  int nx, nt;
//  double* matrix=load_matrix_2D_ascii("/users/jiaxi/matrix/tmp.txt", nx, nt);
double* load_matrix_2D_ascii(char* filename, int& size_1, int& size_2)
{
  FILE* fid = fopen(filename,"r");
  if(!fid) {printf("can not open %s\n",filename); size_1=size_2=0; return NULL;}
  printf("reading %s ... ",filename); 
  fscanf(fid,"%%dimension: %d %d",&size_1, &size_2);
  assert(size_1>0 && size_2>0);
  double* matrix = new double[size_1*size_2];
  for (int i=0; i<size_1*size_2; i++) fscanf(fid, "%lg", matrix+i);
  fclose(fid); printf("done!\n");   
  return matrix;
}
int save_matrix_2D_ascii(char* filename, double* matrix, int size_1, int size_2, char* comment)
{
  FILE* fid = fopen(filename,"w");
  if(!fid) {
    char* copy = strdup(filename);
    char* path_end_pos = strrchr(copy,'/');
    assert(path_end_pos); *path_end_pos='\0';    
    char tmp[256]; sprintf(tmp,"mkdir %s",copy); system(tmp); 
    fid = fopen(filename,"w");
    free(copy);
    if(!fid) printf("can not open %s\n",filename); return 1;
  }
  printf("writing %s ... ",filename);
  //if(comment) fprintf(fid, "%%%s\n", comment);
  fprintf(fid,"%%dimension: %d %d\n",size_1, size_2);
  assert(size_1>0 && size_2>0);
  for (int i=0; i<size_1*size_2; i++){
    fprintf(fid, "%10.15lg ", matrix[i]);
    if((i+1)%size_2==0) fprintf(fid, "\n");
  }
  fclose(fid); printf("done!\n");   
  return 0;
}
