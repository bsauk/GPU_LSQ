#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#include "magma.h"
#include "magma_lapack.h"


#define BLOCK_SIZE 512
#define SWAP 32
/******************************************************************************************
4/25 bsauk
This will be used for best subset selection with dgels as the linear least squares routine.
This method will use a greedy forward approach similiar to the one used by subset.f90.

The point of this is to have a fair comparison between best subset regression with dgels and
with LSQ. 

In this approach, I will solve LLSP varying the number of columns. I will find what the best
first variable is, then add the best second variable and so on repeating until the maximum
variable size.

I will choose to add variables if the sum of squared error is minimized. 
Either we need to copy columns or use tranpose...

*******************************************************************************************/
/*
__global__ void copyColumns(double *dA, double* dTemp, int rows, int col, int cols, int vars) {
  const int tid = threadIdx.x;
  if(tid > rows) return;
  dTemp[vars+tid*cols] = dA[tid*cols+col]; 
}
*/

double l2_norm(double const *u, int n) {
  double accum = 0.;
  for(int i=0; i<n; i++) {
    accum += u[i]*u[i];
  }
  return sqrt(accum);
}

void magma_forwrd(int m, int n, double* A, double* B, int max_size) {

  magma_init();
  magma_queue_t queues;
  magma_queue_create(&queues);
  double sserr = 1.0e100; //Chosen to be large
  int ldda = ((m+31)/32)*32;
  double *dA = NULL, *dB = NULL, *dTemp = NULL, *dPA = NULL, *dX = NULL;
  int nrhs = 1;
  double *hTemp = (double *)malloc(ldda*max_size*sizeof(double));
  double *hPA = (double *)malloc(ldda*max_size*sizeof(double));
  double *X = (double *)malloc(n*sizeof(double));
  double newErr[1]; 
  cudaMalloc((void **)&dA, ldda*n*sizeof(double));
  cudaMalloc((void **)&dB, ldda*sizeof(double));
  cudaMalloc((void **)&dTemp, ldda*max_size*sizeof(double));
  cudaMalloc((void **)&dPA, ldda*max_size*sizeof(double));
  cudaMalloc((void **)&dX, ldda*sizeof(double));

  magma_dsetmatrix(m, n, A, m, dA, ldda);
  magma_dsetmatrix(m, 1, B, m, dB, ldda);

  bool chosen[n];
  for(int i=0; i<n; i++) {
    chosen[i] = false;
  }

  for(int vars=0; vars<max_size; vars++) {
    int nb = magma_get_dgeqrf_nb(m, vars+1);
    int lwork = (m-vars-1+nb)*(1+nb)+nb;
    int info;
    double *hwork = (double *)malloc(lwork*sizeof(double));
    int tempChosen = 0;
    for(int col=0; col<n; col++) {
      if(!chosen[col]) {
	if(vars > SWAP) { // Arbitrarily chosen, can adjust. Point where we swap using LAPACK DGELS to MAGMA DGELS
	  memcpy(hPA, hTemp, sizeof(double)*max_size*m);
	  memcpy(X, B, ldda*sizeof(double));
	  for(int row=0; row<m; row++) {
	    hPA[vars+row*n] = A[col+row*n];
	  }
	  int vars1 = vars+1;
	  const char trans = 'N';
	  lapackf77_dgels(&trans, &m, &vars1, &nrhs, hPA, &ldda, X, &ldda, hwork, &lwork, &info);
	  newErr[0] = l2_norm(X, vars+1);
	  if(newErr[0] < sserr) {
	    tempChosen = col;
	    sserr = newErr[0];
	  }
	  if(col==n-1) {
	    chosen[tempChosen] = true;
	    for(int row=0; row<m; row++) {
	      hTemp[vars+row*n] = A[tempChosen+row*n];
	    }
	  }
	} else { // Use GPU
	  if(vars > 0) cudaMemcpy(dPA, dTemp, sizeof(double)*max_size*m, cudaMemcpyDeviceToDevice); // Reset dPA to be the dTemp 
	  cudaMemcpy(dX, dB, sizeof(double)*ldda, cudaMemcpyDeviceToDevice);
	  cudaMemcpy2D(&dPA[vars], sizeof(double), &dA[col], n*sizeof(double), sizeof(double), m, cudaMemcpyDeviceToDevice);
	  magma_dgels3_gpu(MagmaNoTrans, m, vars+1, 1, dPA, ldda, dX, ldda, hwork, lwork, &info);
	  magmablas_dnrm2_cols(m, vars+1, dPA, ldda, newErr);
	  if(newErr[0] < sserr) {
	    tempChosen = col;
	    sserr = newErr[0];
	  } 
	}
	if(col==n-1) {
	  chosen[tempChosen] = true;
	  cudaMemcpy2D(&dTemp[vars], sizeof(double), &dA[tempChosen], n*sizeof(double), sizeof(double), m, cudaMemcpyDeviceToDevice);
	}
      }   
    }
  }
  // for purposes of testing best subset
  for(int i=0; i<n; i++) {
    if(chosen[i]) std::cout << i << " ";
  }
  std::cout << std::endl;
  free(hTemp); free(hPA); free(X);
  cudaFree(dA); cudaFree(dB); cudaFree(dTemp); cudaFree(dPA); cudaFree(dX);
  magma_finalize();
  cudaDeviceReset();
}
