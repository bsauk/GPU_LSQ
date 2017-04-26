#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#include "magma.h"
#include "magma_lapack.h"

#define BLOCK_SIZE 512

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

__global__ void magmablas_dnrm2_kernel(int m, double *dA, int ldda, double *dxnorm) {
  const int tx = threadIdx.x;
  double *dx = dA+blockIdx.x*ldda;

  __shared__ double sum[BLOCK_SIZE];
  double lsum = 0;
  for(int j=tx; j<m; j+= BLOCK_SIZE) {
    double re = dx[j];
    lsum += re*re;
  }
  sum[tx] = lsum;
  magma_sum_reduce<BLOCK_SIZE>(tx,sum);
  if(tx==0)
    dxnorm[blockIdx.x] = sqrt(sum[0]);

}

void magma_forwrd(int m, int n, double* A, double* B, int max_size) {

  magma_init();
  magma_queue_t queues;
  magma_queue_create(&queues);
  double sserr;
  int ldda = ((m+31)/32)*32;
  double *dA = NULL, *dB = NULL, *dTemp = NULL, *dPA = NULL, *dX = NULL, *dErr = NULL;
  double *X = (double *)malloc(n*sizeof(double));
  cudaMalloc((void **)&dA, ldda*n*sizeof(double));
  cudaMalloc((void **)&dB, ldda*sizeof(double));
  cudaMalloc((void **)&dTemp, ldda*max_size*sizeof(double));
  cudaMalloc((void **)&dPA, ldda*max_size*sizeof(double));
  cudaMalloc((Void **)&dX, ldda*sizeof(double));
  cudaMalloc((void **)&dErr, sizeof(double));

  magma_dsetmatrix(m, n, A, m, dA, ldda);
  magma_dsetmatrix(m, 1, B, m, dB, ldda);

  bool chosen[n];
  for(int i=0; i<n; i++) {
    chosen[i] = false;
  }

  for(int vars=0; vars<max_size; vars++) {
    bool needVar = true;
    int nb = magma_get_dgeqrf_nb(m, vars+1);
    int lwork = (m-vars-1+nb)*(1+nb)+nb;
    int info;

    double *hwork = (double *)malloc(lwork*sizeof(double));
    int tempChosen = 0;
    for(int col=0; col<n; col++) {
      if(needVar && !chosen[col]) {
	//	copyColumns(dA, dTemp, m, col, n, vars); // This is a GPU function to copy column values... Super inefficient
	if(vars > 0) cudaMemcpy(dPA, dTemp, sizeof(double)*max_size*m, cudaMemcpyDeviceToDevice); // Reset dPA to be the dTemp that only contains the currently being used A cols.
	cudaMemcpy(dX, dB, sizeof(double)*ldda, cudaMemcpyDeviceToDevice);
	cudaMemcpyd2D(&dPA[vars], sizeof(double), &dA[col], n*sizeof(double), sizeof(double), m, cudaMemcpyDeviceToDevice);
	magma_dgels3_gpu(MagmaNoTrans, m, vars+1, 1, dA, ldda, dX, ldda, hwork, lwork, &info);
	magmablas_dnrm2_kernel(vars+1, dX, ldda, dErr); 
	cudaMemcpy(newErr, dErr, sizeof(double), cudaMemcpyDeviceToHost);
	if(newErr[0] < sserr) {
	  tempChosen = col;
	  sserr = newErr[0];
	} 
      }
      if(col==n-1) {
	chosen[tempChosen] = true;
	cudaMemcpyd2D(&dTemp[vars], sizeof(double), &dA[tempChosen], n*sizeof(double), sizeof(double), m, cudaMemcpyDeviceToDevice);
      }
    }   
  }
}
