#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <assert.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

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

void dn_dgels(cublasHandle_t cublasH, cusolverDnHandle_t cudenseH, int m, int n, double* dA, int lda, double* dB, double* newErr) {

  double *d_tau;  // linear memory of gpu                                                 
  int *devInfo = NULL; //info in GPU (device copy)                                        
  double *d_work = NULL;                                                                  
  double *dC = NULL; // Intermediate vector on GPU nxnrhs                                
  int lwork = 0;                                                                          
  int info_gpu = 0;                                                                       
  const double alpha = 1;                                                                 
  const double beta = 0;                                                                  
  cudaMalloc((void **)&d_tau, sizeof(double)*m);                              
  cudaMalloc((void **)&devInfo, sizeof(int));                                 
  cudaMalloc((void **)&dC, sizeof(double)*n);                           
  int incx = 1;                                                                           
  int incy = 1;                                                                           
  cublasStatus_t cublas_status = CUBLAS_STATUS_SUCCESS;                                   
  cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;                             

  cublas_status = cublasDgemv(cublasH, CUBLAS_OP_T, m, n, &alpha, dA, lda, dB, incx, &beta, dC, incy);
  assert(cublas_status == CUBLAS_STATUS_SUCCESS);                                         
  
  cusolverDnDgeqrf_bufferSize(cudenseH, m, n, dA, lda, &lwork);
  cudaMalloc((void **)&d_work, sizeof(double)*lwork);                         

  cusolver_status = cusolverDnDgeqrf(cudenseH, m, n, dA, lda, d_tau, d_work, lwork, devInfo);                                                     
  assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);                                     
  cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost);        
  
  cublas_status = cublasDtrsv(cublasH, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, dA, lda, dC, incx);
  assert(CUBLAS_STATUS_SUCCESS == cublas_status);                                         
  cublas_status = cublasDtrsv(cublasH, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, n, dA, lda, dC, incx);
  assert(CUBLAS_STATUS_SUCCESS == cublas_status);  
  
  cublas_status = cublasDnrm2(cublasH, n, dC, incx, newErr);

}

void dn_forwrd(int m, int n, double* A, double* B, int max_size) {
  
  cusolverDnHandle_t cudenseH = NULL;  // defining handle                                 
  cublasHandle_t cublasH = NULL;                                                          
  cublasStatus_t cublas_status = CUBLAS_STATUS_SUCCESS;                                   
  cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;                             

  cusolver_status = cusolverDnCreate(&cudenseH);
  assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);                                     

  cublas_status = cublasCreate(&cublasH);                                                 
  assert(CUBLAS_STATUS_SUCCESS == cublas_status);                                         
                                                                                          
  cudaError_t cudaStat1 = cudaSuccess;                                                    
  cudaError_t cudaStat2 = cudaSuccess;                                                    

  double sserr = 1.0e100; //Chosen to be large
  int ldda = ((m+31)/32)*32;
  double *dA = NULL, *dB = NULL, *dTemp = NULL, *dPA = NULL, *dX = NULL;

  double *hTemp = (double *)malloc(ldda*max_size*sizeof(double));
  double *hPA = (double *)malloc(ldda*max_size*sizeof(double));
  double *X = (double *)malloc(n*sizeof(double));
  double newErr[1]; 
  cudaMalloc((void **)&dTemp, ldda*max_size*sizeof(double));
  cudaMalloc((void **)&dPA, ldda*max_size*sizeof(double));
  cudaMalloc((void **)&dX, ldda*sizeof(double));

  cudaStat1 = cudaMalloc((void **)&dA, m*n*sizeof(double));                              
  cudaStat2 = cudaMalloc((void **)&dB, m*sizeof(double));                           
  assert(cudaStat1 == cudaSuccess);                                                       
  assert(cudaStat2 == cudaSuccess);                                                       
  cudaStat1 = cudaMemcpy(dA, A, sizeof(double)*m*n, cudaMemcpyHostToDevice);
  cudaStat2 = cudaMemcpy(dB, B, sizeof(double)*m, cudaMemcpyHostToDevice);  
  assert(cudaStat1 == cudaSuccess);                                                       
  assert(cudaStat2 == cudaSuccess);                                                       


  bool chosen[n];
  for(int i=0; i<n; i++) {
    chosen[i] = false;
  }

  for(int vars=0; vars<max_size; vars++) {
    bool needVar = true;
    int tempChosen = 0;
    for(int col=0; col<n; col++) {
      if(needVar && !chosen[col]) {
	if(vars > 0) cudaMemcpy(dPA, dTemp, sizeof(double)*max_size*m, cudaMemcpyDeviceToDevice); // Reset dPA to be the dTemp 
	cudaMemcpy(dX, dB, sizeof(double)*ldda, cudaMemcpyDeviceToDevice);
	cudaMemcpy2D(&dPA[vars], sizeof(double), &dA[col], n*sizeof(double), sizeof(double), m, cudaMemcpyDeviceToDevice);
	dn_dgels(cublasH, cudenseH, m, vars+1, dPA, ldda, dX, newErr);
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

  cudaFree(dA); cudaFree(dB); cudaFree(dTemp); cudaFree(dPA); cudaFree(dX);
  free(hTemp); free(hPA); free(X);
  if(cublasH) cublasDestroy(cublasH);                                                     
  if(cudenseH) cusolverDnDestroy(cudenseH);                                               

}
