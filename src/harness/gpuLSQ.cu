#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>

#include "CycleTimer.h"
#include "lsq.h"
#include "sub.h"

#define NB 1
#define COLUMNS 512

static inline int updiv(int n, int d) {
  return (n+d-1)/d;
}


__device__ __inline__ int dmin(int a, int b) {
  if(a>b) {
    return b;
  } else {
    return a;
  }
}

__global__ void includGPU(int rows, int cols, double* dA, double* dY, double* dD, double* dR, double* dRHS, double* dSSERR, double* dWeights, int r_dim) {

  extern __shared__ double dXblock[]; // Used shared memory based on the number of columns in a row, passed in function call.

  const int idx = blockIdx.x*blockDim.x+threadIdx.x; // Maps to rows
  const int jdx = blockIdx.y*blockDim.y+threadIdx.y; // Maps to columns
  double vsmall = 2.225e-307;
  int perRow = dmin(COLUMNS, cols); 
  double w = 0.0, xk = 0.00, di = 0.00, cbar = 0.00, sbar = 0.00, xi = 0.00, tempR = 0.00, RHSi = 0.0, xy = 0.00, yi = 0.00;
  if(idx >= blockDim.x || jdx >= blockDim.y ) return;
  
  for(int i=threadIdx.x; i<rows; i+=blockDim.x) { // Iterate over all rows
    for(int j=threadIdx.y; j<perRow; j+=blockDim.y) { // Iterate over all columns to get values into shared memory for that row
      dXblock[j] = dA[i*cols+j];
    }
    int nextr = 0;
    w = dWeights[i];
    yi = dY[i];
    
    for(int j=0; j<cols; j++) { //Iterate over all columns
      __syncthreads(); // Ensure that the previous iteration finishes before the next one begins. Possibly too tight of a constraint could relax if next value has finished.
      di = dD[j];  // Diagonal matrix
      RHSi = dRHS[j]; // RHS of equation
      if(fabs(w) < vsmall) { // If the weight is less than 1, go to next row.
	dWeights[i] = w;
	dY[i] = yi;
	break;
      } 
      xi = dXblock[j]; 
      if(fabs(xi) >= vsmall) {
	cbar = di/(di+w*xi*xi);
	sbar = w*xi/(di+w*xi*xi);
	di = di+w*xi*xi;
	for(int colBlock=jdx; colBlock<perRow; colBlock+=blockDim.y) { // This is how I have every thread update a value in the row, and then loop through for all values in the row.
	  if(colBlock > j) { // Only update if a value is larger than the current column.
	    tempR = dR[nextr+colBlock-j-1];
	    xk = dXblock[colBlock];
	    dXblock[colBlock] = xk-xi*tempR;
	    dR[nextr+colBlock-j-1] = cbar*tempR+sbar*xk;
	    tempR = cbar*tempR+sbar*xk;
	  }
	}
	// Update values here
	w = cbar*w;
	xy = yi;
	yi = xy-xi*RHSi;
	RHSi = cbar*RHSi+sbar*xy;
	for(int colBlock=threadIdx.y; colBlock<perRow; colBlock+=blockDim.y) { // This ensures that only one thread that has useful information updates the value in global mem.
	  if(colBlock == j) {
	    dD[colBlock] = di;
	    dRHS[colBlock] = RHSi;
	  }
	}
      }
      nextr = nextr+cols-j-1; // Deals with moving to new position for the R matrix.
    }
    if(jdx==0) {
      dWeights[i] = w;
      dY[i] = yi;
    }
  }
  /*  
  // Used to test accuracy of the parallel code
  if(idx==0 && jdx == 0) {
    for(int i=0; i<r_dim; i++) {
      printf("dR[%d]=%f\n", i, dR[i]);
    }
    for(int i=0; i<cols; i++) {
      printf("D[%d]=%f rhs[%d]=%f\n", i, dD[i], i, dRHS[i]);
    }
  }
  */     
  
  if(jdx==0 && idx==0) { // Have to sequentially add the dSSERR values because atomic_dadd doesn't seem to work in CUDA, even though it is in the documentation.
    for(int i=0; i<rows; i++) {
      dSSERR[0] = dSSERR[0]+dWeights[i]*dY[i]*dY[i];
    }
  }
} 

void gpu_lsq(double* A, double* weights, double* y, int rows, int cols, int nbest, int max_size, double** ress, int** lopt, double* bound, int check) {

  int nvar = cols-1, r_dim = cols*(cols-1)/2;
  double sserr[1], rss[cols], rhs[cols], work[cols], tol[cols], D[cols], r[r_dim];
  int vorder[cols], row_ptr[cols], ifault[1], ier[1];

  bool lindep[cols], tol_set[1], rss_set[1];
    
  sserr[0] = 0.0;
  tol_set[0] = false;
  rss_set[0] = false;

  for(int i=0; i<cols; i++) {
    vorder[i] = i;
  }
  row_ptr[0] = 0;
  for(int i=1; i<cols-1; i++) {
    row_ptr[i] = row_ptr[i-1]+cols-i; 
  }
  row_ptr[cols-1] = 0;

  double startCopy = CycleTimer::currentSeconds();

  double* dA = NULL;
  double* dY = NULL;
  double* dD = NULL; // cols
  double* dR = NULL; //r_dim
  double* dRHS = NULL; //cols
  double* dSSERR = NULL; // cols
  double* dWeights = NULL; //rows
  // Allocate and copy values to the GPU, were not timed.
  cudaMalloc((void **)&dA, rows*(cols+1)*sizeof(double));
  cudaMalloc((void **)&dY, rows*sizeof(double));
  cudaMalloc((void **)&dD, cols*sizeof(double));
  cudaMalloc((void **)&dR, r_dim*sizeof(double));
  cudaMalloc((void **)&dRHS, cols*sizeof(double));
  cudaMalloc((void **)&dSSERR, sizeof(double));
  cudaMalloc((void **)&dWeights, rows*sizeof(double));
  double endCopy = CycleTimer::currentSeconds();

  cudaMemcpy(dA, A, rows*(cols+1)*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dY, y, rows*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dWeights, weights, rows*sizeof(double), cudaMemcpyHostToDevice); // May want to consider just assuming 1 for now if this takes too long.
  
  cudaMemset(dD, 0.00, cols*sizeof(double));
  cudaMemset(dR, 0.00, r_dim*sizeof(double));
  cudaMemset(dRHS, 0.00, cols*sizeof(double));
  cudaMemset(dSSERR, 0.00, sizeof(double));

  dim3 threadsPerBlock(NB,updiv(512,NB));
  dim3 blocks(1, 1);
  int shared_size = (cols+1)*NB*sizeof(double);
  cudaDeviceSynchronize();
  //  printf("copyTime= %f ms\n", 1000.f*(endCopy-startCopy)); //Used to determine how long it takes to set up the matrix. Not used for comparison
  double startInclud = CycleTimer::currentSeconds();
  includGPU<<<blocks, threadsPerBlock, shared_size>>>(rows, cols, dA, dY, dD, dR, dRHS, dSSERR, dWeights, r_dim);
  cudaDeviceSynchronize();
  double endInclud = CycleTimer::currentSeconds();
  printf("%f\n", 1000.f*(endInclud-startInclud));
  // Transfer results back to CPU from GPU!
  cudaMemcpy(D, dD, cols*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(r, dR, r_dim*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(rhs, dRHS, cols*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(sserr, dSSERR, sizeof(double), cudaMemcpyDeviceToHost);

  sing(lindep, ifault, cols, D, tol_set, r, tol, row_ptr, rhs, sserr, work);
  /*
  if(check) {
    if(ifault[0] == 0) {
      std::cout << "QR-factorization is not singular" << std::endl;
    } else {
      for(int i=0; i<nvar; i++) {
	if(lindep[i]) 
	  std::cout << vorder[i] << " is exactly linearly related to earlier variables" << std::endl;
      }
    }
  }
  */
  ss(cols, sserr, rss, rss_set, D, rhs);
  
  // Set tolerances and test for singularities
  tolset(cols, work, r, tol, tol_set);
  sing(lindep, ier, cols, D, tol_set, r, tol, row_ptr, rhs, sserr, work);
  /*
  if(check) {
    if(ier[0] != 0) {
      std::cout << ier[0] << " singularities detected in predictor variables" << std::endl;
      std::cout << "These variables are linearly related to earlier ones:" << std::endl;
      for(int i=0; i<cols; i++) {
	if(lindep[i]) {
	  for(int j=0; j<nvar; j++) {
	    if(lindep[j]) {
	      std::cout << vorder[j] << std::endl;
	    }
	  }
	  break;
	}
      }
    }
  }
  */
  // Not sure if these three need to be called again here...
  tolset(cols, work, r, tol, tol_set);
  sing(lindep, ier, cols, D, tol_set, r, tol, row_ptr, rhs, sserr, work);
  ss(cols, sserr, rss, rss_set, D, rhs);
  
  for(int i=0; i<max_size; i++) {
    report(i, rss[i], max_size, bound, nbest, ress, vorder, lopt);
  }
  
  int first = 1;
  int last = cols;

  // The next part is that I will need to implement the different subset selection techniques, pick a few
  // Forward selection
  forwrd(first, last, ifault, cols, max_size, D, rhs, r, nbest, rss, bound, ress, vorder, lopt, rss_set, sserr, row_ptr, tol);
  if(check) {
    for(int i=first; i<max_size; i++) {
      std::cout << "Best subsets found of " << i << " variables" << std::endl;
      std::cout << "     R.S.S.          Variable numbers" << std::endl;
      int pos = (i*i+i)/2;
      for(int j=0; j<nbest; j++) {
	std::cout << ress[i][j] << "    ";
	for(int k=pos; k<pos+i+1; k++) {
	  std::cout << lopt[j][k] << "   ";
	}
	std::cout << std::endl;
      }
    }
  }
  // Free allocated memory afterwards.
  cudaFree(dA);
  cudaFree(dY);
  cudaFree(dD);
  cudaFree(dR);
  cudaFree(dRHS);
  cudaFree(dSSERR);
  cudaFree(dWeights);
}

