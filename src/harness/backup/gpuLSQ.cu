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

#define NB 16
#define COLUMNS 256

/*
static inline int updiv(int n, int d) {
  return (n+d-1)/d;
}
*/

__device__ __inline__ int dmin(int a, int b) {
  if(a>b) {
    return b;
  } else {
    return a;
  }
}

// As of 5:03pm on 5/8 this version calculates dR and dA correctly. However, currently dWeights, and dY are not outputting correctly.
// This also causes dSSERR to output incorrectly because of incorrect inputs. It also appears that dD and dRHS are not currently calculated correctly. 
// Update as of 11:45am 5/9 this version calculates all values correctly when NB = Matrix Size. Now working on cases when NB < Matrix Size.

__global__ void includGPU(int rows, int cols, double* dA, double* dY, double* dD, double* dR, double* dRHS, double* dSSERR, double* dWeights, int r_dim) {

  __shared__ double dXblock[(NB)*COLUMNS];
  __shared__ double sD[COLUMNS];
  __shared__ double sRHS[COLUMNS];

  const int idx = blockIdx.x*blockDim.x+threadIdx.x; // Maps to rows
  const int jdx = blockIdx.y*blockDim.y+threadIdx.y; // Maps to columns
  double vsmall = 2.225e-307;
  int perRow = dmin(COLUMNS, cols); // Currently should not work accurately if COLUMNS < cols
  double w = 0.0, xk = 0.00, di = 0.00, cbar = 0.00, sbar = 0.00, xi = 0.00, tempR = 0.00, RHSi = 0.0, xy = 0.00, yi = 0.00;
  if(idx >= blockDim.x || jdx >= blockDim.y ) return;
  if(threadIdx.x == 0) {
    for(int i=threadIdx.y; i<perRow; i+=blockDim.y) { 
      sD[i] = dD[i];
      sRHS[i] = dRHS[i];
    }
  }

  for(int i=threadIdx.x; i<rows; i+=blockDim.x) { // i<rows
    for(int j=threadIdx.y; j<perRow; j+=blockDim.y) {
      dXblock[threadIdx.x*perRow+j] = dA[i*cols+j];
    }
    int rowsLeft = dmin(NB, rows-i+threadIdx.x);
    int nextr = 0;
    bool smallW = false;
    w = dWeights[i];
    yi = dY[i];
    
    for(int j=0; j<cols; j++) { // j < cols
      __syncthreads();
      if(fabs(w) < vsmall) {
	dWeights[i] = w;
	dY[i] = yi;
	smallW = true;
      } else {	
	di = sD[j];
	RHSi = sRHS[j];
      }
      //      for(int k=i-threadIdx.x; k<i+1; k++) {
      for(int k=0; k<threadIdx.x+1; k++) {
	xi = dXblock[k*perRow+j]; 
	w = dWeights[k+i-threadIdx.x];
	if(fabs(xi) >= vsmall && !smallW && fabs(w) >= vsmall) {
	  yi = dY[k+i-threadIdx.x];
	  cbar = di/(di+w*xi*xi);
	  sbar = w*xi/(di+w*xi*xi);
	  di = di+w*xi*xi;
	  for(int colBlock=jdx; colBlock<perRow; colBlock+=blockDim.y) {
	    if(colBlock > j) {
	      tempR = dR[nextr+colBlock-j-1];
	      xk = dXblock[k*perRow+colBlock];
	      if(k == threadIdx.x) {
		dXblock[k*perRow+colBlock] = xk-xi*tempR;
		dR[nextr+colBlock-j-1] = cbar*tempR+sbar*xk;
	      }
	      tempR = cbar*tempR+sbar*xk;
	    }
	  }
	  w = cbar*w;
	  xy = yi;
	  yi = xy-xi*RHSi;
	  RHSi = cbar*RHSi+sbar*xy;
	}
	__syncthreads();
      }
      nextr = nextr+cols-j-1;
      for(int rowBlock=0; rowBlock<NB; rowBlock++) {
	if(!smallW) {
	  for(int l=jdx; l<perRow; l+=blockDim.y) {
	    if(l == j && rowBlock == threadIdx.x) {
	      sD[l] = di;
	      sRHS[l] = RHSi;
	    }
	  }
	}
      }
      for(int colBlock=threadIdx.y; colBlock<perRow; colBlock+=blockDim.y) {
	if(colBlock == perRow-1 && fabs(xi) >= vsmall) {
	  dWeights[i] = w;
	  dY[i] = yi;
	}
      }	
    }
    for(int colBlock=threadIdx.y; colBlock<perRow; colBlock+=blockDim.y) {
      dA[i*cols+colBlock] = dXblock[threadIdx.x*perRow+colBlock];
    }
    // This will move the values stored in the shared state to the global variables, that I will need later!
    for(int j=threadIdx.y; j<perRow; j+=blockDim.y) {
      dD[j] = sD[j];
      dRHS[j] = sRHS[j];
    }
  }
  __syncthreads();

  /* // Used to test accuracy of the parallel code
  if(idx==0 && jdx == 0) {
    for(int i=0; i<r_dim; i++) {
      printf("dR[%d]=%f\n", i, dR[i]);
    }
    for(int i=0; i<cols; i++) {
      printf("D[%d]=%f rhs[%d]=%f\n", i, dD[i], i, dRHS[i]);
    }
  }
  */

  if(jdx==0 && idx==0) { // Have to sequentially add the dSSERR values because atomic_dadd doesn't seem to work in CUDA.
    for(int i=0; i<rows; i++) {
      dSSERR[0] = dSSERR[0]+dWeights[i]*dY[i]*dY[i];
      //      printf("dSSERR[%d]=%f dWeights[%d]=%f dY[%d]=%f\n", i, dSSERR[0], i, dWeights[i], i, dY[i]);
    }
  }
} 

void gpu_lsq(double* A, double* weights, double* y, int rows, int cols, int nbest, int max_size, double** ress, int** lopt, double* bound) {

  int nvar = cols-1, r_dim = cols*(cols-1)/2;
  double sserr[1], rss[cols], rhs[cols], work[cols], tol[cols], D[cols], r[r_dim];
  int vorder[cols], row_ptr[cols], ifault[1], ier[1];

  bool lindep[cols], tol_set[1], rss_set[1];
  //  double total_sumsq;
    
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

  double* dA = NULL;
  double* dY = NULL;
  double* dD = NULL; // cols
  double* dR = NULL; //r_dim
  double* dRHS = NULL; //cols
  double* dSSERR = NULL; // cols
  double* dWeights = NULL; //rows

  cudaMalloc((void **)&dA, rows*(cols+1)*sizeof(double));
  cudaMalloc((void **)&dY, rows*sizeof(double));
  cudaMalloc((void **)&dD, cols*sizeof(double));
  cudaMalloc((void **)&dR, r_dim*sizeof(double));
  cudaMalloc((void **)&dRHS, cols*sizeof(double));
  cudaMalloc((void **)&dSSERR, sizeof(double));
  cudaMalloc((void **)&dWeights, rows*sizeof(double));

  cudaMemcpy(dA, A, rows*(cols+1)*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dY, y, rows*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dWeights, weights, rows*sizeof(double), cudaMemcpyHostToDevice); // May want to consider just assuming 1 for now if this takes too long :/.
  
  cudaMemset(dD, 0.00, cols*sizeof(double));
  cudaMemset(dR, 0.00, r_dim*sizeof(double));
  cudaMemset(dRHS, 0.00, cols*sizeof(double));
  cudaMemset(dSSERR, 0.00, sizeof(double));

  dim3 threadsPerBlock(NB,NB);
  dim3 blocks(1, 1);
  includGPU<<<blocks, threadsPerBlock>>>(rows, cols, dA, dY, dD, dR, dRHS, dSSERR, dWeights, r_dim);
  cudaDeviceSynchronize();

  // Transfer results back to CPU from GPU!
  cudaMemcpy(D, dD, cols*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(r, dR, r_dim*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(rhs, dRHS, cols*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(sserr, dSSERR, sizeof(double), cudaMemcpyDeviceToHost);

  // Then rest of code is CPU code. Since includ was the most time consuming step (>90% of computation this should be ok)
  /****************************************************************
  // This part gets translated into CUDA device code.
  for(int i=0; i<rows; i++) {
    xrow[0] = 1.0;
    for(int j=1; j<cols; j++) {
      xrow[j] = A[i*cols+j-1];
    }
    includ(weights[i], xrow, y[i], cols, D, r, rhs, sserr);
  }
  *****************************************************************/
  /*  //Used to determine when values are accurate  
  for(int i=0; i<cols; i++) {
    printf("D[%d]=%f rhs[%d]=%f\n", i, D[i], i, rhs[i]);
  }
  for(int i=0; i<r_dim; i++) {
    printf("r[%d]=%f\n", i, r[i]);
  }
  printf("sserr=%f\n", sserr[0]);
  */
  
  //  std::cout << "sserr = " << sserr[0] << std::endl;
  sing(lindep, ifault, cols, D, tol_set, r, tol, row_ptr, rhs, sserr, work);

  if(ifault[0] == 0) {
    std::cout << "QR-factorization is not singular" << std::endl;
  } else {
    for(int i=0; i<nvar; i++) {
      if(lindep[i]) 
	std::cout << vorder[i] << " is exactly linearly related to earlier variables" << std::endl;
    }
  }

  ss(cols, sserr, rss, rss_set, D, rhs);
  
  // Set tolerances and test for singularities
  tolset(cols, work, r, tol, tol_set);
  sing(lindep, ier, cols, D, tol_set, r, tol, row_ptr, rhs, sserr, work);

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

  // Not sure if these three need to be called again here...
  tolset(cols, work, r, tol, tol_set);
  sing(lindep, ier, cols, D, tol_set, r, tol, row_ptr, rhs, sserr, work);
  ss(cols, sserr, rss, rss_set, D, rhs);
  
  for(int i=0; i<max_size; i++) {
    report(i, rss[i], max_size, bound, nbest, ress, vorder, lopt);
  }
  
  //  total_sumsq = rss[0];
  int first = 1;
  int last = cols;

  // The next part is that I will need to implement the different subset selection techniques, pick a few
  // Forward selection
  //  double startForwrd = CycleTimer::currentSeconds();
  forwrd(first, last, ifault, cols, max_size, D, rhs, r, nbest, rss, bound, ress, vorder, lopt, rss_set, sserr, row_ptr, tol);
  //  double endForwrd = CycleTimer::currentSeconds();
  //  std::cout << "Forwrd took " << 1000.f*(endForwrd-startForwrd) << std::endl;

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

