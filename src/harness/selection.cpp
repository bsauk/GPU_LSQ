#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "lsq.h"
#include "sub.h"

/****************************************************************************************************
bsauk 4/18

While developing, this function is used to test the functions I am implementing for subset selection.
These functions will serve as an accuracy comparison later on when developing parallel algorithms.
Initial comparison is with the Fortran algorithms that Alan Miller provided on his website. 

I will try to use simple examples that test the conditionals of every algorithm being developed.

******************************************************************************************************/

#define MAXVAR 30 // Problem specific, number of variables in equation can adjust these values
#define MAXCASES 500 // Problem specific, number of trials

void subset_gold(double** A, double* weights, double* y, int rows, int cols, int nbest, int nvar_max) {
  int nvar = cols-1, nobs = 0, in, r_dim = cols*(cols-1)/2, max_cdim = MAXVAR*(MAXVAR+1)/2;
  double sserr[1], D[cols], r[r_dim], rss[cols], rhs[cols], xrow[cols+1], work[cols], tol[cols], ycorr[MAXVAR], cormat[max_cdim], beta[MAXVAR], xx[MAXVAR];
  double sterr[MAXVAR];
  int vorder[cols], row_ptr[cols], ifault[1], list[4], i0, ier[1], max_size, lopt_dim1;

  bool lindep[cols], tol_set[1], rss_set[1], fit_const[1];
  double vsmall = 2.225e-307, total_sumsq;
  double t[MAXVAR], fitted, resid[MAXCASES], std_resid[MAXCASES], hii[1], std_err_pred, covmat[max_cdim];
  double eps = 1e-14;
  double vlarge = 3.4028e38;
  
  fit_const[0] = true;
  sserr[0] = 0.0;
  tol_set[0] = false;
  rss_set[0] = false;

  for(int i=0; i<cols; i++) {
    vorder[i] = i;
  }
  row_ptr[0] = 0;
  for(int i=1; i<cols; i++) {
    row_ptr[i] = row_ptr[i-1] + cols - i; 
  }

  if(fit_const) {
    max_size = nvar_max+1;
  } else {
    max_size = nvar_max;
  }
  lopt_dim1 = max_size*(max_size+1)/2;
  
  double bound[max_size];
  double **ress = new double*[max_size]; // Input matrix
  for(int i=0; i<max_size; i++) {
    bound[i] = vlarge;
    ress[i] = new double[nbest];
    for(int j=0; j<nbest; j++) {
      ress[i][j] = vlarge;
    }
  }
  
  int **lopt = new int*[lopt_dim1];
  for(int i=0; i<lopt_dim1; i++) {
    lopt[i] = new int[nbest];
    for(int j=0; j<nbest; j++) {
      lopt[i][j] = 0;
    }
  }
  
  for(int i=0; i<rows; i++) {
    for(int j=0; j<cols; j++) {
      if(j==0) 
	xrow[j] = 1.0;
      else
	xrow[j] = A[i][j-1];
    }
    includ(weights[i], xrow, y[i], cols, D, r, rhs, sserr);
  }
  nobs = rows;
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
  for(int i=0; i<cols; i++) {
    if(lindep[i]) {
      std::cout << ier << "singularities detected in predictor variables" << std::endl;
      std::cout << "These variables are linearly related to earlier ones:" << std::endl;
      for(int j=0; j<nvar; j++) {
	if(lindep[j]) {
	  std::cout << vorder[j] << std::endl;
	}
      }
    }
  }

  for(int i=0; i<max_size; i++) {
    report(i, rss[i], max_size, bound, nbest, ress, vorder, lopt);
  }
  
  i0 = 0;
  total_sumsq = rss[0];
  std::cout << total_sumsq << std::endl;
  // The next part is that I will need to implement the different subset selection techniques, pick a few
  // Forward selection
  int first = 0;
  int last = cols;
  forwrd(first, last, ifault, cols, max_size, D, rhs, r, nbest, rss, bound, ress, vorder, lopt, rss_set, sserr, row_ptr, tol);

  // Exhaustive search, which I will compare performance against on GPU implementation!
  

}

int main(int argc, char* argv[]) {
    std::cout.precision(16); //16 digit precision to match up with Fortran implementation move to header

  // Error handling for # of inputs  
  if(argc == 1) {
    std::cout << "Please provide .dat file and number of rows, cols, nbest, and max. subset size!" << std::endl;
    return 0;
  }
  if(argc == 2) {
    std::cout << "Please provide number of rows, cols, nbest, and max. subset size!" << std::endl;
    return 0;
  }
  if(argc == 3){
    std::cout << "Please provide number of cols, nbest, and max. subset size!" << std::endl;
    return 0;
  }
  if(argc == 4) {
    std:: cout << "Please provide nbest, and max. subset size!" << std::endl;
  }
  if(argc == 5) {
    std:: cout << "Please provide max. subset size!" << std::endl;
  }
  if(argc > 6) {
    std:: cout << "Too many arguments please only provide file name, rows, cols, nbest, and max. subset size" << std::endl;
  }

  // This part reads in the input file using the # of rows and columns
  const int rows = atoi(argv[2]);
  const int cols = atoi(argv[3])+1;
  const int nbest = atoi(argv[4]);
  const int nvar_max = atoi(argv[5]);
  // If fit_constant true, which i'm assuming for now add 1 to cols
  double **A = new double*[rows]; // Input matrix
  double A2[rows][cols]; //Input matrix for other methods
  double y[rows];
  double weights[rows];
  double b[rows]; // Output vector for Ax = b for other algorithms  
  
  std::ifstream file;
  file.open(argv[1]);
  for(int i=0; i<rows; i++) {
    A[i] = new double[cols];
    for(int j=0; j<cols+1; j++) {
      if(j<cols-1) {
	file >> A[i][j];
	A2[i][j] = A[i][j];
      } else if(j==cols-1) {
	file >> weights[i];
      } else {
	file >> y[i];
	b[i] = y[i];
      }
    } 
  }
  file.close();
  //Subset
  subset_gold(A, weights, y, rows, cols, nbest, nvar_max);
  return 0;
}



