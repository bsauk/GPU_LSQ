/*********************************************************************************************************************
This file accepts a .dat file or .mtx file that has a matrix of data for a linear least squares problem.
This version is a recreation of the Alan Miller Fortran implementation of lsq.f90. 
LSQ_Gold is a sequential version used to compare parallel versions against for accuracy purposes.


Accept argv[1] = file name, argv[2] = # of rows, argv[3] = # of cols
Input file should have two extra columns at the end of 
1) weight 
2) b data

Weights need to be specified, if solving an unweighted linear least squares problem use weights of 1.0.
LSQ_Gold is designed to generate sequential values. 
This method will give values to compare parallel version with to ensure we are always accurate.
After ensuring accuracy, can disable to improve performance of testing.
For developing ensure that the solutions to gold and parallel are always the same.

MAXVAR is a parameter specifying the maximum number of variables allowed to be used in the problem, adjusted for memory
MAXCASES is the maximum number of input rows that can be accepted, again set to limit space but can be increased.

Developed by bsauk on 4/17/17
***********************************************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "lsq.h"

#define MAXVAR 30 // Problem specific, number of variables in equation can adjust these values
#define MAXCASES 500 // Problem specific, number of trials

void lsq_gold(double** A, double* weights, double* y, int rows, int cols) {
  int nvar = cols, n = 4, pos1 = 2, nobs = 0, in, r_dim = cols*(cols-1)/2, max_cdim = MAXVAR*(MAXVAR+1)/2;
  double sserr[1], D[cols], r[r_dim], rss[cols], rhs[cols], xrow[cols+1], work[cols], tolerances[cols], ycorr[MAXVAR], cormat[max_cdim], beta[MAXVAR], xx[MAXVAR];
  double sterr[MAXVAR];
  int vorder[cols], row_ptr[cols], ifault[1], list[4];
  bool lindep[cols], tol_set[1], rss_set[1];
  double vsmall = 2.225e-307;
  double t[MAXVAR], fitted, resid[MAXCASES], std_resid[MAXCASES], hii[1], std_err_pred, covmat[max_cdim];
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
  sing(lindep, ifault, cols, D, tol_set, r, tolerances, row_ptr, rhs, sserr, work);
  
  if(ifault[0] == 0) {
    std::cout << "QR-factorization is not singular" << std::endl;
  } else {
    for(int i=0; i<nvar; i++) {
      if(lindep[i]) 
	std::cout << vorder[i] << " is exactly linearly related to earlier variables" << std::endl;
    }
  }
  in = 1;
  partial_corr(cols, ifault, in, max_cdim, D, r, sserr, rhs, cormat, ycorr);

  // Need to implement how to use list, n, and pos1. If I am testing all combinations that will be how I do this!
  // For now use examples from demo.
  list[0] = 2;
  list[1] = 4;
  list[2] = 5;
  list[3] = 7; 
  reordr(list, n, pos1, ifault, cols, vorder, rss_set, rss, row_ptr, D, r, sserr, tolerances, rhs);
  tolset(cols, work, r, tolerances, tol_set);
  int nreq = 5;
  regcf(beta, nreq, ifault, cols, work, r, tolerances, tol_set, D, rhs, row_ptr);
  ss(cols, sserr, rss, rss_set, D, rhs);
  double var = rss[nreq-1] / (nobs - nreq);

  cov(nreq, var, covmat, max_cdim, sterr, ifault, cols, D, vsmall, nobs, rss, rss_set, rhs, row_ptr, sserr, r);

  //  std::cout << " Variable     Regn.Coeff    Std.Error   t-value   Res.Sum of Squares" << std::endl;
  // Calculate t values
  //  for(int i=0; i<nreq; i++) {
  //    t[i] = beta[i]/sterr[i];
  //    std::cout << "i = " << i << ", vo = " << vorder[i] << ", beta = " << beta[i] << ", sterr = " << sterr[i] << ", t = " << t[i] << ", rss = " << rss[i] << std::endl;
  //  }
  /*
  std::cout << "Covariances of parameter estimates" << std::endl;
  int i2 = nreq, i1;
  for(int i=0; i<nreq; i++) {
    i1 = i2+1;
    i2 = i2+nreq-i;
    std::cout << vname[vorder[i]] << std::endl;
    for(int j=i1; j<i2; j++) {
      std::cout << covmat[j] << std::endl;
    }
    std::cout << // write output is what fortran code does
  */  
  // Now delete variable with the smallest t-value by moving it to position 5 and then repeating the calculation for the next 3 variables!!
  // This part is where I will have to implement how I am modifying values and such ***TODO***

  int counter=0;
  double tmin = fabs(t[0]);
  for(int i=1; i<nreq; i++) {
    if(fabs(t[i]) < tmin) {
      counter=i;
      tmin = fabs(t[i]);
    }
  }
  //  std::cout << "Removing variable in position " << counter << std::endl;
  vmove(counter, nreq, ifault, cols, rss_set, rss, sserr, row_ptr, D, rhs, r, vorder, tolerances);
  nreq--;
  regcf(beta, nreq, ifault, cols, work, r, tolerances, tol_set, D, rhs, row_ptr);
  ss(cols, sserr, rss, rss_set, D, rhs);
  cov(nreq, var, covmat, max_cdim, sterr, ifault, cols, D, vsmall, nobs, rss, rss_set, rhs, row_ptr, sserr, r);

  //  std::cout << " Variable     Regn.Coeff    Std.Error   t-value   Res.Sum of Squares" << std::endl;
  // Calculate t values
  for(int i=0; i<nreq; i++) {
    t[i] = beta[i]/sterr[i];
    //    std::cout << vorder[i] << beta[i] << sterr[i] << t[i] << rss[i] << std::endl;
  }

  var = rss[nreq]/(nobs-nreq);

  double stdev_res = sqrt(var);
  double r2 = 1.0-rss[nreq]/rss[0];
  fitted = 0;
  //  std::cout << "State   Actual   Fitted   Residual   Std.Resid   SE(prediction)" << std::endl;
  for(int i=0; i<nobs; i++) {
    xx[0] = 1.0;
    for(int j=0; j<nreq; j++) {
      xx[j] = A[i][vorder[j]]; //Terrible efficiency :/
    }
    for(int k=0; k<nreq; k++) {
      fitted = fitted+beta[k]*xx[k];
    }
    resid[i] = y[i]-fitted;
    hdiag(xx, nreq, hii, ifault, cols, D, tolerances, r);
    std_resid[i] = resid[i]/sqrt(var*(1-hii[0]));
    std_err_pred = sqrt(varprd(xx, nreq, cols, sserr, nobs, D, tolerances, r));
  }
}

// Will need to include some type of header file with all of my other functions defined in it
int main(int argc, char* argv[]) {
  std::cout.precision(16); //16 digit precision to match up with Fortran implementation move to header

  // Error handling for # of inputs  
  if(argc == 1) {
    std::cout << "Please provide .dat file and number of rows and cols!\n" << std::endl;
    return 0;
  }
  if(argc == 2) {
    std::cout << "Please provide number of rows and cols!\n" << std::endl;
    return 0;
  }
  if(argc == 3){
    std::cout << "Please provide number of cols!\n" << std::endl;
    return 0;
  }
  if(argc > 4) {
    std:: cout << "Too many arguments please only provide file name, rows, cols" << std::endl;
  }

  // This part reads in the input file using the # of rows and columns
  const int rows = atoi(argv[2]);
  const int cols = atoi(argv[3])+1;
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
  //LSQ implementation
  lsq_gold(A, weights, y, rows, cols);
  
}
