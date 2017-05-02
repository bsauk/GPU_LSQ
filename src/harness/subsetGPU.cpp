#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "lsq.h"
#include "sub.h"

#include "CycleTimer.h"

/****************************************************************************************************
bsauk 4/18

While developing, this function is used to test the functions I am implementing for subset selection.
These functions will serve as an accuracy comparison later on when developing parallel algorithms.
Initial comparison is with the Fortran algorithms that Alan Miller provided on his website. 

I will try to use simple examples that test the conditionals of every algorithm being developed.
Timing will only consider the best subset selection techniques and not the generation part.

Inputs for this function are: 
1) File name for best subset selection and linear regression
2) Number of rows in input file
3) Number of variables in input file (should be number of columns -2, last columns are weight, and y values)
4) How many best subset selection results should be stored for each number of variables.
5) Maximum number of variables to be used in the regression.
6) Either 0 or 1 to determine if checking should be used. 1 for comparison against gold version 0 otherwise.

******************************************************************************************************/

//#define MAXVAR 1500 // Problem specific, number of variables in equation can adjust these values
//#define MAXCASES 7500 // Problem specific, number of trials

void compare_results(int first, int max_size, int nbest, int lopt_dim1, double** ressGold, double** ressGPU, int** loptGold, int** loptGPU);

void gpu_lsq(double* A, double* weights, double* y, int rows, int cols, int nbest, int max_size, double** ress, int** lopt, double* bound);
/*
void subset_gold(double* A, double* weights, double* y, int rows, int cols, int nbest, int max_size, double** ress, int** lopt, double* bound) {

  int nvar = cols-1, nobs = 0, r_dim = cols*(cols-1)/2, max_cdim = max_size*(max_size+1)/2;
  double sserr[1], D[cols], r[r_dim], rss[cols], rhs[cols], xrow[cols+1], work[cols], tol[cols];
  double sterr[max_size];
  int vorder[cols], row_ptr[cols], ifault[1], list[4], ier[1];

  bool lindep[cols], tol_set[1], rss_set[1];
  double vsmall = 2.225e-307, total_sumsq;
  double eps = 1e-14;
  
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

  
  for(int i=0; i<rows; i++) {
    xrow[0] = 1.0;
    for(int j=1; j<cols; j++) {
      xrow[j] = A[i*cols+j-1];
    }
    double startInclud = CycleTimer::currentSeconds();
    includ(weights[i], xrow, y[i], cols, D, r, rhs, sserr);
    double endInclud = CycleTimer::currentSeconds();
    std::cout << "Includ for " << i << " = " << 1000.f*(endInclud-startInclud) << " ms" << std::endl;
  }
  //  std::cout << "sserr = " << sserr[0] << std::endl;
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
  total_sumsq = rss[0];
  int first = 1;
  int last = cols;

  // The next part is that I will need to implement the different subset selection techniques, pick a few
  // Forward selection
  forwrd(first, last, ifault, cols, max_size, D, rhs, r, nbest, rss, bound, ress, vorder, lopt, rss_set, sserr, row_ptr, tol);
  /*
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
*/

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
  if(argc > 8) {
    std:: cout << "Too many arguments please only provide file name, rows, cols, nbest, and max. subset size" << std::endl;
  }
  int check = 0;

  //  if(argc == 7) {
  //    check = atoi(argv[6]);
  //  }

  // This part reads in the input file using the # of rows and columns
  const int rows = atoi(argv[2]);
  const int cols = atoi(argv[3])+1;
  const int nbest = atoi(argv[4]);
  const int nvar_max = atoi(argv[5]);
  // If fit_constant true, which i'm assuming for now add 1 to cols
  double Agold[rows*cols]; // Input matrix
  double Amagma[rows*cols];
  double Adn[rows*cols];

  double yGold[rows];
  double yMagma[rows];
  double yDN[rows];

  double weights[rows];
  bool fit_const = true;
  int max_size, lopt_dim1;
  double vlarge = 1.79769e308;
  
  if(fit_const) {
    max_size = nvar_max+1;
  } else {
    max_size = nvar_max;
  }
  lopt_dim1 = max_size*(max_size+1)/2;
  double boundGold[max_size];
  double boundGPU[max_size];
  
  double **ressGold = new double*[max_size]; // Input matrix
  for(int i=0; i<max_size; i++) {
    boundGold[i] = vlarge;
    ressGold[i] = new double[nbest];
    for(int j=0; j<nbest; j++) {
      ressGold[i][j] = vlarge;
    }
  }

  int **loptGold = new int*[nbest];
  for(int i=0; i<nbest; i++) {
    loptGold[i] = new int[lopt_dim1];
    for(int j=0; j<lopt_dim1; j++) {
      loptGold[i][j] = 0;
    }
  }
  
  double **ressGPU = new double*[max_size]; // Input matrix
  for(int i=0; i<max_size; i++) {
    boundGPU[i] = vlarge;
    ressGPU[i] = new double[nbest];
    for(int j=0; j<nbest; j++) {
      ressGPU[i][j] = vlarge;
    }
  }
  
  int **loptGPU = new int*[nbest];
  for(int i=0; i<nbest; i++) {
    loptGPU[i] = new int[lopt_dim1];
    for(int j=0; j<lopt_dim1; j++) {
      loptGPU[i][j] = 0;
    }
  }
  
  std::ifstream file;
  file.open(argv[1]);
  for(int i=0; i<rows; i++) {
    for(int j=0; j<cols+1; j++) {
      if(j<cols-1) {
	file >> Agold[i*cols+j];
	//	Amagma[i*cols+j] = A[i*cols+j];
	//	Adn[i*cols+j] = A[i*cols+j];
      } else if(j==cols-1) {
	file >> weights[i];
      } else {
	file >> yGold[i];
	//	b[i] = y[i];
      }
    } 
  }
  file.close();

  //  memcpy(Amagma, Agold, rows*cols*sizeof(double));
  //  memcpy(Adn, Agold, rows*cols*sizeof(double));
  //  memcpy(yMagma, yGold, rows*sizeof(double));
  //  memcpy(yDN, yGold, rows*sizeof(double));
  // CPU sequential subset methodology copied from Fortran implementation
  double startGold = CycleTimer::currentSeconds();
  for(int i=0; i<1; i++) {
    gpu_lsq(Agold, weights, yGold, rows, cols, nbest, max_size, ressGold, loptGold, boundGold);
  }
  double endGold = CycleTimer::currentSeconds();
  //  if(check) compare_results(1, max_size, nbest, lopt_dim1, ressGold, ressGPU, loptGold, loptGPU);
  
  std::cout << "Overall Time: " << 1000.f*(endGold-startGold) << " ms" << std::endl;

  /****************************************************************************************************************************************
  // The following parts below are for the MAGMA and cuSolverDN comparisons. Currently commented out while I work on parallelizing my code.


  for(int i=0; i<1; i++) {
    dn_forwrd(rows, cols, Adn, yDN, max_size);
  }
  std::cout << "cuSolverDN didn't seg fault?" << std::endl;

  for(int i=0; i<1; i++) {
    magma_forwrd(rows, cols,  Amagma, yMagma, max_size);
  }
  std::cout << "MAGMA didn't seg fault?" << std::endl;
  ****************************************************************************************************************************************/  

  
  return 0;
}



