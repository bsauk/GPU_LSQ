// Accept argv[1] = file name, argv[2] = # of rows, argv[3] = # of cols
// Input file should have two extra columns of weights and b data at the end

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

// This will need to be parallelized
void includ(double weight, double* xrow, double y, int cols, double* D, double* r, double* rhs, double* sserr) {
  double w = weight;
  int nextr = 0;
  double xk = 0.00;
  double di = 0.00;
  double wxi = 0.00;
  double dpi = 0.00;
  double cbar = 0.00;
  double sbar = 0.00;
  double xi = 0.0;
  for(int i=0; i<cols; i++) {
    if(fabs(w) == 0.0)
      return;
    xi = xrow[i];
    if(fabs(xi) == 0.0) {
      nextr = nextr+cols-i-1;
    } else {
      di = D[i];
      wxi = w * xi;
      dpi = di + wxi*xi;
      cbar = di / dpi;
      sbar = wxi / dpi;
      w = cbar*w;
      D[i] = dpi;
      for(int k=i+1; k<cols; k++) {
	xk = xrow[k];
	xrow[k] = xk - xi*r[nextr];
	r[nextr] = cbar*r[nextr]+sbar*xk;
	nextr++;
      }
      xk = y;
      y = xk-xi*rhs[i];
      rhs[i] = cbar*rhs[i]+sbar*xk;
    }
  }
  sserr[0] = sserr[0]+w*y*y;
}    

void lsq(double** A, double* weights, double* y, int rows, int cols) {
  int nvar = cols;
  int nobs = 0;
  double sserr[1];
  *sserr = 0.0;
  double D[cols];
  int r_dim = cols*(cols-1)/2;
  double r[r_dim];
  double tol[cols];
  double rss[cols];
  int vorder[cols];
  int row_ptr[cols];
  double rhs[cols]; // Vector of RHS projections after scaling by sqrt(D).
  double xrow[cols+1]; // May need to check this part of the implementation for speed later
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
  
}

// Will need to include some type of header file with all of my other functions defined in it
int main(int argc, char* argv[]) {
  std::cout.precision(16); //16 digit precision to match up with Fortran implementation move to header

  // Error handling for # of inputs  
  if(argc == 1) {
    std::cout << "Please provide .dat file and number of rows and cols!" << std::endl;
    return 0;
  }
  if(argc == 2) {
    std::cout << "Please provide number of rows and cols!" << std::endl;
    return 0;
  }
  if(argc == 3){
    std::cout << "Please provide number of cols!" << std::endl;
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
  lsq(A, weights, y, rows, cols);
  
}
