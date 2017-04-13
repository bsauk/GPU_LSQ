// Accept argv[1] = file name, argv[2] = # of rows, argv[3] = # of cols
// Input file should have two extra columns of weights and b data at the end

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#define MAXVAR 30 // Problem specific, number of variables in equation can adjust these values
#define MAXCASES 500 // Problem specific, number of trials

// This will need to be parallelized
void includ(double weight, double* xrow, double y, int cols, double* D, double* r, double* rhs, double* sserr) {
  double vsmall = 1e-300;
  int nextr = 0;
  double w = weight;
  double xk = 0.00;
  double di = 0.00;
  double wxi = 0.00;
  double dpi = 0.00;
  double cbar = 0.00;
  double sbar = 0.00;
  double xi = 0.0;
  for(int i=0; i<cols; i++) {
    if(fabs(w) < vsmall)
      return;
    xi = xrow[i];
    if(fabs(xi) < vsmall) {
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

void ss(int cols, double* sserr, double* rss, bool* rss_set, double* D, double* rhs) {
  double total;
  total = sserr[0];
  rss[cols-1] = sserr[0];
  for(int i=cols-1; i>0; i--) {
    total = total + D[i]*pow(rhs[i],2);
    rss[i-1] = total;
  }
  rss_set[0] = true;
}

void vmove(int from, int to, int* ifault, int cols, bool* rss_set, double* rss, double* sserr, int* row_ptr, double* D, double* rhs, double* r, int* vorder, double* tolerances) {
  double d1, d2, x, d1new, d2new, cbar, sbar, y;
  int first, last, inc, m1, m2, pos;
  double vsmall = 1e-300;
  ifault[0] = 0;
  if(from < 0 || from >= cols) ifault[0] = ifault[0]+4;
  if(to < 0 || to >= cols) ifault[0] = ifault[0]+8;
  if(ifault[0]!=0) return;
  if(from==to) return;

  if(!rss_set[0])
    ss(cols, sserr, rss, rss_set, D, rhs);

  if(from < to) {
    first = from;
    last = to-1;
    inc = 1;
  } else {
    first = from-1;
    last = to;
    inc = -1;
  }

  for(int i=first; i!=last+inc; i+=inc) {
    m1 = row_ptr[i];
    m2 = row_ptr[i+1];
    d1 = D[i];
    d2 = D[i+1];
    //****************** Block of code is 40 ********
    if(fabs(d1) < vsmall && fabs(d2) < vsmall) {
      pos = i; // Standard swap
      for(int j=0; j<i-1; j++) { 
	x = r[pos];  // Store first as temp
	r[pos] = r[pos-1];
	r[pos-1] = x;
	pos = pos+cols-j-2;
      }
      m1 = vorder[i];
      vorder[i] = vorder[i+1];
      vorder[i+1] = m1;
      x = tolerances[i];
      tolerances[i] = tolerances[i+1];
      tolerances[i+1] = x;
      rss[i] = rss[i+1]+D[i+1]*pow(rhs[i+1],2);
    }
    //****************** Block of code is 40 ********
    x = r[m1];
    if(fabs(x)*sqrt(d1) < tolerances[i+1]) 
      x=0.0;
    if(fabs(d1) < vsmall || fabs(x) < vsmall) {
      D[i] = d2;
      D[i+1] = d1;
      r[m1] = 0.0;
      for(int j=i+1; j<cols; j++) {
	m1++;
	x = r[m1];
	r[m1] = r[m2];
	r[m2] = x;
	m2++;
      }
      x = rhs[i];
      rhs[i] = rhs[i+1];
      rhs[i+1] = x;
      pos = i; // Standard swap
      for(int j=0; j<i-1; j++) { 
	x = r[pos];  // Store first as temp
	r[pos] = r[pos-1];
	r[pos-1] = x;
	pos = pos+cols-j-2;
      }
      m1 = vorder[i];
      vorder[i] = vorder[i+1];
      vorder[i+1] = m1;
      x = tolerances[i];
      tolerances[i] = tolerances[i+1];
      tolerances[i+1] = x;
      rss[i] = rss[i+1]+D[i+1]*pow(rhs[i+1],2);
    } else if(fabs(d2) < vsmall) {
      D[i] = d1*pow(x,2);
      r[m1] = 1.0/x;
      for(int j=m1; j<m1+cols-i-2; j++) {
	r[j] = r[j]/x;
      }
      rhs[i] = rhs[i]/x;
      pos = i; // Standard swap
      for(int j=0; j<i; j++) { 
	x = r[pos];  // Store first as temp
	r[pos] = r[pos-1];
	r[pos-1] = x;
	pos = pos+cols-j-2;
      }
      m1 = vorder[i];
      vorder[i] = vorder[i+1];
      vorder[i+1] = m1;
      x = tolerances[i];
      tolerances[i] = tolerances[i+1];
      tolerances[i+1] = x;
      rss[i] = rss[i+1]+D[i+1]*pow(rhs[i+1],2);
    } 

    d1new = d2+d1*pow(x,2);
    cbar = d2/d1new;
    sbar = x*d1/d1new;
    d2new = d1*cbar;
    D[i] = d1new;
    D[i+1] = d2new;
    r[m1] = sbar;
    for(int j=i+1; j<cols; j++) {
      m1++;
      y = r[m1];
      r[m1] = cbar*r[m2]+sbar*y;
      r[m2] = y-x*r[m2];
      m2++;
    }
    y = rhs[i];
    rhs[i] = cbar*rhs[i+1]+sbar*y;
    rhs[i+1] = y-x*rhs[i+1];

    pos = i; // Standard swap
    for(int j=0; j<i; j++) { 
      x = r[pos];  // Store first as temp
      r[pos] = r[pos-1];
      r[pos-1] = x;
      pos = pos+cols-j-2;
    }
    m1 = vorder[i];
    vorder[i] = vorder[i+1];
    vorder[i+1] = m1;
    x = tolerances[i];
    tolerances[i] = tolerances[i+1];
    tolerances[i+1] = x;
    rss[i] = rss[i+1]+D[i+1]*pow(rhs[i+1],2);
  }
}

void reordr(int* list, int n, int pos1, int* ifault, int cols, int* vorder, bool* rss_set, double* rss, int* row_ptr, double* D, double* r, double* sserr, double* tolerances, double* rhs) {
  int next, tmp, val, val_order;
  ifault[0] = 0;
  if(n<0 || n>cols-pos1) 
    ifault[0] = ifault[0]+4;
  if(ifault[0] != 0)
    return;
  next = pos1-1;
  tmp = pos1-1;
  for(int i=0; i<cols; i++) {
    val = vorder[tmp];
    for(int j=0; j<n; j++) {
      if(val == list[j]) {
	if(tmp > next) {
	  vmove(tmp, next, ifault, cols, rss_set, rss, sserr, row_ptr, D, rhs, r, vorder, tolerances); // Need to implement still
	  next++;
	  if(next >= n+pos1) {
	    return;
	  } else {
	    break;
	  }
	} else if(next >= n+pos1) {
	  return;
	}
      }
    }
    tmp++;
    if(tmp > cols) {
      ifault[0] = 8;
      return;
    }
  }
}

void partial_corr(int cols, int* ifault, int in, int dimc, double* D, double* r, double* sserr, double* rhs, double* cormat, double* ycorr) {
  int base_pos, pos, pos1, pos2;
  double rms[cols-in];
  double sumxx, sumxy, sumyy;
  double work[cols-in];
  ifault[0] = 0;
  if(in < 0 || in > cols-1)
    ifault[0] = ifault[0] +4;
  if(dimc < (cols-in)*(cols-in-1)/2)
    ifault[0]=ifault[0]+8;
  if(ifault[0]!=0)
    return;
  base_pos = in*cols-(in+1)*(in+2)/2;
  
  if(D[in] > 0.0)
    rms[0] = 1.0/sqrt(D[in]);
  for(int i=in+1; i<cols; i++) {
    pos = base_pos+i+1;
    sumxx = D[i];
    for(int j=in; j<i-1; j++) {
      sumxx = sumxx+D[j]*pow(r[pos],2);
      pos = pos+cols-j-2;
    }
    if(sumxx > 0.0) {
      rms[i-1] = 1.0/sqrt(sumxx);
    } else {
      rms[i-1] = 0.0;
      ifault[0] = -i;
    }
  }
  
  sumyy = sserr[0];
  for(int i=in; i<cols; i++) {
    sumyy = sumyy + D[i]*pow(rhs[i],2);
  }
  if(sumyy > 0.0)
    sumyy = 1.0/sqrt(sumyy);
  
  pos = 0;
  for(int i=in; i<cols; i++) {
    sumxy = 0.0;
    for(int j=i-1; j<cols-in; j++) {
      work[j] = 0.0;
    }
    pos1 = base_pos + i;
    for(int j=in; j<i; j++) {
      pos2 = pos1+1;
      for(int k=i; k<cols; k++) {
	work[k-1] = work[k-1]+D[j]*r[pos1]*r[pos2];
	pos2++;
      }
      sumxy = sumxy+D[j]*r[pos1]*rhs[j];
      pos1 = pos1+cols-j-2;
    }
    pos2 = pos1+1;
    for(int j=i; j<cols; j++) {
      work[j-1] = work[j-1]+D[i]*r[pos2];
      pos2++;
      cormat[pos] = work[j]*rms[i]*rms[j];
      pos++;
    }
    sumxy = sumxy+D[i]*rhs[i];
    ycorr[i] = sumxy*rms[i-1]*sumyy;
  }
  for(int i=0; i<in; i++) {
    ycorr[i] = 0.0;
  }    
}

void tolset(int cols, double* work, double* r, double* tolerances, bool* tol_set) {
  double eps1 = 2.22e-15;
  int pos = 0;
  double total = 0;
  for(int i=0; i<cols; i++) {
    pos = i-1;
    total = work[i];
    for(int j=0; j<i; j++) {
      total = total + fabs(r[pos])*work[j];
      pos = pos+cols-j-2;
    }
    tolerances[i] = eps1 * total;
  }
  tol_set[0] = true;
}

void sing(bool* lindep, int* ifault, int cols, double* D, bool* tol_set, double* r, double* tolerances, int* row_ptr, double* rhs, double* sserr, double* work) {
  double temp, y, weight;
  double x[cols];
  int pos, pos2;
  ifault[0] = 0;
  for(int i=0; i<cols; i++) {
    work[i] = sqrt(fabs(D[i]));
  }
  if(!tol_set[0])
    tolset(cols, work, r, tolerances, tol_set);
  for(int i=0; i<cols; i++) {
    temp = tolerances[i];
    pos = row_ptr[i];
    lindep[i] = false;
    if(work[i] <= temp) {
      lindep[i] = true;
      ifault[0] = ifault[0] - 1;
      if(i<cols) {
	pos2 = pos+cols-i-2;
	x[0] = 0.0;
	for(int j=i+1; j<cols; j++) {
	  x[j] = r[pos+j-i-1];
	  r[pos+j-i-1] = 0.0;
	}
	y = rhs[i];
	weight = D[i];
	D[i] = 0.0;
	rhs[i] = 0.0;
	includ(weight, x, y, cols, D, r, rhs, sserr);
      } else {
	sserr[0] = sserr[0]+D[i]*pow(rhs[i],2);
      }
    }
  }
}

void regcf(double* beta, int nreq, int* ifault, int cols, double* work, double* r, double* tolerances, bool* tol_set, double* D, double* rhs, int* row_ptr) {
  int nextr;

  ifault[0] = 0;
  if(nreq < 1 || nreq > cols) ifault[0] = ifault[0]+4;
  if(ifault[0] != 0) return;

  if(!tol_set) tolset(cols, work, r, tolerances, tol_set);

  for(int i=nreq; i>0; i--) {
    if(sqrt(D[i] < tolerances[i])) {
      beta[i] = 0.0;
      D[i] = 0.0;
      ifault[0] = -i;
    } else {
      beta[i] = rhs[i];
      nextr = row_ptr[i];
      for(int j=i+1; j<nreq; j++) {
	beta[i] = beta[i] - r[nextr]*beta[j];
	nextr++;
      }
    }
  }
}

void inv(int nreq, double* rinv, int* row_ptr, double* r) {
  int pos, start, pos1, pos2;
  double total;
  pos = nreq*(nreq-1)/2;
  for(int i=nreq-1; i>=0; i--) {
    start = row_ptr[i];
    for(int j=nreq-2; i; j--) {
      pos1 = start;
      pos2 = pos;
      total = 0.0;
      for(int k=i; k<j-2; k++) {
	pos2 = pos2+nreq-k;
	total = total-r[pos1]*rinv[pos2];
	pos1++;
      }
      rinv[pos] = total-r[pos1];
      pos++;
    }
  }
}

void cov(int nreq, double var, double* covmat, int dimcov, double* sterr, int* ifault, int cols, double* D, double vsmall, int nobs, double* rss, bool* rss_set, double* rhs, int* row_ptr, double* sserr, double* r) {
  int dim_rinv, pos, start, pos2, pos1;
  dim_rinv = nreq*(nreq-1)/2;
  double rinv[dim_rinv];
  double total;
  ifault[0] = 0;
  if(dimcov < nreq*(nreq+1)/2) {
    ifault[0] = 1;
    return;
  }

  for(int i=0; i<nreq; i++) {
    if(fabs(D[i]) < vsmall) ifault[0] = -i;
  }
  if(ifault[0] != 0) return;

  if(nobs > nreq) {
    if(!rss_set) ss(cols, sserr, rss, rss_set, D, rhs);
    var = rss[nreq]/(nobs-nreq);
  } else {
    ifault[0] = 2;
    return;
  }

  inv(nreq, rinv, row_ptr, r);
    
  pos = 0;
  start = 0;
  for(int i=0; i<nreq; i++) {
    pos2 = start;
    for(int j=i; j<nreq; j++) {
      pos1 = start+j-i;
      if(j==i) {
	total = 1.0/D[j];
      } else {
	total = rinv[pos1-1] / D[j];
      }
      for(int k=j; k<nreq; k++) {
	total = total + rinv[pos1]*rinv[pos2]/D[k];
	pos1++;
	pos2++;
      }
      covmat[pos] = total*var;
      if(i==j)
	sterr[i] = sqrt(covmat[pos]);
      pos++;
    }
    start = start+nreq-i;
  }
}

void lsq(double** A, double* weights, double* y, int rows, int cols) {
  int nvar = cols, n = 4, pos1 = 2, nobs = 0, in, r_dim = cols*(cols-1)/2, max_cdim = MAXVAR*(MAXVAR+1)/2;
  double sserr[1], D[cols], r[r_dim], rss[cols], rhs[cols], xrow[cols+1], work[cols], tolerances[cols], ycorr[MAXVAR], cormat[max_cdim], beta[MAXVAR], xx[MAXVAR];
  double sterr[MAXVAR];
  int vorder[cols], row_ptr[cols], ifault[1], list[4];
  bool lindep[cols], tol_set[1], rss_set[1];
  double vsmall = 1e-300;
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
    std::cout << "QR-factorization is not singular\n" << std::endl;
  } else {
    for(int i=0; i<nvar; i++) {
      if(lindep[i]) 
	std::cout << vorder[i] << " is exactly linearly related to earlier variables\n" << std::endl;
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
  // regcf other hard part of the implementation I believe but almost the last!
  regcf(beta, nreq, ifault, cols, work, r, tolerances, tol_set, D, rhs, row_ptr);
  ss(cols, sserr, rss, rss_set, D, rhs);
  double var = rss[nreq] / (nobs - nreq);
  // Implement cov
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
  lsq(A, weights, y, rows, cols);
  
}
