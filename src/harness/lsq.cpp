/********************************************************************************************************
This file accepts a .dat file or .mtx file that has a matrix of data for a linear least squares problem.
This version is a recreation of the Alan Miller Fortran implementation of lsq.f90. 
LSQ_Gold is a sequential version used to compare parallel versions against for accuracy purposes.


Accept argv[1] = file name, argv[2] = # of rows, argv[3] = # of cols
Input file should have two extra columns at the end of 
1) weight 
2) b data

Weights need to be specified, if solving an unweighted linear least squares problem use weights of 1.0.

Developed by bsauk on 4/17/17
*********************************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "CycleTimer.h"

//#define MAXVAR 1500 // Problem specific, number of variables in equation can adjust these values
//#define MAXCASES 5000 // Problem specific, number of trials

// This will need to be parallelized
void includ(double weight, double* xrow, double y, int cols, double* D, double* r, double* rhs, double* sserr) {
  double vsmall = 2.225e-307;
  int nextr = 0;
  double w = weight;
  double xk = 0.00;
  double di = 0.00;
  double cbar = 0.00;
  double sbar = 0.00;
  double xi = 0.0;
  for(int i=0; i<cols; i++) {
    if(fabs(w) < vsmall) {
      return;
    }
    xi = xrow[i];
    if(fabs(xi) < vsmall) { // Case for ill-conditioned matrices
      nextr = nextr+cols-i-1;
    } else {
      di = D[i];
      cbar = di/(di+w*xi*xi);
      sbar = w*xi/(di+w*xi*xi);
      D[i] = di+w*xi*xi;
      w = cbar*w;
      for(int k=i+1; k<cols; k++) {
	xk = xrow[k];
	xrow[k] = xk-xi*r[nextr];
	r[nextr] = cbar*r[nextr]+sbar*xk;
	nextr++;
      }
      xk = y;
      y = xk-xi*rhs[i];
      rhs[i] = cbar*rhs[i]+sbar*xk;
    }
  }
  std::cout << "w = " << w << std::endl;
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
  double vsmall = 2.225e-307;
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
      for(int j=i+2; j<cols; j++) {
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

    // Planar rotation in regular case
    
    d1new = d2+d1*pow(x,2);
    cbar = d2/d1new;
    sbar = x*d1/d1new;
    d2new = d1*cbar;
    D[i] = d1new;
    D[i+1] = d2new;
    r[m1] = sbar;
    for(int j=i+2; j<cols; j++) {
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
	  vmove(tmp, next, ifault, cols, rss_set, rss, sserr, row_ptr, D, rhs, r, vorder, tolerances); 
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
      ifault[0] = ifault[0] + 1;
      if(i<cols-1) {
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

  for(int i=nreq; i>=0; i--) {
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
  pos = nreq*(nreq-1)/2-1;
  for(int i=nreq-2; i>=0; i--) {
    start = row_ptr[i];
    for(int j=nreq-1; j>i; j--) {
      pos1 = start;
      pos2 = pos;
      total = 0.0;
      for(int k=i+1; k<j; k++) {
	pos2 = pos2+nreq-k-1;
	total = total-r[pos1]*rinv[pos2];
	pos1++;
      }
      rinv[pos] = total-r[pos1];
      pos--;
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
    var = rss[nreq-1]/(nobs-nreq);
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
      for(int k=j+1; k<nreq; k++) {
	total = total + rinv[pos1]*rinv[pos2]/D[k];
	pos1++;
	pos2++;
      }
      covmat[pos] = total*var;
      if(i==j)  {
	sterr[i] = sqrt(covmat[pos]);
      }
      pos++;
    }
    start = start+nreq-i-1;
  }
}

void hdiag(double* xrow, int nreq, double* hii, int* ifault, int cols, double* D, double* tolerances, double* r) {
  double total, wk[cols];
  int pos;
  ifault[0] = 0;
  if(nreq > cols) ifault[0] = ifault[0]+4;
  if(ifault[0] != 0) return;

  hii[0] = 0.0;
  for(int i=0; i<nreq; i++) {
    if(sqrt(D[i]) <= tolerances[i]) {
      wk[i] = 0.0;
    } else {
      pos = i-1;
      total = xrow[i];
      for(int j=0; j<i; j++) {
	total = total-wk[j]*r[pos];
	pos = pos+cols-j-1;
      }
      wk[i] = total;
      hii[0] = hii[0]+pow(total,2)/D[i];
    }
  }
}

void bksub2(double* x, double* b, int nreq, int cols, double* r) {
  int pos;
  double temp;

  for(int i=0; i<nreq; i++) {
    pos = i-1;
    temp=x[i];
    for(int j=0; j<i; j++) {
      temp = temp-r[pos]*b[j];
      pos = pos+cols-j-1;
    }
    b[i] = temp;
  }
}


double varprd(double* x, int nreq, int cols, double* sserr, int nobs, double* D, double* tolerances, double* r) {
  double fn_val, var, wk[nreq];
  int ifault;

  fn_val = 0.0;
  ifault=0;
  if(nreq < 1 || nreq > cols) ifault=ifault+4;
  if(nobs <= nreq) ifault=ifault+8;
  if(ifault != 0) {
    std::cout << "Error in function VARPRD: ifault = " << ifault << std::endl;
    return 0;
  }
  var = sserr[0]/(nobs-nreq);
  bksub2(x, wk, nreq, cols, r);
  for(int i=0; i<nreq; i++) {
    if(D[i] > tolerances[i]) fn_val = fn_val+pow(wk[i],2)/D[i];
  }
  fn_val = fn_val*var;
  return fn_val;
}

