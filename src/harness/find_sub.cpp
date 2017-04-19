#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "lsq.h"
#include "sub.h"
/***********************************************************************
 Exhaustive search algorithm, using leaps and bounds, applied to
   the variables in position FIRST, ..., LAST. If FIRST > 1, variables
   in positions prior to this are forced in. If LAST < cols, variables
   in positions after this are forced out.   
************************************************************************/
void xhaust(int first, int last, int* ifault, int cols, int max_size, double* D, double* tol, double* rss, double* bound, int nbest, double** ress, int* vorder, int** lopt, double* sserr, double* rhs, double* r, int* row_ptr, bool* rss_set) {
  int row, i, jmax, ipt, newpos, iwk[max_size];
  double ss[last], smax, temp;

  ifault[0] = 0;
  if(first >= cols) ifault[0] = 1;
  if(last <= 1) ifault[0] = ifault[0]+2;
  if(first < 1) ifault[0] = ifault[0]+4;
  if(last > cols) ifault[0] = ifault[0]+8;
  if(ifault[0] != 0) return;

  // Record subsets contained in the initial ordering, including check for variables which are linearly related to earlier variables.
  // This should be redundant if the user has first called SING and init_subsets, which is included in the gold_subset function.

  for(int i=first; i<max_size; i++) {
    if(D[i] <= tol[i]) {
      ifault[0] = -999;
      return;
    }
    report(i, rss[i], max_size, bound, nbest, ress, vorder, lopt);
  }

  // IWK[i] contains the upper limit for the I-th simulated for loop for
  // I = FIRST, ...., max_size-1. IPT points to the current for loop.

  for(int i=first; i<max_size; i++) {
    iwk[i] = last;
  }
  // The following the inner most loop
  // *******************************************************************************************************
  bool outer_loop = true;
  while(true) {
    if(outer_loop) {
      add1(max_size, iwk[max_size], ss, smax, jmax, ifault, cols, D, rhs, r, tol, row_ptr);
      exadd1(max_size, smax, jmax, ss, iwk[max_size], max_size, rss, bound,  nbest, ress,  vorder, lopt);
      
      // Move to next lower numbered loop which has not been exhausted.
      ipt = max_size-1;
    }
    if(ipt >= iwk[ipt]) {
      ipt--;
      if(ipt >= first) {
	outer_loop = false;
	continue;
      } else {
	return;
      }
    }
    // Lower variable from position IPT to position IWK[ipt].
    // Record any good new subsets found by the move.
	
    newpos = iwk[ipt];
    vmove(ipt, newpos, ifault, cols, rss_set, rss, sserr, row_ptr, D, rhs, r, vorder, tol);
    int last_val = std::min(max_size, newpos-1);
    for(int i=ipt; i<last_val; i++) {
      report(i, rss[i], max_size, bound, nbest, ress, vorder, lopt);
    }
    
    // Reset all ends of loops for i>=IPT.
    
    for(int i=ipt; i<max_size; i++) {
      iwk[i] = newpos-1;
    }
	
    // If residual sum of squares for all variables above position NEWPOS
    // is greater than bound[i], no better subsets of size i can be found
    // inside the current loop.
	
    temp = rss[newpos-1];
    for(int i=ipt; i<max_size; i++) {
      if(temp > bound[i]) {
	ipt = i-1;
	if(ipt < first) return;
	outer_loop = false;
      }
    }
    if(iwk[max_size] > max_size) {
      outer_loop = true;
      continue;
    }
    ipt = max_size-1;
    outer_loop = false;
  }
// *********************************************************************************************************************************
}
	


void exadd1(int ivar, int sm, int jm, double* ss, int last, int max_size, double* rss, double* bound, int nbest, double** ress, int* vorder, int** lopt) {
  int ltemp;
  double ssbase, wk[last], temp;

  if(jm == 0) return;
  if(ivar <= 0) return;
  if(ivar > max_size) return;
  ltemp = vorder[ivar];
  if(ivar > 1) ssbase = rss[ivar-1];
  if(ivar == 1) ssbase = rss[ivar]+ss[0];
  for(int i=ivar; i<last; i++) {
    wk[i] = ss[i];
  }
  for(int i=0; i<nbest; i++) {
    temp = std::max(ssbase-sm, 0.0);
    if(temp >= bound[ivar]) break;
    vorder[ivar] = vorder[jm];
    if(jm == ivar) vorder[ivar] = ltemp;
    report(ivar, temp, max_size, bound, nbest, ress, vorder, lopt);
    if(i >= nbest) return;
    wk[jm] = 0.0;
    sm = 0.0;
    jm = 0.0;
    for(int j=ivar; j<last; j++) {
      if(wk[j] > sm) {
	jm = j;
	sm = wk[j];
      }
    }
    if(jm == 0) break;
  }
  vorder[ivar] = ltemp;
}

void add1(int first, int last, double* ss, int smax, int jmax, int* ifault, int cols, double* D, double* rhs, double* r, double* tol, int* row_ptr) {
  double sxx[cols], sxy[cols], diag, dy, ssqx;
  int inc, pos, row, col;
  jmax = 0;
  smax = 0.0;
  ifault[0] = 0;

  if(first > cols) ifault[0] = 1;
  if(last < first) ifault[0] = ifault[0]+2;
  if(first < 1) ifault[0] = ifault[0]+4;
  if(last > cols) ifault[0] = ifault[0]+8;
  if(ifault[0] != 0) return;

  for(int i=first; i<last; i++) {
    sxx[i] = 0.0;
    sxy[i] = 0.0;
  }
  inc = cols-last;
  pos = row_ptr[first];
  for(int i=first; i<last; i++) {
    diag = D[i];
    dy = diag*rhs[i];
    sxx[i] = sxx[i]+diag;
    sxy[i] = sxy[i]+dy;
    for(int j=i+1; j<last; j++) {
      sxx[j] = sxx[j]+diag*pow(r[pos],2);
      sxy[j] = sxy[j]+dy*r[pos];
      pos++;
    }
    pos = pos+inc;
  }
  for(int j=first; j<last; j++) {
    ssqx = sxx[j];
    if(sqrt(ssqx) > tol[j]) {
      ss[j] = pow(sxy[j], 2)/sxx[j];
      if(ss[j] > smax) {
	smax = ss[j];
	jmax = j;
      }
    } else {
      ss[j] = 0.0;
    }
  }
}

void forwrd(int first, int last, int* ifault, int cols, int max_size, double* D, double* rhs, double* r, int nbest, double* rss, double* bound, double** ress, int* vorder, int** lopt, bool* rss_set, double* sserr, int* row_ptr, double* tol) {
  int jmax;
  double ss[last], smax;

  ifault[0] = 0;
  if(first >= cols) ifault[0] = 1;
  if(last <= 1) ifault[0] = ifault[0]+2;
  if(first < 1) ifault[0] = ifault[0]+4;
  if(last > cols) ifault[0] = ifault[0]+8;
  if(ifault[0] != 0) return;

  for(int i=first; i<max_size; i++) {
    add1(first, last, ss, smax, jmax, ifault, cols, D, rhs, r, tol, row_ptr);
    if(nbest > 0) exadd1(i, smax, jmax, ss, last, max_size, rss, bound, nbest, ress, vorder, lopt);
    
    if(jmax > i) vmove(jmax, i, ifault, cols, rss_set, rss, sserr, row_ptr, D, rhs, r, vorder, tol);		    
  }
}

bool same_vars(int* list1, int* list2, int n) {
  bool same=true;
  for(int i=0; i<n; i++) {
    if(list1[i] != list2[i]) {
      same = false;
      return same;
    }
  }
  return same;
}

void shell(int* l, int n) {
  int finish, temp, new_value, i1, i2, incr, it;
  incr = n+2;
  while(incr > 1) {
    incr = n;
    incr = incr/3;
    if(incr == 2*(incr/2)) incr++;
    for(int i=0; i<incr+1; i++) {
      finish = n;
      while(finish > incr) {
	i1 = i;      
	temp = l[i1];
	it = i1;
	for(i2=i1+incr; i2<finish+incr; i2+=incr) {
	  new_value = l[i2];
	  if(temp > new_value) {
	    l[i1] = new_value;
	    i1 = i2;
	    continue;
	  }
	  if(i1 > it) l[i1] = temp;  
	  i1 = i2;
	  temp = new_value;
	  it = i1;
	}
	if(i1 > it) l[i1] = temp;
	finish = finish-incr;
      }
    }
  }
}

void report(int nv, double ssq, int max_size, double* bound, int nbest, double** ress, int* vorder, int** lopt) {
  int pos1, lists[nv], lists2[nv];
  double under1 = 0.99999999, above1 = 1.00000001;

  if(nv > max_size) return;
  if(ssq >= bound[nv]) return;
  pos1 = (nv*(nv-1))/2;

  // Find rank of the new subset
  for(int rank=0; rank<nbest; rank++) {
    if(ssq < (ress[nv][rank])*above1) {
      for(int j=0; j<nv; j++) {
	lists[j] = vorder[j];
      }
      shell(lists, nv);
      if(ssq > ress[nv][rank]*under1) {
	for(int i=0; i<nv; i++) {
	  lists2[i] = lopt[pos1+i][nv];
	}
	if(same_vars(lists, lists2, nv)) break;
      }
      for(int j=nbest-1; j>rank; j--) {
	ress[nv][j] = ress[nv][j-1];
	for(int k=pos1; k<pos1+nv-2; k++) {
	  lopt[k][j+1] = lopt[k][j];
	}
      }
      ress[nv][rank] = ssq;
      for(int k=0; k<nv; k++) {
	lopt[k+pos1][rank] = lists[k];
	bound[nv] = ress[nv][nbest];
      }
      break;
    }
  }
}
/*
void init_subsets(int cols, int nvar_max, bool fit_const, int nbest, double* work, double* r, double* tol, bool* tol_set, double* D, int* row_ptr, double* rhs, double* sserr, bool* rss_set, double* rss, int* vorder) {
  double eps = 1e-14;
  double vlarge = 3.4028e38;
  bool lindep[cols];
  int ier, max_size, lopt_dim1, ifault[1];
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
  tolset(cols, work, r, tol, tol_set);
  sing(lindep, ifault, cols, D, tol_set, r, tol, row_ptr, rhs, sserr, work);
  
  ss(cols, sserr, rss, rss_set, D, rhs);
  for(int i=0; i<max_size; i++) {
    report(i, rss[i], max_size, bound, nbest, ress, vorder, lopt);
   }
}
*/
