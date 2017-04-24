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
  int row, i, jmax[1], ipt, newpos, iwk[max_size];
  double ss[last], smax[1], temp;

  ifault[0] = 0;
  if(first >= cols) ifault[0] = 1;
  if(last <= 0) ifault[0] = ifault[0]+2;
  if(first < 0) ifault[0] = ifault[0]+4;
  if(last >= cols) ifault[0] = ifault[0]+8;
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

void exadd1(int ivar, double* smax, int* jmax, double* ss, int last, int max_size, double* rss, double* bound, int nbest, double** ress, int* vorder, int** lopt) {
  int ltemp, jm;
  double ssbase, wk[last], temp, sm;
  
  if(jmax[0] == 0) return;
  if(ivar <= 0) return;
  if(ivar > max_size) return;
  ltemp = vorder[ivar];
  jm = jmax[0];
  sm = smax[0];
  if(ivar > 0) ssbase = rss[ivar-1];
  if(ivar == 0) ssbase = rss[ivar]+ss[0];
  for(int i=ivar; i<last; i++) {
    wk[i] = ss[i];
  }
  for(int i=0; i<nbest; i++) {
    temp = std::max(ssbase-sm, 0.0);
    if(temp >= bound[ivar]) break;
    vorder[ivar] = vorder[jm];
    if(jm == ivar) vorder[ivar] = ltemp;
    report(ivar, temp, max_size, bound, nbest, ress, vorder, lopt);
    if(i >= nbest) break;
    wk[jm] = 0.0;
    sm = 0.0;
    jm = 0;
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

void add1(int first, int last, double* ss, double* smax, int* jmax, int* ifault, int cols, double* D, double* rhs, double* r, double* tol, int* row_ptr) {
  double sxx[cols], sxy[cols], diag, dy, ssqx;
  int inc, pos, row, col;
  jmax[0] = 0;
  smax[0] = 0.0;
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
      if(ss[j] > smax[0]) {
	smax[0] = ss[j];
	jmax[0] = j;
      }
    } else {
      ss[j] = 0.0;
    }
  }
}

void forwrd(int first, int last, int* ifault, int cols, int max_size, double* D, double* rhs, double* r, int nbest, double* rss, double* bound, double** ress, int* vorder, int** lopt, bool* rss_set, double* sserr, int* row_ptr, double* tol) {
  int jmax[1];
  double ss[last], smax[1];

  ifault[0] = 0;
  if(first >= cols) ifault[0] = 1;
  if(last <= 1) ifault[0] = ifault[0]+2;
  if(first < 1) ifault[0] = ifault[0]+4;
  if(last > cols) ifault[0] = ifault[0]+8;
  if(ifault[0] != 0) return;

  for(int i=first; i<max_size; i++) {
    add1(i, last, ss, smax, jmax, ifault, cols, D, rhs, r, tol, row_ptr);
    if(nbest > 0) exadd1(i, smax, jmax, ss, last, max_size, rss, bound, nbest, ress, vorder, lopt);
    if(jmax[0] > i) vmove(jmax[0], i, ifault, cols, rss_set, rss, sserr, row_ptr, D, rhs, r, vorder, tol);
    for(int j=0; j<cols;j++) {
    }
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
  incr = n;
  while(true) {
    incr = incr/3;
    if(incr == 2*(incr/2)) incr++;
    for(int i=0; i<incr; i++) {
      finish = n;
      while(true) {
	i1 = i;      
	temp = l[i1];
	it = i1;
	while(true) {
	  i2 = i1+incr;
	  if(i2 > finish) {
	    if(i1 > it) l[i1] = temp;
	    finish = finish - incr;
	    break;
	  }
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
	if(finish <= incr) break;
      }
    }
    if(incr <= 1) break;
  }
}

void report(int nv, double ssq, int max_size, double* bound, int nbest, double** ress, int* vorder, int** lopt) {
  int pos1, lists[nv], lists2[nv];
  double under1 = 0.99999999, above1 = 1.00000001;

  if(nv > max_size) return;
  if(ssq >= bound[nv]) return;

  pos1 = (nv*nv+nv)/2;

  // Find rank of the new subset
  for(int rank=0; rank<nbest; rank++) {
    if(ssq < (ress[nv][rank])*above1) {
      std::copy(vorder, vorder+nv+1, lists);
      shell(lists, nv);
      if(ssq > ress[nv][rank]*under1) {
	for(int i=0; i<nv; i++) {
	  lists2[i] = lopt[nv][pos1+i];
	}
	if(same_vars(lists, lists2, nv)) break;
      }
      
      for(int j=nbest-2; j>rank; j--) {
	ress[nv][j+1] = ress[nv][j];
	for(int k=pos1; k<pos1+nv; k++) {
	  lopt[j+1][k] = lopt[j][k];
	}
      }
      ress[nv][rank] = ssq;
      for(int k=0; k<nv+1; k++) {
	lopt[rank][pos1+k] = lists[k];
      }
      bound[nv] = ress[nv][nbest-1];
      break;
    }
  }
}
