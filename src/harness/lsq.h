#ifndef LSQ_H
#define LSQ_H

void includ(double weight, double* xrow, double y, int cols, double* D, double* r, double* rhs, double* sserr);

void ss(int cols, double* sserr, double* rss, bool* rss_set, double* D, double* rhs);

void vmove(int from, int to, int* ifault, int cols, bool* rss_set, double* rss, double* sserr, int* row_ptr, double* D, double* rhs, double* r, int* vorder, double* tolerances);

void reordr(int* list, int n, int pos1, int* ifault, int cols, int* vorder, bool* rss_set, double* rss, int* row_ptr, double* D, double* r, double* sserr, double* tolerances, double* rhs);

void partial_corr(int cols, int* ifault, int in, int dimc, double* D, double* r, double* sserr, double* rhs, double* cormat, double* ycorr);

void tolset(int cols, double* work, double* r, double* tolerances, bool* tol_set);

void sing(bool* lindep, int* ifault, int cols, double* D, bool* tol_set, double* r, double* tolerances, int* row_ptr, double* rhs, double* sserr, double* work);

void regcf(double* beta, int nreq, int* ifault, int cols, double* work, double* r, double* tolerances, bool* tol_set, double* D, double* rhs, int* row_ptr);

void inv(int nreq, double* rinv, int* row_ptr, double* r);

void cov(int nreq, double var, double* covmat, int dimcov, double* sterr, int* ifault, int cols, double* D, double vsmall, int nobs, double* rss, bool* rss_set, double* rhs, int* row_ptr, double* sserr, double* r);

void hdiag(double* xrow, int nreq, double* hii, int* ifault, int cols, double* D, double* tolerances, double* r);

void bksub2(double* x, double* b, int nreq, int cols, double* r);

double varprd(double* x, int nreq, int cols, double* sserr, int nobs, double* D, double* tolerances, double* r);

#endif


  
