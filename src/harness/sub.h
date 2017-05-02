#ifndef SUB_H
#define SUB_H

void xhaust(int first, int last, int* ifault, int cols, int max_size, double* D, double* tol, double* rss, double* bound, int nbest, double** ress, int* vorder, int** lopt, double* sserr, double* rhs, double* r, int* row_ptr, bool* rss_set);

void add1(int first, int last, double* ss, double* smax, int* jmax, int* ifault, int cols, double* D, double* rhs, double* r, double* tol, int* row_ptr);

void exadd1(int ivar, double* sm, int* jm, double* ss, int last, int max_size, double* rss, double* bound, int nbest, double** ress, int* vorder, int** lopt);

void forwrd(int first, int last, int* ifault, int cols, int max_size, double* D, double* rhs, double* r, int nbest, double* rss, double* bound, double** ress, int* vorder, int** lopt, bool* rss_set, double* sserr, int* row_ptr, double* tol);

bool same_vars(int* list1, int* list2, int n);

void shell(int* l, int n);

void report(int nv, double ssq, int max_size, double* bound, int nbest, double** ress, int* vorder, int** lopt);

void magma_forwrd(int m, int n, double* A, double* B, int max_size);

void dn_forwrd(int m, int n, double* A, double* B, int max_size);


#endif
