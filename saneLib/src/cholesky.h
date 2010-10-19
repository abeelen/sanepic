#ifndef CHOLESKY_H_
#define CHOLESKY_H_

void cholesky(long n, double **a, double **l);
void system_triang(double **A, long n, long m, double *x, double *b, bool sup_inf);
void solve_cholesky(double **a, double *b, double **l, double *x, long n);

#endif /* CHOLESKY_H_ */
