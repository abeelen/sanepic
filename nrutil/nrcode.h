#ifndef NRCODE_H_
#define NRCODE_H_


#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void dgaussj(double **a, int n, double **b, int m); //NR a verifier
void dcovsrt(double **covar, int ma, int ia[], int mfit); //NR a verifier

void dlfit(double x[], double y[], double sig[], int ndat, double a[], int ia[],
	   int ma, double **covar, double *chisq, void (*funcs)(double, double [], int)); // NR

void dcholdc(double **a, long n, double p[]); // cholesky decomposition : NR
void dcholsl(double **a, long n, double p[], double b[], double x[]); // solve cholesky linear system : NR

#endif /* NRCODE_H_ */
