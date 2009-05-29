#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#include "nrutil.h"
#include "nrcode.h"

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

//numrec code, changed to double
void dcovsrt(double **covar, int ma, int ia[], int mfit)
{
        int i,j,k;
        double swap;

        for (i=mfit+1;i<=ma;i++)
                for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
        k=mfit;
        for (j=ma;j>=1;j--) {
                if (ia[j]) {
                        for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
                        for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
                        k--;
                }
        }
}


void dlfit(double x[], double y[], double sig[], int ndat, double a[], int ia[],
	int ma, double **covar, double *chisq, void (*funcs)(double, double [], int))
{
	int i,j,k,l,m,mfit=0;
	double ym,wt,sum,sig2i,**beta,*afunc;

	afunc = dvector(1,ma);
	beta = dmatrix(1,ma,1,1);

	for (j=1;j<=ma;j++)
		if (ia[j]) mfit++;
	if (mfit == 0) printf("lfit: no parameters to be fitted \n");
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=mfit;k++) covar[j][k]=0.0;
		beta[j][1]=0.0;
	}

	for (i=1;i<=ndat;i++) {
		(*funcs)(x[i],afunc,ma);
		ym=y[i];
		if (mfit < ma) {
			for (j=1;j<=ma;j++)
				if (!ia[j]) ym -= a[j]*afunc[j];
		}
		sig2i=1.0/SQR(sig[i]);
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				wt=afunc[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) covar[j][++k] += wt*afunc[m];
				beta[j][1] += ym*wt;
			}
		}
	}



	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++)
			covar[k][j]=covar[j][k];
	dgaussj(covar,mfit,beta,1);

	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) a[l]=beta[++j][1];
	*chisq=0.0;

	for (i=1;i<=ndat;i++) {
		(*funcs)(x[i],afunc,ma);
		for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
		*chisq += SQR((y[i]-sum)/sig[i]);
	}

	dcovsrt(covar,ma,ia,mfit);

	free_dvector(afunc,1,ma);
	free_dmatrix(beta,1,ma,1,1);
}


// Gauss-Jordan elimination with full pivoting NR 2.1
// NOT USED
void dgaussj(double **a, int n, double **b, int m)
{
        int *indxc,*indxr,*ipiv;
        int i,icol,irow,j,k,l,ll;
        double big,dum,pivinv,swap;

        indxc=ivector(1,n);
        indxr=ivector(1,n);
        ipiv=ivector(1,n);
        for (j=1;j<=n;j++) ipiv[j]=0;
	icol=0;
	irow=0;
        for (i=1;i<=n;i++) {
	  big=0.0;
	  for (j=1;j<=n;j++)
	    if (ipiv[j] != 1)
	      for (k=1;k<=n;k++) {
		if (ipiv[k] == 0) {
		  if (fabs(a[j][k]) >= big) {
		    big=fabs(a[j][k]);
		    irow=j;
		    icol=k;
		  }
		} else if (ipiv[k] > 1) printf("gaussj: Singular Matrix-1 \n");
	      }
	  ++(ipiv[icol]);
	  if (irow != icol) {
	    for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
				 for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
						      }
	  indxr[i]=irow;
	  indxc[i]=icol;
	  if (a[icol][icol] == 0.0) printf("gaussj: Singular Matrix-2 \n");
	  pivinv=1.0/a[icol][icol];
	  a[icol][icol]=1.0;
	  for (l=1;l<=n;l++) a[icol][l] *= pivinv;
	  for (l=1;l<=m;l++) b[icol][l] *= pivinv;
	  for (ll=1;ll<=n;ll++)
	    if (ll != icol) {
	      dum=a[ll][icol];
	      a[ll][icol]=0.0;
	      for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
	      for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
	    }
        }
        for (l=n;l>=1;l--) {
	  if (indxr[l] != indxc[l])
	    for (k=1;k<=n;k++)
	      SWAP(a[k][indxr[l]],a[k][indxc[l]]);
        }
        free_ivector(ipiv,1,n);
        free_ivector(indxr,1,n);
        free_ivector(indxc,1,n);

}


void dcholdc(double **a, long n, double p[])
{

  long i,j,k;
  double sum;

  for (i=0;i<n;i++) {
    for (j=i;j<n;j++) {
      for (sum=a[i][j],k=i-1;k>=0;k--) sum -= a[i][k]*a[j][k];
      if (i == j) {
	if (sum <= 0.0) printf("Error in dcholdc: matrix is probably not positive definite\n");
	p[i]=sqrt(sum);
      } else a[j][i]=sum/p[i];
    }
  }
}


void dcholsl(double **a, long n, double p[], double b[], double x[])
{

  long i,k;
  double sum;

  for (i=0;i<n;i++) {
    for (sum=b[i],k=i-1;k>=0;k--) sum -= a[i][k]*x[k];
    x[i]=sum/p[i];
  }
  for (i=n-1;i>=0;i--) {
    for (sum=x[i],k=i+1;k<n;k++) sum -= a[k][i]*x[k];
    x[i]=sum/p[i];
  }

}


#endif /* ANSI */
