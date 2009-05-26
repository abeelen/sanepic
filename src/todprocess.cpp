#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "todprocess.h"
#include <fftw3.h>
#include <time.h>
#include <unistd.h>

#define NR_END 1
#define FREE_ARG char*

using namespace std;



//arrays init
void init1D_double(double *A, long im, long n, double val){

  long ii;

  for (ii=im;ii<n+im;ii++)
    A[ii] = val;

}

void init1D_long(long *A, long im, long n, long val){

  long ii;

  for (ii=im;ii<n+im;ii++)
    A[ii] = val;

}


void init2D_double(double **A, long im, long jm, long nx, long ny, double val){

  long ii, jj;

  for (ii=im;ii<nx+im;ii++)
    for (jj=jm;jj<ny+jm;jj++)
      A[ii][jj] = val;

}



void init2D_long(long **A, long im, long jm, long nx, long ny, long val){

  long ii, jj;

  for (ii=im;ii<nx+im;ii++)
    for (jj=jm;jj<ny+jm;jj++)
      A[ii][jj] = val;

}



// matrix init
double ** dma(int nrl, int nrh, int ncl, int nch) {

  int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double ** m;

  m = new double* [nrow+NR_END];
  m += NR_END;
  m-=nrl;

  m[nrl]= new double[(long)nrow*ncol+NR_END];
  m[nrl] += NR_END;
  m[nrl]-=ncl;

  for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  return m;

}



long ** lma(int nrl, int nrh, int ncl, int nch) {

  int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  long ** m;

  m = new long* [nrow+NR_END];
  m += NR_END;
  m-=nrl;

  m[nrl]= new long[(long)nrow*ncol+NR_END];
  m[nrl] += NR_END;
  m[nrl]-=ncl;

  for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  return m;

}


char ** charma(int nrl, int nrh, int ncl, int nch) {

  int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  char ** m;

  m = new char* [nrow+NR_END];
  m += NR_END;
  m-=nrl;

  m[nrl]= new char[(long)nrow*ncol+NR_END];
  m[nrl] += NR_END;
  m[nrl]-=ncl;

  for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  return m;

}



void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dma() */
{
  //free((FREE_ARG) (m[nrl]+ncl-NR_END));
  //free((FREE_ARG) (m+nrl-NR_END));
  delete[]  (m[nrl]+ncl-NR_END);
  delete[] (m+nrl-NR_END);

}


void free_lmatrix(long **m, long nrl, long nrh, long ncl, long nch)
/* free a matrix of long allocated by lma() */
{
  //free((FREE_ARG) (m[nrl]+ncl-NR_END));
  //free((FREE_ARG) (m+nrl-NR_END));
  delete[] (m[nrl]+ncl-NR_END);
  delete[] (m+nrl-NR_END);

}


int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v = new int[nh-nl+1+NR_END];

	return v-nl+NR_END;
}



long *lvector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	long *v;

	v = new long[nh-nl+1+NR_END];

	return v-nl+NR_END;
}


double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;
	long nrow;

	nrow = nh-nl+1+NR_END;

	v = new double[nrow];

	return v-nl+NR_END;
}


void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}


void free_lvector(long *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}


void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}




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


void  minmax(double* data, int ndata, double *min, double *max, int *posmin, int *posmax, unsigned char *flag)
{
  int k;

  *min = data[0];
  *max = data[0];
  *posmin = 0;
  *posmax = 0;
  k=0;
  if (flag != NULL){
    while (flag[k] != 0){
      *min = data[k+1];
      *max = data[k+1];
      *posmin = k+1;
      *posmax = k+1;
      k++;
    }
  }

  for(k = 1; k < ndata; k++) {
    if ((flag == NULL || flag[k] == 0) && (isnan(data[k]) == 0)){
      if(data[k] < *min){
	*min = data[k];
	*posmin = k;
      } else if(data[k] > *max){
	*max = data[k];
	*posmax = k;
      }
    }
  }
}




#define M 7
#define NSTACK 50
#define SWAP_A(a,b) itemp=(a);(a)=(b);(b)=itemp;

void iindexx(long n, long arr[], long indx[])
{
	long i,indxt,ir=n,itemp,j,k,l=1;
	long jstack=0,*istack;
	long a;

	istack=lvector(1,NSTACK);
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=1;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP_A(indx[k],indx[l+1]);
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP_A(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP_A(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[l]]) {
				SWAP_A(indx[l+1],indx[l])
			}
			i=l+1;
			j=ir;
			indxt=indx[l];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP_A(indx[i],indx[j])
			}
			indx[l]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack > NSTACK) printf("NSTACK too small in iindexx.\n");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_lvector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP_A

/* (C) Copr. 1986-92 Numerical Recipes Software 5"#. */






void dlfit(double x[], double y[], double sig[], int ndat, double a[], int ia[],
	int ma, double **covar, double *chisq, void (*funcs)(double, double [], int))
{
	int i,j,k,l,m,mfit=0;
	double ym,wt,sum,sig2i,**beta,*afunc;

	afunc = dvector(1,ma);
	beta = dma(1,ma,1,1);

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



void dlinfit(double y[], double sig[], int ndat, double a[], int ia[],
	int ma, double **covar, double *chisq, double **arrays)
{
	int i,j,k,l,m,mfit=0;
	double ym,wt,sum,sig2i,**beta;

	beta = dma(1,ma,1,1);


	for (j=1;j<=ma;j++)
		if (ia[j]) mfit++;
	if (mfit == 0) printf("lfit: no parameters to be fitted \n");
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=mfit;k++) covar[j][k]=0.0;
		beta[j][1]=0.0;
	}


	for (i=1;i<=ndat;i++) {
		ym=y[i];
		if (mfit < ma) {
			for (j=1;j<=ma;j++)
			  if (!ia[j]) ym -= a[j]*arrays[i-1][j-1];
		}
		sig2i=1.0/SQR(sig[i]);
		for (j=0,l=1;l<=ma;l++) {
		        if (ia[l]) {
				wt=arrays[i-1][l-1]*sig2i;
				for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) covar[j][++k] += wt*arrays[i-1][m-1];
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
		for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*arrays[i-1][j-1];
		*chisq += SQR((y[i]-sum)/sig[i]);
	}

	dcovsrt(covar,ma,ia,mfit);

	free_dmatrix(beta,1,ma,1,1);
}



void dmrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[],
	int ma, double **covar, double **alpha, double *chisq,
	void (*funcs)(double, double [], double *, double [], int), double *alamda)
{

	int j,k,l,m;
	static int mfit;
	static double ochisq,*atry,*beta,*da,**oneda;

	if (*alamda < 0.0) {
		atry=dvector(1,ma);
		beta=dvector(1,ma);
		da=dvector(1,ma);
		for (mfit=0,j=1;j<=ma;j++)
			if (ia[j]) mfit++;
		oneda=dma(1,mfit,1,1);
		*alamda=0.001;

		dmrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs);
		ochisq=(*chisq);
		for (j=1;j<=ma;j++) atry[j]=a[j];
	}
	for (j=0,l=1;l<=ma;l++) {
		if (ia[l]) {
			for (j++,k=0,m=1;m<=ma;m++) {
				if (ia[m]) {
					k++;
					covar[j][k]=alpha[j][k];
				}
			}
			covar[j][j]=alpha[j][j]*(1.0+(*alamda));
			oneda[j][1]=beta[j];
		}
	}

	dgaussj(covar,mfit,oneda,1);
	for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
	if (*alamda == 0.0) {
		dcovsrt(covar,ma,ia,mfit);
		free_dmatrix(oneda,1,mfit,1,1);
		free_dvector(da,1,ma);
		free_dvector(beta,1,ma);
		free_dvector(atry,1,ma);
		return;
	}
	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) atry[l]=a[l]+da[++j];
	dmrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);
	if (*chisq < ochisq) {
		*alamda *= 0.1;
		ochisq=(*chisq);
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				for (j++,k=0,m=1;m<=ma;m++) {
					if (ia[m]) {
						k++;
						alpha[j][k]=covar[j][k];
					}
				}
				beta[j]=da[j];
				a[l]=atry[l];
			}
		}
	} else {
		*alamda *= 10.0;
		*chisq=ochisq;
	}
}



void dmrqcof(double x[], double y[], double sig[], int ndata, double a[], int ia[],
	int ma, double **alpha, double beta[], double *chisq,
	void (*funcs)(double, double [], double *, double [], int))
{
	int i,j,k,l,m,mfit=0;
	double ymod,wt,sig2i,dy,*dyda, *temp;


	dyda=dvector(1,ma);
	for (j=1;j<=ma;j++)
		if (ia[j]) mfit++;
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {

                temp = dyda+1;
		(*funcs)(x[i],a+1,&ymod,temp,ma);
		dyda=temp-1; //// to match indices ...

		sig2i=1.0/(sig[i]*sig[i]);
		dy=y[i]-ymod;
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				wt=dyda[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) alpha[j][++k] += wt*dyda[m];
				beta[j] += dy*wt;
			}
		}
		*chisq += dy*dy*sig2i;
	}
	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
	free_dvector(dyda,1,ma);
}



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


#undef NRANSI


void polynomia(double x, double y[], int dma)
{
  int i;

  for(i=1;i<=dma;i++){
    y[i] = pow(x,i-1);
  }

}


void dpolyfit(double x[], double y[], int ndata, int norder, double *a)
{
  int i;
  int ma;
  double chisq;
  double *sig, *b;
  int *ia;
  double** covar;

  ma = norder+1;

  sig = new double[ndata];
  ia  = new int[ma];
  covar = dma(1,ma,1,ma);

  //initialize sig to 1
  for (i=0;i<ndata;i++){
    sig[i] = 1.0;
  }
  //set to estimate all parameters
  for (i=0;i<ma;i++){
    ia[i]  = 1;
  }

  b=a-1;
  dlfit(x-1,y-1,sig-1,ndata,b,ia-1,ma,covar,&chisq,polynomia);
  a = b+1;

  delete(sig);
  delete(ia);
  free_dmatrix(covar,1,ma,1,ma);

}


void remove_poly(double y[], int ndata, int norder, double* yout, unsigned char* flag)
{
  int i, j;
  int ndint;
  double *sx, *sy;
  double* a;

  sx = new double[ndata];
  sy = new double[ndata];
  a = new double[norder+1];

  j=0;
  for (i = 0; i < ndata; i++) {
    if(flag != NULL && flag[i]) continue;
    sx[j] = i;
    sy[j] = y[i];
    j++;
  }
  ndint = j;

  dpolyfit(sx,sy,ndint,norder,a);

  //remove best fit poly
  for (i=0;i<ndata;i++) yout[i] = y[i];
  for (i=0;i<ndata;i++)
    for (j=0;j<=norder;j++)
      yout[i] -= a[j]*pow((double)i,j);

  delete(sx);
  delete(sy);
  delete(a);

}


double dgaussfit(double y[], double *a, int *ia, double rms, int ndata, int ma, unsigned char *flag, int vfl)
{
  int i, j, iter, nfdata;
  double chisq, alamda;
  double *sig, *b, *x, *y2;
  double **covar, **alpha;

  sig = new double[ndata];
  x   = new double[ndata];
  y2  = new double[ndata];
  covar = dma(1,ma,1,ma);
  alpha = dma(1,ma,1,ma);

  j=0;
  for (i=0;i<ndata;i++){
    sig[i] = rms;   //initialize sig to constant
    if (flag != NULL && (flag[i] & vfl)) continue;
    x[j] = (double)i;
    y2[j] = y[i];
    j++;
  }
  nfdata = j;
  alamda = -1.0;

  b=a-1;
  iter=0;
  while (iter<20 && alamda != 0.0){
    dmrqmin(x-1,y2-1,sig-1,nfdata,b,ia-1,ma,covar,alpha,&chisq,gausspoly,&alamda);
    iter++;
  }
  // unalocate memory for static arrays defined by dmrqmin
  alamda = 0.0;
  dmrqmin(x-1,y2-1,sig-1,nfdata,b,ia-1,ma,covar,alpha,&chisq,gausspoly,&alamda);

  a = b+1;

  delete(sig);
  delete(x);
  delete(y2);
  free_dmatrix(covar,1,ma,1,ma);
  free_dmatrix(alpha,1,ma,1,ma);

  return chisq;

}


void gausspoly(double x, double a[], double *y, double dyda[], int na)
{
        int i;
        double fac,ex,arg;

        *y=0.0;

	arg=(x-a[0])/a[1];
	ex=exp(-arg*arg/2.0);
	fac=a[2]*ex*arg;
	*y = a[2]*ex;
	dyda[2]=ex;
	dyda[0]=fac/a[1];
	dyda[1]=fac*arg/a[1];

	if (na > 3){
	  for (i=3;i<=na-1;i++){
	    *y+=a[i]*pow(x,i-3);
	    if (i == 3) dyda[i]=1;
	    else dyda[i]=(i-3)*a[i]*pow(x,i-4);
	  }
	}

}



void dlubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum;

	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii-1;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i+1;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}


#define NRANSI
#define TODTINY 1.0e-20;

int dludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=dvector(0,n-1);
	imax=0;
	*d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) return -1;
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TODTINY;
		if (j != n-1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
	free_dvector(vv,0,n-1);
	return 0;
}
#undef TODTINY
#undef NRANSI



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




void butterworth(double y[], int ndata, double f_lp, int orderB, double *yout,
		 double *bfilter, bool apodize, int napod, bool overwrite)
{

  int ii;
  double *apodwind;

  fftw_complex  *fdata;
  fftw_plan fftplan;

  //fftw_complex fdata[ns/2+1];
  fdata = new fftw_complex[ndata/2+1];

  //apodize if asked, and define plan for fft
  if (apodize){
    apodwind = apodwindow(ndata,napod);
    if (overwrite){
      for (ii=0;ii<ndata;ii++) y[ii] = apodwind[ii] * y[ii];
      fftplan = fftw_plan_dft_r2c_1d(ndata, y, fdata, FFTW_ESTIMATE);
    } else{
      for (ii=0;ii<ndata;ii++) yout[ii] = apodwind[ii] * y[ii];
      fftplan = fftw_plan_dft_r2c_1d(ndata, yout, fdata, FFTW_ESTIMATE);
    }
    delete[] apodwind;
  } else{
    fftplan = fftw_plan_dft_r2c_1d(ndata, y, fdata, FFTW_ESTIMATE);
  }


  fftw_execute(fftplan);
  fftw_destroy_plan(fftplan);

  //filter
  for (ii=0;ii<ndata/2+1;ii++){
    bfilter[ii] = pow(double(ii)/f_lp, 2*orderB) /(1.0+pow(double(ii)/f_lp, 2*orderB));
    fdata[ii][0] = fdata[ii][0]*bfilter[ii]/ndata;
    fdata[ii][1] = fdata[ii][1]*bfilter[ii]/ndata;
  }

  if (overwrite){
    fftplan = fftw_plan_dft_c2r_1d(ndata, fdata, y, FFTW_ESTIMATE);
  }else{
    fftplan = fftw_plan_dft_c2r_1d(ndata, fdata, yout, FFTW_ESTIMATE);
  }
  fftw_execute(fftplan);
  fftw_destroy_plan(fftplan);

  delete [] fdata;

}




double* apodwindow(int ns, int nn)
{

  int ii;
  double *apodis;

  apodis = new double[ns];

  for (ii=0;ii<ns;ii++){
    apodis[ii] = 1.0;
  }

  if (nn){
    for (ii=0;ii<nn;ii++){
      apodis[ii] = (sin(double(ii)/(nn-1.0)*M_PI - M_PI/2.0) + 1.0)/2.0;
    }
    for (ii=ns-nn;ii<ns;ii++){
      apodis[ii] = (sin(double(ns-ii-1)/(nn-1.0)*M_PI - M_PI/2.0) + 1.0)/2.0;
    }
  }

  return apodis;

}




//*************************** linear Prediction ****************************//

void memcof(double data[], int n, int m, double d[])
{
  int k,j,i;
  double *wk1, *wk2, *wkm;


  wk1 = dvector(1,n+1);
  wk2 = dvector(1,n+1);
  wkm = dvector(1,m+1);



  wk1[1] = data[1];
  wk2[n - 1] = data[n];

  for (j = 2; j <= n - 1; j++) {
    wk1[j] = data[j];
    wk2[j - 1] = data[j];
  }

  for (k = 1; k <= m; k++) {
    double num = 0.0, denom = 0.0;

    for (j = 1; j <= (n - k); j++) {
      num += wk1[j] * wk2[j];
      denom += wk1[j] * wk1[j] + wk2[j] * wk2[j];
    }

    d[k]=2.0 * num / denom;
    for (i = 1; i <= (k - 1); i++)
      d[i] = wkm[i] - d[k] * wkm[k - i];


    if (k == m) {
      return;
    }

    for (i=1; i <= k; i++)
      wkm[i] = d[i];

    for (j=1; j <= (n - k - 1); j++) {
      wk1[j] -= wkm[k] * wk2[j];
      wk2[j] = wk2[j + 1] - wkm[k] * wk1[j + 1];
    }
  }

  free_dvector(wk1,1,n+1);
  free_dvector(wk2,1,n+1);
  free_dvector(wkm,1,m+1);


}


/* d = data vector */
/* window = index of first sample of data */
/* samples = # of samples of data */
/* run = length of d == samples + 2 * window */

#define SAMPLELEN 5000   //200
#define D_RANGE 100
#define NCOEF 4500            //160
void Pad(double* d, int window, int samples, int run) {

  int i, j;
  double coeff[NCOEF];

  memcof(&d[window - 1], SAMPLELEN, NCOEF, coeff - 1);

  for (j = 0; j <= window; ++j) {
    d[window - j] = 0;
    for (i = 0; i < NCOEF; ++i)
      d[window - j] += coeff[i] * d[window - j + 1 + i];
  }

  memcof(&d[window + samples - SAMPLELEN - 1], SAMPLELEN, NCOEF,
      coeff - 1);
  for (j = window + samples; j < run; ++j) {
    d[j] = 0;
    for (i = 0; i < NCOEF; ++i)
      d[j] += coeff[i] * d[j - 1 - i];
  }
}







void binnedSpectrum2log_interpol(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode)
{

  /////////////////
  //
  // ell is an array of double, units are Hz
  //
  ////////////////


  int ii, k, counttemp, f_lp;
  double ellmin, ellmax, kmin, kmax, a, b;
  double *ellm;
  double N_flp;

  // interpolate logarithmically the noise power spectrum

  ellm = new double[nbins];
  for (ii=0;ii<nbins;ii++)
    ellm[ii] = exp((log(ell[ii+1])+log(ell[ii]))/2.0);

  counttemp = 0;
  ellmin = ellm[0];
  ellmax = ellm[1];
  kmin = ellmin*ns/fsamp;
  kmax = ellmax*ns/fsamp;

  a = (log(SpN[1]) - log(SpN[0]))/(log(kmax)-log(kmin));
  b = log(SpN[0]);


  if (mode == NULL){
    for (k=0;k<ns/2+1;k++){
      while (double(k) > ellm[counttemp]*ns/fsamp && counttemp < nbins-1){
	counttemp++;
      }
      if (counttemp > 0){
	ellmin = ellm[counttemp-1];
	ellmax = ellm[counttemp];
	kmin = ellmin*ns/fsamp;
	kmax = ellmax*ns/fsamp;
	if ((abs(SpN[counttemp]) > 0) || (SpN[counttemp] > 0)){
	  a = (log(SpN[counttemp]) -
	       log(SpN[counttemp-1]))/(log(kmax)-log(kmin));
	  b = log(SpN[counttemp-1]);
	} else {
	  a = 0;
	  b = 0;
	}
      }
      if ((SpN[counttemp] > 0) || (SpN[counttemp] > 0))
	Nk[k] = exp(a*(log((double)k)-log(kmin))+b)/double(ns);
      else {
	Nk[k] = 0.0;
      }
      //apply filter
      Nk[k] *= pow(bfilter[k],2);

    }
    Nk[0] = Nk[1];



    f_lp = 0;
    while (bfilter[f_lp] <= 0.5) f_lp++;
    f_lp++;


    //give a lower limit to the spectrum
    for (k=0;k<f_lp;k++) if (Nk[k] < Nk[f_lp]) Nk[k] = Nk[f_lp];
    //for (k=0;k<f_lp;k++) Nk[k] = Nk[f_lp];

    // suppress effect of aafilter on the Noise Sp
    for (k=0;k<ns/2+1;k++) if (Nk[k] < Nk[ns/5]) Nk[k] = Nk[ns/5];




  } else {////// compute noise power spectrum for a given mode
    while (*mode > ellm[counttemp]*ns/fsamp && counttemp < nbins-1){
	counttemp++;
      }
    if (counttemp > 0){
      ellmin = ellm[counttemp-1];
      ellmax = ellm[counttemp];
      kmin = ellmin*ns/fsamp;
      kmax = ellmax*ns/fsamp;
      a = (log(SpN[counttemp]) -
	     log(SpN[counttemp-1]))/(log(kmax)-log(kmin));
      b = log(SpN[counttemp-1]);
    }
    *Nk = exp(a*(log(*mode)-log(kmin))+b)/double(ns);
    *Nk *= pow(bfilter[int(*mode)],2);


    f_lp = 0;
    while (bfilter[f_lp] <= 0.5) f_lp++;
    f_lp++;


    if (*mode < f_lp){

      counttemp = 0;
      while (f_lp > ellm[counttemp]*ns/fsamp && counttemp < nbins-1){
	counttemp++;
      }
      if (counttemp > 0){
	ellmin = ellm[counttemp-1];
	ellmax = ellm[counttemp];
	kmin = ellmin*ns/fsamp;
	kmax = ellmax*ns/fsamp;
	a = (log(SpN[counttemp]) -
	     log(SpN[counttemp-1]))/(log(kmax)-log(kmin));
	b = log(SpN[counttemp-1]);
      }
      N_flp = exp(a*(log((double)f_lp)-log(kmin))+b)/double(ns);
      N_flp *= pow(bfilter[f_lp],2);



      //give a lower limit to the spectrum
      if (*Nk < N_flp) *Nk = N_flp;

      // suppress effect of aafilter on the Noise Sp
      for (k=0;k<ns/2+1;k++) if (Nk[k] < Nk[ns/5]) Nk[k] = Nk[ns/5];
    }
  }




  delete [] ellm;

}






void InvbinnedSpectrum2log_interpol(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode)
{

  /////////////////
  //
  // ell is an array of double, units are Hz
  //
  ////////////////


  int ii, ibin, k, counttemp, f_lp;
  double ellmin, ellmax, kmin, kmax, lkmin, lkmax, a, b;
  double *ellm, *logSpN;
  double N_flp;



  /*
  counttemp = 0;

  if (mode == NULL){
    for (k=0;k<ns/2+1;k++){
      while (double(k) > ell[counttemp]*ns/fsamp && counttemp < nbins){
	counttemp++;
      }
      if (counttemp > 0){
	Nk[k] = SpN[counttemp-1]/double(ns);
      } else {
	Nk[k] = SpN[0]/double(ns);
      }
      //apply filter
      Nk[k] *= pow(bfilter[k],2);

    }
    Nk[0] = Nk[1];



    f_lp = 0;
    while (bfilter[f_lp] <= 0.5) f_lp++;
    f_lp++;

    //give a lower limit to the spectrum
    //for (k=0;k<f_lp;k++) if (Nk[k] < Nk[f_lp]) Nk[k] = Nk[f_lp];
    for (k=0;k<f_lp;k++) Nk[k] = Nk[f_lp];


  } else {////// compute noise power spectrum for a given mode
    while (*mode > ell[counttemp]*ns/fsamp && counttemp < nbins){
      counttemp++;
    }
    if (counttemp > 0){
      *Nk = SpN[counttemp-1]/double(ns);
    } else {
      *Nk = SpN[0]/double(ns);
    }

    f_lp = 0;
    while (bfilter[f_lp] <= 0.5) f_lp++;
    f_lp++;

    if (*mode < f_lp){

      counttemp = 0;
      while (f_lp > ell[counttemp]*ns/fsamp && counttemp < nbins){
	counttemp++;
      }
      if (counttemp > 0){
	N_flp = SpN[counttemp-1]/double(ns);
      } else {
	N_flp = SpN[0]/double(ns);
      }
      N_flp *= pow(bfilter[f_lp],2);



      //give a lower limit to the spectrum
      if (*Nk < N_flp) *Nk = N_flp;

    }
  }

  */





  ellm = new double[nbins];
  logSpN = new double[nbins];


  for (ii=0;ii<nbins;ii++)
    ellm[ii] = exp((log(ell[ii+1])+log(ell[ii]))/2.0);


  counttemp = 0;
  ellmin = ellm[0];
  ellmax = ellm[1];
  kmin = ellmin*ns/fsamp;
  kmax = ellmax*ns/fsamp;
  lkmin = log(kmin);
  lkmax = log(kmax);


  if (mode == NULL){



    for (k=1;k<=(long)kmin;k++){

      Nk[k] = SpN[0]/double(ns);

    }




    if (0){

      for (ii=0;ii<nbins;ii++)
	logSpN[ii] = log(abs(SpN[ii]));
      //////////////////////////////////////////  log-log interpolation
      for (ibin=0;ibin<nbins-1;ibin++){
	ellmin = ellm[ibin];
	ellmax = ellm[ibin+1];
	kmin = ellmin*ns/fsamp;
	kmax = ellmax*ns/fsamp;
	lkmin = log(kmin);
	lkmax = log(kmax);
	if (abs(SpN[ibin]) > 0){
	  a = (logSpN[ibin+1] - logSpN[ibin])/(lkmax-lkmin);
	  b = logSpN[ibin];
	} else {
	  a = 0;
	  b = 0;
	}
	for (k=(long)kmin;k<long(kmax);k++){
	  if (SpN[ibin] > 0)
	    Nk[k] = exp(a*(log((double)k)-lkmin)+b)/double(ns);
	  else {
	    if (SpN[ibin] < 0){
	      Nk[k] = -exp(a*(log((double)k)-lkmin)+b)/double(ns);
	    } else {
	      Nk[k] = 0.0;
	    }
	  }
	  //apply filter
	  Nk[k] *= pow(bfilter[k],2);

	}
      }
    }





    /////////////  linear interpolation
    for (ibin=0;ibin<nbins-1;ibin++){
      ellmin = ellm[ibin];
      ellmax = ellm[ibin+1];
      kmin = ellmin*ns/fsamp;
      kmax = ellmax*ns/fsamp;
      if (abs(SpN[ibin]) > 0){
	a = (SpN[ibin+1] - SpN[ibin])/(kmax-kmin)/double(ns);
	b = SpN[ibin]/double(ns);
      } else {
	a = 0.0;
	b = 0.0;
      }
      for (k=long(kmin+1);k<=long(kmax);k++){
	Nk[k] = (a*((double)k-kmin)+b);

	//apply filter
	Nk[k] *= pow(bfilter[k],2);

      }
    }




    for (k=long(kmax);k<ns/2+1;k++){

      Nk[k] = SpN[nbins-1]*pow(bfilter[k],2)/double(ns);

    }


    Nk[0] = Nk[1];








    f_lp = 0;
    while (bfilter[f_lp] > 2.0) f_lp++;
    f_lp++;



    //give a lower limit to the spectrum
    //for (k=0;k<f_lp;k++) if (Nk[k] < Nk[f_lp]) Nk[k] = Nk[f_lp];
    for (k=0;k<f_lp;k++) Nk[k] = Nk[f_lp];

    // suppress effect of aafilter on the Noise Sp
    for (k=ns/20;k<ns/2+1;k++) Nk[k] = Nk[ns/20];














  } else {////// compute noise power spectrum for a given mode
    while (*mode > ellm[counttemp]*ns/fsamp && counttemp < nbins-1){
	counttemp++;
      }
    if (counttemp > 0){
      ellmin = ellm[counttemp-1];
      ellmax = ellm[counttemp];
      kmin = ellmin*ns/fsamp;
      kmax = ellmax*ns/fsamp;
      a = (log(SpN[counttemp]) -
	     log(SpN[counttemp-1]))/(log(kmax)-log(kmin));
      b = log(SpN[counttemp-1]);
    }
    *Nk = exp(a*(log(*mode)-log(kmin))+b)/double(ns);
    *Nk *= pow(bfilter[int(*mode)],2);


    f_lp = 0;
    while (bfilter[f_lp] > 2.0) f_lp++;
    f_lp++;


    if (*mode < f_lp){

      counttemp = 0;
      while (f_lp > ellm[counttemp]*ns/fsamp && counttemp < nbins-1){
	counttemp++;
      }
      if (counttemp > 0){
	ellmin = ellm[counttemp-1];
	ellmax = ellm[counttemp];
	kmin = ellmin*ns/fsamp;
	kmax = ellmax*ns/fsamp;
	a = (log(SpN[counttemp]) -
	     log(SpN[counttemp-1]))/(log(kmax)-log(kmin));
	b = log(SpN[counttemp-1]);
      }
      N_flp = exp(a*(log((double)f_lp)-log(kmin))+b)/double(ns);
      N_flp *= pow(bfilter[f_lp],2);



      //give a lower limit to the spectrum
      if (*Nk < N_flp) *Nk = N_flp;

      // suppress effect of aafilter on the Noise Sp
      //for (k=0;k<ns/2+1;k++) if (Nk[k] < Nk[ns/100]) Nk[k] = Nk[ns/100];
    }
  }





  delete [] ellm;
  delete [] logSpN;

}










/*
void deconv_antialias(double y[], int ndata, double f_lp, double* yout, bool apodize)
{

  int k, j;
    double omega, filter, ofilter;

  double omega_nyquist = 2.0*M_PI*0.5*HUNDRED_HZ_CLK;
  double omega_lp = 2.0*M_PI*f_lp;
  double omega_max = omega_lp + 0.2*omega_nyquist;

  const double stage[] = {208.0/ADC_CLK, 175.0/ADC_CLK,
			 147.0/ADC_CLK, 123.0/ADC_CLK};

  Vector v(0, ndata-1, 0.0);

  for(j = 0; j <= ndata-1; j++) {
    v[j] = y[j];

    // apodize with a Hanning window function
    if(apodize) v[j] *= 1.0 - cos(2.0*M_PI * (double)j/(double)ndata);
  }

  FFT(v);

  // Filter here
  for (ofilter = 1.0,  k = 0; k <= ndata/2 ; k++) {
    omega = 2.0*(double)k/double(ndata) * omega_nyquist;

    if (omega < omega_lp) {
      for(filter = 1.0, j = 0; j < 4; j++) {
	if (omega != 0) filter /= sin(0.5*omega*stage[j])/(0.5*omega*stage[j]);
      }
      ofilter = filter;


    } else if (omega < omega_max){
      // Low pass filter
      filter = ofilter *
	(0.5 + 0.5*cos(M_PI*(omega - omega_lp)/(omega_max - omega_lp)));
    } else {
      filter = 0.0;
    }


    // negative frequency treatment
    if(k == 0 || k == ndata/2) {
      //cout << omega/(2*M_PI) << " " << fabs(v[k]) << " ";
      v[k] *= filter;
      //cout <<  fabs(v[k]);
    } else {
      //cout << omega/(2*M_PI) << " " << sqrt(SQR(v[k]) + SQR(v[out_size - k])) << " ";
      v[k] *= filter;
      v[ndata - k] *= filter;
      //cout << sqrt(SQR(v[k]) + SQR(v[out_size - k]));
    }
    //cout <<  " " << filter << endl;
  }

  InverseFFT(v);
  v /= (double)ndata;


  for(j=0 ; j <= ndata-1; j++) {
    if(apodize) v[j] /= 1.0 - cos(2.0*M_PI * (double)j/(double)ndata) + 1E-20;
    yout[j] = v[j];
  }
}
*/





int compare( const void *arg1, const void *arg2 )
{

  double diff;
  diff = (double*)arg1 - (double*)arg2;
  return (diff < 0) ? -1 : (diff == 0) ? 0 : 1;

}

int compare_doubles (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}

int compare_long (const void *a, const void *b)
{
  const long *da = (const long *) a;
  const long *db = (const long *) b;

  return (*da > *db) - (*da < *db);
}




void filter_despike(double y[], int ndata, double* yout, unsigned char *flag)
{

  int ii, j, k, l;

  long nzone = 10000;

  double *apodwind, *youttemp, *youttemp2, *youtsort, *yfilt;

  //double TAU=3.5;
  //double S_SPK=4.0;
  //int L_SPK = ((int)(4*TAU));
  //double GAIN=6.0;
  double critring = 5.0;
  int nc = 7;

  //double s = S_SPK*(2.0*M_PI);
  //double nu = 1./TAU;
  double hsig, hmu;


  youtsort = new double[ndata];

  apodwind = apodwindow(ndata,100000);

  for(j=0;j<ndata;j++){
    yout[j] = y[j];
    // apodize
    yout[j] *= apodwind[j];
  }

  delete[] apodwind;







  /*

  //compute fft

  fftw_complex  *fdata;
  fftw_plan fftplan;

  //fftw_complex fdata[ns/2+1];
  fdata = new fftw_complex[ndata/2+1];


  fftplan = fftw_plan_dft_r2c_1d(ndata, y, fdata, FFTW_ESTIMATE);
  fftw_execute(fftplan);
  fftw_destroy_plan(fftplan);

  for ( k = 0; k <= ndata/2 ; k++) {
    f = (double)k/double(ndata) - nu;
    filter = exp(-0.5*f*f*s*s);
    fdata[k][0] *= filter;
    fdata[k][1] *= filter;
  }

  //Inverse fft
  fftplan = fftw_plan_dft_c2r_1d(ndata, fdata, yout, FFTW_ESTIMATE);
  fftw_execute(fftplan);
  fftw_destroy_plan(fftplan);
  delete[] fdata;


  i_s = (int)s;
  for(k=0;k<i_s;k++) yout[k] = 0.0;



  // loop over some intervals in data
  ypartsort = new double[nzone];

  xzone = 0;
  while (xzone < (ndata-nzone)){

    for (j=0;j<nzone;j++)
      youtsort[j] = yout[j+xzone];

    qsort(youtsort,ndata,sizeof(double),compare_doubles);
    hmean  = youtsort[ndata/2];
    hsigma = youtsort[int((double)ndata*0.84)]-hmean;
    hvariance= pow(hsigma,2);




  //despike
  for(k=0;k<nzone;k++)
    if(fabs(yout[k])>GAIN*hsigma)


      for (j=-5;j<30;j++)
	if (k+j >= 0 && k+j < ndata)
	  flag[k+j] |= 1;

  xzone += nzone-1000;

  }

   */



  //for (j=0;j<ndata;j++)
  //  yout[j] = youtsort[j];





  printf("ok\n");



  // produce a timeline taking differencies
  yfilt = new double[ndata];
  for (j=0;j<ndata;j++)
    yfilt[j] = 0.0;

  for (j=2;j<ndata-2;j++)
    yfilt[j] = y[j] - 0.8*y[j-1] - 0.8*y[j+1] + 0.3*y[j-2] + 0.3*y[j+2];




  youttemp = new double[nzone];
  youttemp2 = new double[nzone];

  for (ii=0;ii<ndata/nzone*2-2;ii++){

    k = ii*nzone/2;

    cutdata(yfilt,k,k+nzone-1,youttemp);
    for (j=0;j<nzone;j++)
      youttemp2[j] = youttemp[j];
    qsort(youttemp2,nzone,sizeof(double),compare_doubles);

    hmu  = youttemp2[nzone/2];
    hsig = youttemp2[int((double)nzone*0.84)]-hmu;

    for (j=0;j<nzone;j++)
      if (fabs(youttemp[j]) > critring*hsig)
	for (l=-nc;l<nc;l++)
	  if ((j+k+l > 0) && (j+k+l < ndata))
	    flag[j+k+l] |= 2;

  }


  for (j=0;j<ndata;j++)
    yout[j] = yfilt[j];


  delete [] youtsort;
  delete [] yfilt;
  delete [] youttemp;
  delete [] youttemp2;

}












double* CR_tf(int ndata)
{

  int i;

  const double stage[] = {208.0/ADC_CLK, 175.0/ADC_CLK,
			 147.0/ADC_CLK, 123.0/ADC_CLK};

  int nn = ndata*200;

  double aa;
  double *data, *ldata, *f1, *f2, *f3, *f4;

  data = new double[nn];
  ldata = new double[ndata];
  f1 = new double[nn];
  f2 = new double[nn];
  f3 = new double[nn];
  f4 = new double[nn];

  //simulate a spike
  for (i=0;i<nn;i++)
    data[i] = 0.0;
  data[nn/2-3*200] = 1.0;
  for (i=0;i<nn;i++){
    f1[i] = 0.0;
    f2[i] = 0.0;
    f3[i] = 0.0;
    f4[i] = 0.0;
  }

  for (i=int(stage[0]*20000.0);i<nn;i++){
    f1[i] = f1[i-1] + data[i] - data[i-int(stage[0]*20000.0)];
    f2[i] = f2[i-1] + f1[i] - f1[i-int(stage[1]*20000.0)];
    f3[i] = f3[i-1] + f2[i] - f2[i-int(stage[2]*20000.0)];
    f4[i] = f4[i-1] + f3[i] - f3[i-int(stage[3]*20000.0)];
  }

  delete[] f1;
  delete[] f2;
  delete[] f3;
  delete[] data;

  aa = f4[0];
  for (i=1;i<nn;i++)
    if (aa < f4[i])
      aa = f4[i];

  for (i=0;i<ndata;i++){
    ldata[i] = f4[i*200]/aa;
  }

  delete[] f4;

  return ldata;
}





void cutdata(double y[], int indm, int indp, double *yout)
{
  int i;

  for (i=indm;i<=indp;i++){
    yout[i-indm] = y[i];
  }
}


void cutdata(unsigned char y[], int indm, int indp, unsigned char *yout)
{
  int i;

  for (i=indm;i<=indp;i++){
    yout[i-indm] = y[i];
  }
}


void mergedata(double y1[], int ndata1, double y2[], int ndata2, double *yout)
{
  int i;

  for (i=0;i<ndata1;i++)
    yout[i] = y1[i];
  for (i=0;i<ndata2;i++)
    yout[i+ndata1] = y2[i];

}


void dindgen(int nn, double *y)
{
  int i;

  for (i=0;i<nn;i++)
    y[i]=(double)i;

}


void sort(double y[], int nn, double *yout, int *nrel, unsigned char *flag)
{

  int i, j;
  double temp;

  if (flag != NULL){
    *nrel=0;
    for (i=0;i<nn;i++){
      if (flag[i] == 0){
	yout[*nrel] = y[i];
	*nrel = *nrel+1;
      }
    }
  }
  else{
    *nrel=nn;
    for (i=0;i<nn;i++) yout[i]=y[i];
  }


  for (i=0;i<*nrel;i++){
    for (j=1;j<*nrel;j++){
      if (yout[j-1] > yout[j]){
	temp = yout[j];
	yout[j] = yout[j-1];
	yout[j-1] = temp;
      }
    }
    //printf("inside sort %d\n",i);
  }

}



void findspike(double y[], int ndata, double transfer[], const int ntr, double* yout, unsigned char* flag, unsigned char* pflag, double** pplanet, int *xerror, int *pos, int *xcosm, int *xspike, double *chisqs, double *amplspike, double *allhsig, int *posp, double *amplp, double *allhsigp)
{
  // data are assumed to vary linearly in every window

  int i, iter, xorig, count, pl_count, ierr, countspike, icosm, countp;
  int posmin, posmax, posmm, nrel;
  int temp;

  int nterm = 6;
  int npb = 400;
  double sigbeam = 4.0;
  //double criterion = 3.5;
  double criterion = 4.0;
  //double criterion = 4.5;
  double critspike = 8.0;

  //// interval to cut when an event is detected
  int margcut_m = 5;
  int margcut_p = 30;
  int margcutpos_m = 5;
  int margcutpos_p = 5;
  int critsource = 10;

  int nmaxiter=10;  //// maximum number of iteration in every window
  int nwindowsig=50;

  bool balt = 0;

  double intervcut;
  double valmin, valmax, hmean, hsig, chisq, tempvar, hsig_;
  double *hsigall, *strsig;

  double *y2, *smdata, *tempdata, *strval, *smsmdata;
  double *params, *a, *paramscr, *tempsig, *temppar;
  double **tempcovar, **arrfilter;
  int *wparams, *wpcr;
  unsigned char *smflag;

  hsigall = new double[nwindowsig+1];
  strsig = new double[nwindowsig+1];

  y2       = new double[ndata];
  smdata   = new double[npb];
  tempdata = new double[npb];
  strval   = new double[npb];
  smflag   = new unsigned char[npb];

  params = new double[nterm];
  wparams = new int[nterm];
  a = new double[2];

  if (ntr){
    tempsig = new double[ntr];
    tempcovar = dma(1,2,1,2);
    arrfilter = dma(0,ntr-1,0,1);
    wpcr = new int[2];
    paramscr = new double[2];
    smsmdata = new double[ntr];
  }

  xorig=0;
  count=0;
  pl_count=0;
  ierr=0;
  countspike=0;
  icosm=0;
  countp=0;
  posmm=0;

  ////copy data
  for (i=0;i<ndata;i++){
    y2[i]=y[i];
  }

  iter=0;

  for (i=0;i<nwindowsig;i++)
    hsigall[i] = 1.16993e-05;

  while ((xorig+npb+npb/2) < ndata){

    cutdata(y2,xorig,xorig+npb-1,tempdata);
    cutdata(flag,xorig,xorig+npb-1,smflag);
    remove_poly(tempdata,npb,1,smdata,smflag);

    minmax(smdata,npb,&valmin,&valmax,&posmin,&posmax,smflag);

    if (abs(valmin) > abs(valmax)) posmm = posmin;
    if (abs(valmin) <= abs(valmax)) posmm = posmax;

    //center on the extremum if well located in the window
    if ((xorig > npb/2-posmm) && (posmm < npb-int(4.0*sigbeam))){

      xorig = xorig-npb/2+posmm;

      cutdata(y2,xorig,xorig+npb-1,tempdata);
      cutdata(flag,xorig,xorig+npb-1,smflag);
      remove_poly(tempdata,npb,1,smdata,smflag);

      // estimate mean and standard deviation of the data in the
      // small interval (may be improved)

      sort(smdata,npb,strval,&nrel,smflag);

      hmean  = strval[nrel/2];
      //hsig = MIN(strval[int((double)nrel*0.84)]-hmean,hmean-strval[int((double)nrel*0.16)]); //it is slightly biased, and not optimal, but who cares?
      hsig = strval[int((double)nrel*0.84)]-hmean;


      hsig_ = hsig;
      hsigall[nwindowsig] = hsig;
      sort(hsigall,nwindowsig+1,strsig,&temp);
      hsig = strsig[(nwindowsig+1)/2];

      //// 1st test: spike larger than threshold
      if (abs(smdata[npb/2]-hmean) > criterion*hsig && iter < nmaxiter){


	  //// identify single corrupted pixels
	if (abs(smdata[npb/2]-hmean) > critspike*hsig && abs(smdata[npb/2+1]-hmean) < 4.0*hsig && abs(smdata[npb/2-1]-hmean) <4.0*hsig){
	    xspike[countspike] = xorig+npb/2;
	    countspike++;

	    //yout[xorig+npb/2] = (yout[xorig+npb/2+1]+yout[xorig+npb/2-1])/2.0;
	    if ((flag[xorig+npb/2] & 1) == 0){
	      flag[xorig+npb/2] += 1;
	      pflag[xorig+npb/2] += 1;
	    }
	    if ((flag[xorig+npb/2] & 2) == 0){
	      flag[xorig+npb/2] += 2;
	      pflag[xorig+npb/2] += 2;
	    }
	  }


	    ////////// if group of pixels ////////////
	else{

	    // first guess of parameters for Gaussian fit
	    params[0]=double(npb/2);
	    params[1]=sigbeam;
	    params[2]=(smdata[npb/2]-hmean);///abs(smdata[npb/2]-hmean)*hsig;
	    params[3]=hmean;
	    params[4]=hsig;
	    params[5]=hsig;
	    for (i=0;i<nterm;i++)
	      wparams[i]=1;


	    chisq = dgaussfit(smdata,params,wparams,hsig,npb,4,smflag,1);

	      //fp = fopen("/home/patanch/test.txt","w");
	      //for (i=0;i<npb;i++)
	      //fprintf(fp,"%15.10g \n",smdata[i]);
	      //fprintf(fp," \n");
	      //fclose(fp);

	    printf("xorig = %d \r",xorig);
	    if (count % 300 == 0) printf(" \n");


	    /// remove gaussian if relevant
	    if (params[1] > sigbeam  && params[1] < 5.0*sigbeam  && abs(params[2]) < critsource*hsig){
	      for (i=xorig+int(params[0]-6.0*sigbeam);i<xorig+(int(params[0])+6.0*sigbeam);i++)
		y2[i] = y2[i]-params[2]*exp(-double(i-params[0]-xorig)*double(i-params[0]-xorig)/2.0/params[1]/params[1]);
	    }



	    /// detect fitting error: if maximum location is far from the center of the gaussian
	    //or if maximum amplitude is very different than gaussian amplitude,
	    //and if center is not flagged
	    if ( (  abs(params[0]-double(npb/2)) > 4.0*sigbeam  | (params[2]+hmean-smdata[npb/2] > 4.0*hsig)) && !flag[npb/2]){
	      xerror[ierr] = xorig+npb/2;
	      ierr++;
	      for (i=xorig+npb/2-margcut_m;i<xorig+npb/2+margcut_p;i++)
		flag[i] = flag[i] | 64;
	      pflag[xorig+npb/2] = pflag[xorig+npb/2] | 64;


	      /// should recompute minmax within 4.0*sigbeam

	      //// redo the fit fixing the location and the amplitude of the event
	      params[0]=double(npb/2);
	      params[1]=sigbeam/10.0;
	      params[2]=smdata[npb/2]-hmean;
	      wparams[0]=0;
	      wparams[2]=0;

	      chisq = dgaussfit(smdata,params,wparams,hsig,npb,4,smflag,1);

	    }


	    //////////////////////////////  fit CR transfer function if set
	    if (ntr) {
	      //initialize sig to 1
	      for (i=0;i<ntr;i++){
		tempsig[i] = 1.0;
	      }
	      //set to estimate all parameters
	      for (i=0;i<2;i++){
		wpcr[i]  = 1;
	      }
	      cutdata(smdata,npb/2-ntr/2,npb/2+ntr/2-1,smsmdata);

	      for (i=0;i<ntr;i++){
		arrfilter[i][0] = 1.0;
		arrfilter[i][1] = transfer[i];
	      }

	      temppar = paramscr-1;
	      dlinfit(smsmdata-1,tempsig-1,ntr,temppar,wpcr-1,2,tempcovar,&tempvar,arrfilter);
	      paramscr = temppar+1;
	    }
	    /////////////////////////////



	    //// 2nd test: CR transfer funct matches data AND no strong planet
	    //// if !nrt sigma Gauss should be smaller than sigbeam

	    if (((ntr && abs(paramscr[1]) > criterion*hsig  && !(params[1] > sigbeam && -params[2] > criterion*hsig)) || !ntr && params[1] < sigbeam) && abs(params[0]-double(npb/2)) < 4.0*sigbeam){

	      pos[count] = xorig+int(params[0]);
	      chisqs[count] = chisq;
	      count++;

	      //// 3rd test (for cosmic rays detection): negative amplitude
	      if ((ntr && paramscr[1] < 0) | (!ntr && params[2] < 0)){

		xcosm[icosm] = xorig+int(params[0]);
		amplspike[icosm] = smdata[npb/2]-hmean;
		allhsig[icosm] = hsig;
		icosm++;
		for (i=xorig+npb/2-margcut_m;i<xorig+npb/2+margcut_p;i++)
		  flag[i] = flag[i] | 5;
		pflag[xorig+npb/2] = pflag[xorig+npb/2] | 5;

	      }
     	      else{
		/// positive for some reason

		posp[countp] = xorig+int(params[0]);
		amplp[countp] = smdata[npb/2]-hmean;
		allhsigp[countp] = hsig;
		countp++;
		for (i=xorig+npb/2-margcutpos_m;i<xorig+npb/2+margcutpos_p;i++){
		  if ((flag[i] & 16) == 0)
		    flag[i] += 16;
		  /// flag large positive data point
		  if (abs(smdata[npb/2]-hmean) > critspike*hsig && (flag[i] & 1) == 0)
		    flag[i] += 1;
		}
		pflag[xorig+npb/2] = pflag[xorig+npb/2] | 16;

		/// flag large positive data point
		if (abs(smdata[npb/2]-hmean)  > critspike*hsig && (pflag[xorig+npb/2] & 1) == 0)
		    pflag[xorig+npb/2] += 1;
   	      }
	    }



	    // temporary bit to mask the event for next iteration
	    if (ntr && abs(paramscr[1]) < criterion*hsig){
	      for (i=xorig+int(npb/2-sigbeam);i<xorig+int(npb/2+sigbeam);i++)
		if ((flag[i] & 128) == 0)
		  flag[i] += 128;
	    }


	    ////3rd test: source??
	    if ((params[1] > sigbeam | abs(params[0]-double(npb/2)) > 4.0*sigbeam) && abs(params[2]) > critsource*hsig){
	      //// planets??
	      //if (params[2] < 0){
	      for (i=1;i<nterm;i++)
		pplanet[pl_count][i] = params[i];
	      pplanet[pl_count][0] = params[0]+double(xorig);
	      pl_count++;

	      //if (params[2] <= 25) intervcut = 2.0*params[2];
	      if (1 || params[2] > 25 ) intervcut = 50.0;
	      for (i=xorig+int(params[0]-intervcut);i<xorig+(int(params[0])+intervcut);i++){
		//y2[i] = y2[i]-params[2]*exp(-double(i-params[0]-xorig)*double(i-params[0]-xorig)/2.0/params[1]/params[1]);
		if ((flag[i] & 32) == 0)
		  flag[i] += 32;
	      }
	      pflag[xorig+int(params[0])] = pflag[xorig+int(params[0])] | 32;
	      //}
	    }


	}
	xorig += npb/2 - posmm;
	iter++;
	if (balt) xorig -= npb/2-int(4.0*sigbeam);
	balt = 0;
      }
      else{
	xorig += npb - posmm;
	iter=0;
	if (balt) xorig -= npb/2-int(4.0*sigbeam);
	balt = 0;

	for (i=0;i<nwindowsig-1;i++)
	  hsigall[i] = hsigall[i+1];
	hsigall[nwindowsig-1] = hsig_;
      }
    }
    else{
      xorig += npb/2-int(4.0*sigbeam);
      balt = 1;
    }
  }

  pplanet[pl_count][0] = -1;
  pos[count] = -1;
  xerror[ierr]=-1;
  xcosm[icosm]=-1;
  xspike[countspike]=-1;
  chisqs[count]=-1;
  amplspike[icosm]=-1;
  allhsig[icosm]=-1;
  posp[countp]=-1;
  amplp[countp]=-1;
  allhsigp[countp]=-1;

  delete(hsigall);
  delete(strsig);

  delete(y2);
  delete(smdata);
  delete(tempdata);
  delete(strval);
  delete(smflag);

  delete(params);
  delete(wparams);
  delete(a);

  if (ntr){
    delete(tempsig);
    free_dmatrix(tempcovar,1,2,1,2);
    free_dmatrix(arrfilter,0,ntr-1,0,1);
    delete(wpcr);
    delete(paramscr);
    delete(smsmdata);
  }

}



void fillgaps(double y[], int ndata, double* yout, unsigned char* flag, double sign)
{
  // data are assumed to vary linearly in every window

  int i, j, count, countp, countm, seedpass;
  bool sp;

  seedpass = 0;

  //// interval to cut when an event is detected
  //int margcut_m = 15;
  //int margcut_p = 30;
  int margfit   = 20;

  double *xx, *yy, *seriep, *seriem, *tempdata1, *tempdata2, *xx2;
  double *a;
  double *valtemp;

  a = new double[2];


  //init random generator
  valtemp = randg(1,0);


  ////copy data
  for (i=0;i<ndata;i++)
    yout[i] = y[i];

  count = 0;
  sp = 0;
  countm = 0;
  while ((countm<margfit) && (flag[countm] & 0))
    countm++;


  for (i=0;i<ndata;i++){
    if (flag[i] & 1){
      count++;
      sp = 0;
    }
    else{
      sp = 1;
    }

    if (sp && count){
      countp = 0;
      while ((countp < margfit) && (countp+i<ndata-1) && (flag[i+countp] & 1) == 0){
	countp++;
	if (i+countp >= ndata) printf("SDHFIDF\n");
      }

      xx = new double[countp+countm];
      yy = new double[countp+countm];
      xx2 = new double[count];

      if (countm > 0){
	seriem = new double[countm];
	tempdata1 = new double[countm];
	dindgen(countm,seriem);
	cutdata(y,i-count-countm,i-count-1,tempdata1);
      }
      if (countp > 0){
	seriep = new double[countp];
	tempdata2 = new double[countp];
	dindgen(countp,seriep);
	cutdata(y,i,i+countp-1,tempdata2);
      }



      if (countm && countp){
	mergedata(seriem,countm,seriep,countp,xx);
	mergedata(tempdata1,countm,tempdata2,countp,yy);
      } else {
	if (countm){
	  for (j=0;j<countm;j++){
	    xx[j] = seriem[j];
	    yy[j] = tempdata1[j];
	  }
	}
	if (countp){
	  for (j=0;j<countp;j++){
	    xx[j] = seriep[j];
	    yy[j] = tempdata2[j];
	  }
	}
      }


      if (countp){
	for (j=0;j<countp;j++){
	  xx[countm+j] += double(countm+count);
	}
      }

      dpolyfit(xx,yy,countp+countm,1,a);


      dindgen(count,xx2);
      for (j=0;j<count;j++){
	xx2[j] += double(countm);
      }
      for (j=0;j<count;j++){
	valtemp = randg(1,-1);
	yout[i+j-count] = a[0]+a[1]*xx2[j] + sign*valtemp[0];
	delete [] valtemp;
      }


      if (countp){
	delete[] seriep;
	delete[] tempdata2;
      }
      if (countm){
	delete[] seriem;
	delete[] tempdata1;
      }
      delete[] xx;
      delete[] yy;
      delete[] xx2;

      countm = countp;

      countp = 0;
      count = 0;
      sp = 0;

    }
  }

  delete[] a;

}



foffset* read_mapoffsets(string fname, float *scoffsets, int *nfoff)
{
  char buffer[256];
  float p, y;
  string line, word;
  int f, s0, s1, i, fcount;
  foffset *foffsets;

  ifstream FILE (fname.c_str());
  if (! FILE.is_open()) {
    cerr << "Error opening bolometer offset file '" << fname << "'.\n";
    exit(1);
  }

  // get overall offsets
  while (! FILE.eof()) {
    FILE.getline(buffer, 255);
    line = buffer;

    s0 = 20;
    if (line.substr(0,s0) != "StarCamera2BoreSight") continue;

    // extract 6 fields from line
    i = 0;
    while (i < 6) {
      // find beginning of word
      s0 = line.find_first_not_of(" \t", s0);

      // find end of word
      s1 = line.find_first_of(" \t", s0);

      // get and storeword
      word = line.substr(s0, s1-s0);
      scoffsets[i++] = atof(word.c_str());

      // shift placeholder
      s0 = s1;
    }

    break;
  }

  // find "Begin" tag
  while (! FILE.eof()) {
    FILE.getline(buffer, 255);
    line = buffer;

    if (line.substr(0,6) == "Begin:") break;
  }

  // count frame lines
  fcount = 0;
  while (! FILE.eof()) {
    FILE.getline(buffer, 255);
    line = buffer;

    if (line[0] == '#') continue;
    if (line.substr(0,4) == "End:") break;

    fcount++;
  }

  // allocate memory
  foffsets = new foffset [fcount];
  *nfoff = fcount;

  // reset pointer
  FILE.seekg(0);

  // find "Begin" tag
  while (! FILE.eof()) {
    FILE.getline(buffer, 255);
    line = buffer;

    if (line.substr(0,6) == "Begin:") break;
  }

  // store data
  fcount = 0;
  while (! FILE.eof()) {
    FILE.getline(buffer, 255);
    line = buffer;

    if (line[0] == '#') continue;
    if (line.substr(0,4) == "End:") break;

    sscanf(line.c_str(), "%d%f%f",  &f, &p, &y);

    (foffsets[fcount]).frame = f;
    (foffsets[fcount]).pitch = p;
    (foffsets[fcount]).yaw   = y;

    fcount++;
  }

  return(foffsets);
}



int formatdata(string dir, int ns, char *namefield, double *data, unsigned int *data_out, string &bolofield)
{

  int ii;

   FILE *fpformat;
   char lineformat[100];
   char tempchar[100];
   char boloname[100];
   string formatname = "0";

   double fact1, scal1;

   formatname = dir+string("format");

   fpformat = fopen(formatname.c_str(),"r");

   while (string(lineformat) != string(namefield) && getc(fpformat) != EOF){
     fscanf(fpformat,"%s",lineformat);
   }
   if (getc(fpformat) == EOF){
     printf("ERROR: field does not exist in the format file");
     return 0;
   }

   fscanf(fpformat,"%s",lineformat);
   fscanf(fpformat,"%s",tempchar);
   fscanf(fpformat,"%s",boloname);
   fscanf(fpformat,"%lf%lf",&fact1,&scal1);

   // conversion
   for (ii=0;ii<ns;ii++)
     data_out[ii] = (unsigned int)((data[ii]-scal1)/fact1+0.5);

   bolofield=boloname;

   return 1;

}



// int read_data(string fname, int frame, int fs, int ns,
//               void* data, string field, char type)
// {
//   int error_code, nread;

//   char ffname[100];
//   strcpy(ffname,fname.c_str());

//   nread = GetData(ffname, field.c_str(),
//                     frame,fs, /* 1st sframe, 1st samp */
//                     0, ns, /* num sframes, num samps */
//                     type, data,
//                     &error_code);

//   if (error_code != GD_E_OK) {
//     cerr << "    GetData Error while reading "<< field
// 	 << " from " << fname <<":\n";
//     cerr << GD_ERROR_CODES[error_code] << "\n";
//     cerr << " Frame: " << frame << "\n";

//     //exit(0);
//   }

//   if(nread == 0) {
//     cerr << "Warning: nread = 0\n";
//     //exit(0);
//   }

//   return nread;
// }




int read_data_std(string fname, int frame, int fs, int ns,
              void* data, string field, char type)
{

  int sizetype;
  char test[2];
  test[0] = type;
  test[1] = '\0';
  string typestr = string(test);
  //  printf("type = %s\n",test);


  FILE *fp;

  if (typestr == "d") sizetype = 8;
  if (typestr == "c") sizetype = 1;

  string filename = fname + field;
  fp = fopen(filename.c_str(),"r");
  fseek(fp,(20*frame+fs)*sizetype,SEEK_SET);
  fread(data,sizetype,ns,fp);
  fclose(fp);

  return 1;
}



int write_data(string& fname, int frame, int fs, int ns, void* data, string& field, char type, int samples_per_frame = 1)

{
  int file;
  int offs;
  size_t size = 0;
  size_t s_type;
  string out_file;
  string format_field;

  if(samples_per_frame != 1 && samples_per_frame != 20) {
    cerr << "write_data: samples per frame not valid: ";
    cerr << samples_per_frame << " (must be either 1 or 20)\n";
    exit(0);
  }

  switch (type) {
  case 'c':
    s_type = 1;
    break;
  case 's': case 'u':
    s_type = 2;
    break;
  case 'S': case 'U': case 'f': case 'i':
    s_type = 4;
    break;
  case 'd':
    s_type = 8;
    break;
  default:
    cerr << "Error:write_data. type format not valid" << endl;
    exit(0);
  }

  out_file = fname + "/" + field;

  file = open(out_file.c_str(),  O_RDWR | O_CREAT, 00644);
  if (file < 0) {
    cerr << "Write() : Error opening file " << out_file << endl;
    return 0;
  }

  offs = (frame * samples_per_frame + fs) * s_type;
  lseek(file, offs, SEEK_SET);

  size = ns * s_type;
  ns = write(file, data, size);

  close(file);

  return  ns;
}





double* randg(long nombre, int seedpass) {

  double* nombre_hasard;
  time_t temps;
  temps = time(NULL);

  unsigned int seed = 0;

  if (seedpass == 0) seed = (unsigned int) temps;
  if (seedpass != 0 && seedpass != -1) seed = (unsigned int) seedpass;
  if (seedpass != -1) srandom(seed);

  nombre_hasard= new double[nombre];

  for (long i=0;i<nombre/2;i++) {
    double t1 = (double(rand())/RAND_MAX);
    double t2 = (double(rand())/RAND_MAX);
    nombre_hasard[2*i]=sqrt(-2*log(t1))*cos(2*M_PI*t2);
    nombre_hasard[2*i+1]=sqrt(-2*log(t1))*sin(2*M_PI*t2);
  }

  if (nombre/2!=nombre/2.) {
    double t1 = (double(rand())/RAND_MAX);
    double t2 = (double(rand())/RAND_MAX);
    nombre_hasard[nombre-1]=sqrt(-2*log(t1))*cos(2*M_PI*t2);
  }


  return nombre_hasard;
}


double* rand(long nombre, int seed) {

  double* nombre_hasard;
  time_t temps;
  temps = time(NULL);
  if (seed == 0) seed = (unsigned int) temps;
  if (seed != -1) srandom(seed);

  nombre_hasard= new double[nombre];

  for (long i=0;i<nombre;i++) nombre_hasard[i]=(double(rand())/RAND_MAX);

  return nombre_hasard;
}


