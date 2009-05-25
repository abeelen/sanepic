#ifndef TODPROCESS
#define TODPROCESS

#include <iostream>
#include <cmath>
#include <string>
//#include <matpack/matpack.h>
// #include <getdata.h>
//#include "/fir_data/patanch/numrec/inc/nrutil.h"

#include <fcntl.h>



#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


#define BLAST_MAIN_CLK (32.0E6) // Hz
#define ADC_CLK        (BLAST_MAIN_CLK/3072.0)
#define HUNDRED_HZ_CLK (BLAST_MAIN_CLK/319488.0 )

#define DT (1.0/HUNDRED_HZ_CLK)

#define FRAME_FREQUENCY (HUNDRED_HZ_CLK / 20.0)
#define FAST_PER_SLOW 20


#ifndef SQR
#define SQR(a) ((a)*(a))
#endif




using namespace std;


struct foffset {
  long frame;
  float pitch;
  float yaw;
};



void init1D_double(double *A, long im, long n, double val);
void init1D_long(long *A, long im, long n, long val);
void init2D_double(double **A, long im, long jm, long nx, long ny, double val);
void init2D_long(long **A, long im, long jm, long nx, long ny, long val);


double ** dma(int nrl, int nrh, int ncl, int nch);
long ** lma(int nrl, int nrh, int ncl, int nch);
char ** charma(int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_lmatrix(long **m, long nrl, long nrh, long ncl, long nch);
int *ivector(long nl, long nh);
long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_lvector(long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);

void iindexx(long n, long arr[], long indx[]);

void dlfit(double x[], double y[], double sig[], int ndat, double a[], int ia[],
	   int ma, double **covar, double *chisq, void (*funcs)(double, double [], int));
void dlinfit(double y[], double sig[], int ndat, double a[], int ia[],
	     int ma, double **covar, double *chisq, double **arrays);
void dmrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[],
	int ma, double **covar, double **alpha, double *chisq,
	    void (*funcs)(double, double [], double *, double [], int), double *alamda);
void dmrqcof(double x[], double y[], double sig[], int ndata, double a[], int ia[],
	int ma, double **alpha, double beta[], double *chisq,
	    void (*funcs)(double, double [], double *, double [], int));
void dcovsrt(double **covar, int ma, int ia[], int mfit);
void dgaussj(double **a, int n, double **b, int m);
void polynomia(double x, double y[], int dma);
void dpolyfit(double x[], double y[], int ndata, int norder, double *a);
void remove_poly(double y[], int ndata, int norder, double* yout, unsigned char *flag = NULL);
void gausspoly(double x, double a[], double *y, double dyda[], int na);
double dgaussfit(double y[], double *a, double rms, int ndata, int ma, unsigned char *flag, int vfl);

void dlubksb(double **a, int n, int *indx, double b[]);
int dludcmp(double **a, int n, int *indx, double *d);
void dcholdc(double **a, long n, double p[]);
void dcholsl(double **a, long n, double p[], double b[], double x[]);

void cutdata(double y[], int indm, int indp, double *yout);
void cutdata(unsigned char y[], int indm, int indp, unsigned char *yout);
void mergedata(double y1[], int ndata1, double y2[], int ndata2, double *yout);
void dindgen(int nn, double *y);
void sort(double y[], int nn, double *yout, int *nrel, unsigned char* flag = NULL);

void findspike(double y[], int ndata, double transfer[], const int ntr, double* yout, unsigned char* flag, unsigned char* pflag, double** pplanet, int *xerror, int *pos, int *xcosm, int *xspike, double *chisqs, double *amplspike, double *allhsig, int *posp, double *amplp, double *allhsigp);
void fillgaps(double y[], int ndata, double* yout, unsigned char* flag, double sign);
//void deconv_antialias(double y[], int ndata, double f_lp, double* yout, bool apodize);
double* CR_tf(int ndata);

void butterworth(double y[], int ndata, double f_lp, int orderB, double *yout, double *bfilter, bool apodize, int napod, bool overwrite);
double* apodwindow(int ns, int nn);
void memcof(double data[], int n, int m, double d[]);
void Pad(double* d, int window, int samples, int run);

void binnedSpectrum2log_interpol(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode = NULL);
void InvbinnedSpectrum2log_interpol(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode = NULL);

int compare( const void *arg1, const void *arg2 );
int compare_doubles (const void *a, const void *b);
int compare_long (const void *a, const void *b);

void filter_despike(double y[], int ndata, double* yout, unsigned char *flag);

foffset* read_mapoffsets(string fname, float *scoffsets, int *nfoff);
int formatdata(string dir, int ns, char *namefield, double *data, unsigned int *data_out, string &bolofield);
int read_data_std(string fname, int frame, int fs, int ns, void* data, string field, char type);
int read_data(string fname, int frame, int fs, int ns, void* data, string field, char type);
int write_data(string& fname, int frame, int fs, int ns, void* data, string& field, char type, int samples_per_frame);

double* randg(long nombre, int seedpass);
double* rand(long nombre, int seed);

void minmax(double* data, int ndata, double *min, double *max, int *posmin, int *posmax, unsigned char *flag = NULL);




#endif