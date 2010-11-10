#ifndef TODPROCESS
#define TODPROCESS


#include <string>
#include <vector>



extern "C" {
#include "nrutil.h"
}

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void init2D_double(double **A, long im, long jm, long nx, long ny, double val);

void remove_poly(double y[], long ndata, int norder, double* yout, int *flag = NULL);
void fillgaps2(double data[], long ns, double* yout,  int* flag, int taille);

void butterworth(double y[], int ndata, double f_lp, int orderB, double *yout, double *bfilter, bool apodize, int napod, bool overwrite);
double* apodwindow(int ns, int nn);

void binnedSpectrum2log_interpol(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode = NULL); // Patanchon
void InvbinnedSpectrum2log_interpol(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode = NULL); // Patanchon
void InvbinnedSpectrum2bis(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode = NULL);

//double* randg(long nombre, int seedpass);
//double* rand(long nombre, int seed);

#endif
