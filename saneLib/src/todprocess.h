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

void cutdata(double y[], int indm, int indp, double *yout);
void cutdata(unsigned char y[], int indm, int indp, unsigned char *yout);
void mergedata(double y1[], int ndata1, double y2[], int ndata2, double *yout);
void dindgen(int nn, double *y);

void fillgaps2(double data[], long ns, double* yout,  int* flag, int taille);

void butterworth(double y[], int ndata, double f_lp, int orderB, double *yout, double *bfilter, bool apodize, int napod, bool overwrite);
double* apodwindow(int ns, int nn);


void binnedSpectrum2log_interpol(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode = NULL); // Patanchon
void InvbinnedSpectrum2log_interpol(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode = NULL); // Patanchon
void InvbinnedSpectrum2bis(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode = NULL);


int read_data_std(std::string fname, int frame, int fs, int ns, void* data, std::string field, char type);
int read_data(std::string fname, int frame, int fs, int ns, void* data, std::string field, char type);


double randg_value(long nombre, int seedpass);
double* randg(long nombre, int seedpass);
double* rand(long nombre, int seed);

#endif
