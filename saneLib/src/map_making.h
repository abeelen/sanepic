#ifndef MAP_MAKING
#define MAP_MAKING


#include <string>
#include <fftw3.h>
#include "struct_definition.h"


void compute_PtNmd(double *data, double *Nk, long ndata, long NAXIS1, long NAXIS2,
		long long *indpix, long long *samptopix, long long npix, double *PNd);


void compute_diagPtNP(double *Nk, long long *samptopix, long ndata,
		long NAXIS1, long NAXIS2, long long *indpix,
		long npix, double f_lppix, double *dPtNP);


void compute_diagPtNPCorr(double *Nk, long long *samptopix, long ndata,
		long NAXIS1, long NAXIS2, long long *indpix,
		long long npix, double f_lppix, double *dPtNP);

void MapMakePreProcessData(double *data,  int *flag, long ns, struct param_sanePre proc_param, double f_lppix, double *data_lp, double *Ps);


void noisepectrum_estim(double *data, long ns, double *ell, int nbins, double fsamp,
		double *bfilter, double *Nell, double *Nk);


void noisecrosspectrum_estim(fftw_complex *fdata1, fftw_complex *fdata2, int ns, double *ell, int nbins, double fsamp, double *bfilter, double *Nell, double *Nk);


int readNSpectrum(std::string nameSpfile, double *bfilter, long ns, double fsamp, double *Nk);


void deproject(double *S, long long *indpix, long long *samptopix, long long ndata, long NAXIS1, long NAXIS2,
		long long npix, double *Ps, int flgdupl, int factdupl, long ntotscan = 0,
		long long *indpsrc = NULL, long long npixsrc = 0);


int compare_global_array_long_long (const void *array_1, const void *array_2);
int compare_long_long (const void *a, const void *b);


#endif
