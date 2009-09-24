#ifndef TODPROCESS
#define TODPROCESS

#include <iostream>
//#include <cmath>
#include <cstdlib>
#include <string>
//#include <fcntl.h>
#include <fstream>
#include <cmath>
#include <fftw3.h>
//#include <time.h>
//#include <unistd.h>
#include <vector>
#include <gsl/gsl_math.h>

extern "C" {
#include "nrutil.h"
#include "nrcode.h"
}



using namespace std;



//void init1D_double(double *A, long im, long n, double val); // on peut garder
//void init1D_long(long *A, long im, long n, long val); // on peut garder
void init2D_double(double **A, long im, long jm, long nx, long ny, double val); // on peut garder


template<class T> void vector2array(std::vector<T> vect, T* a)
{
	// copy list of type T to array of type T
	typename std::vector<T>::iterator iter;
	int i;

	for (iter=vect.begin(), i=0; iter != vect.end(); iter++, i++) {
		a[i] = *iter;
	}
}


/*
double ** dmatrix(int nrl, int nrh, int ncl, int nch); // NR
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch); // NR
int *ivector(long nl, long nh); // NR
double *dvector(long nl, long nh); // NR
void free_ivector(int *v, long nl, long nh); // NR
void free_dvector(double *v, long nl, long nh); // NR
 */



void polynomia(double x, double y[], int dma); // on peut garder
void remove_poly(double y[], int ndata, int norder, double* yout, short *flag = NULL); // viens de NR ?? TODO : verifier que ca marche avec flag short
void dpolyfit(double x[], double y[], int ndata, int norder, double *a); // Nr a verifier


void cutdata(double y[], int indm, int indp, double *yout); // on peut garder
void cutdata(unsigned char y[], int indm, int indp, unsigned char *yout); // on peut garder
void mergedata(double y1[], int ndata1, double y2[], int ndata2, double *yout); // on peut garder
void dindgen(int nn, double *y); // on peut garder


void fillgaps(double y[], int ndata, double* yout,  short* flag, double sign); // fait par Patanchon


void butterworth(double y[], int ndata, double f_lp, int orderB, double *yout, double *bfilter, bool apodize, int napod, bool overwrite); // on garde ?
double* apodwindow(int ns, int nn); // on peut garder


void binnedSpectrum2log_interpol(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode = NULL); // Patanchon
void InvbinnedSpectrum2log_interpol(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode = NULL); // Patanchon
void InvbinnedSpectrum2bis(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode = NULL);

int compare_long (const void *a, const void *b); // on peut garder


int read_data_std(string fname, int frame, int fs, int ns, void* data, string field, char type); // on garde
int read_data(string fname, int frame, int fs, int ns, void* data, string field, char type); // on garde


double randg_value(long nombre, int seedpass); // on garde
double* rand(long nombre, int seed);	// on garde

//void minmax(double* data, int ndata, double *min, double *max, int *posmin, int *posmax, unsigned char *flag = NULL); // on garde




#endif
