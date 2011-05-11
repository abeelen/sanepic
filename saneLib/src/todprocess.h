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

//! Partially or totally initialize an array of double values with "val"
/*!
 * To fill "A" totally, im/jm are set to 0 and nx/ny are set to "A" dimensions
 * \param A is a double matrix to be initializes
 * \param im is the first line to be filled indice
 * \param jm is the first column to be filled indice
 * \param nx is the number of lines to be filled
 * \param ny is the number of columns to be filled
 * \param val is the value that is used to fill "A"
 */
void init2D_double(double **A, long im, long jm, long nx, long ny, double val);

//! Fit and Remove a polynomial value from an array of double
/*!
 * \param y is the input array
 * \param ndata is y size
 * \param norder is the polynomial order
 * \param yout is the output array (after polynomial subtraction) : y - poly
 * \param flag is an optional field : y flagged samples are not taken into account for polynomial fit
 */
void remove_poly(double y[], long ndata, int norder, double* yout, int *flag = NULL);

//! Fill timeline gaps using a locally fitted polynomia
/*!
 * polynomial order used to fill data is 1
 * \param data The timeline which gaps will be filled
 * \param ns the number of samples in data
 * \param yout is the ouput data : filled data array is stored in yout
 * \param flag is an optional field : flagged samples are to be filled
 * \param taille is the maximum number of "good" samples that will be used to locally fit the polynomia
 */
void fillgaps2(double data[], long ns, double* yout,  int* flag, int taille);

//! Compute a butterworth filter and apply it to the data "y"
/*!
 * \param y is the input data array to be filtered
 * \param ndata the number of samples in y
 * \param f_lp frequency of this high pass filter applied to the data
 * \param orderB
 * \param taille is the maximum number of "good" samples that will be used to locally fit the polynomia
 */
void butterworth(double y[], int ndata, double f_lp, int orderB, double *yout, double *bfilter, bool apodize, int napod, bool overwrite);


double* apodwindow(int ns, int nn);


void binnedSpectrum2log_interpol(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode = NULL); // Patanchon


void InvbinnedSpectrum2log_interpol(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode = NULL); // Patanchon

#endif
