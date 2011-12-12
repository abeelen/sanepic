#ifndef TODPROCESS
#define TODPROCESS


#include <string>
#include <vector>



extern "C" {
#include "nrutil.h"
}

//! todprocess.cpp macros for faster simple comparison between 2 variables
/*! SWAP 2 variables "a" and "b" : a will be equal to old b value and b to old a value */
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

//! todprocess.cpp macros for faster simple comparison between 2 variables
/*! MIN is equal to the minimum value between the 2 input variables "a" and "b" */
#define MIN(a,b) (((a)<(b))?(a):(b))

//! todprocess.cpp macros for faster simple comparison between 2 variables
/*! MAX is equal to the maximum value between the 2 input variables "a" and "b" */
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
 * \param f_hp frequency of this high pass filter applied to the data
 * \param orderB Butterworth filter order
 * \param yout The output data : butterworth filter applied to y
 * \param bfilter Butterworth filter values in Fourier space
 * \param apodize A boolean specifying if an apodization has to be done before applying the filter
 * \param napod Number of samples to apodize
 * \param overwrite A boolean specifying whether y is overwritten (yout is not used), or not (output array is yout)
 */
void butterworth_filter(int ndata, double f_hp, int orderB, double *bfilter);
void butterworth(double y[], int ndata, double *yout, double *bfilter, bool apodize, int napod, bool overwrite);

//! Compute an apodization window
/*!
 * \param ns Window size (also data size)
 * \param nn For [0,nn] and [nn,ns], data are apodized, otherwise, window value is set to 1 and data are not modified
 * \return Apodization window values as an array of double, which size is ns
 */
double* apodwindow(int ns, int nn);

//! Interpolate the input binary spectrum to obtain a continous spectrum over ns samples
/*!
 * Interpolate logarithmically the noise power spectrum
 * This routine can compute noise power spectrum for a given mode (not used in Sanepic actually), or for all modes
 \param ell The spectra bins
 \param SpN binned Spectrum
 \param bfilter Butterworth filter values in Fourier space
 \param nbins Number of spectra bins
 \param ns Considered scan's number of samples
 \param fsamp Instrument Sampling frequency
 \param Nk Interpolated Spectrum binned, computed by this routines
 \param mode If mode != NULL, Nk is computed for the given mode value
 */
void binnedSpectrum2log_interpol(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode = NULL);

//! Interpolate the input binary spectrum to obtain a continous spectrum over ns samples
/*!
 * Linear interpolation, because, by convention, some spectrum can be negatives !!!
 * This routine can compute noise power spectrum for a given mode (not used in Sanepic actually), or for all modes
 \param ell The spectra bins
 \param SpN binned Spectrum
 \param bfilter Butterworth filter values in Fourier space
 \param nbins Number of spectra bins
 \param ns Considered scan's number of samples
 \param fsamp Instrument Sampling frequency
 \param Nk Interpolated Spectrum binned, computed by this routines
 \param mode If mode != NULL, Nk is computed for the given mode value
 */
void InvbinnedSpectrum2log_interpol(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode = NULL);

#endif
