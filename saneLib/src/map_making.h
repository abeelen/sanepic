#ifndef MAP_MAKING
#define MAP_MAKING


#include <string>
#include <fftw3.h>
#include "struct_definition.h"


//! long long *data_compare



void compute_PtNmd(double *data, double *Nk, long ndata, long NAXIS1, long NAXIS2,
		long long *indpix, long long *samptopix, long long npix, double *PNd);


void compute_diagPtNP(double *Nk, long long *samptopix, long ndata,
		long NAXIS1, long NAXIS2, long long *indpix,
		long npix, double fhp_pix, double *dPtNP);


void compute_diagPtNPCorr(double *Nk, long long *samptopix, long ndata,
		long NAXIS1, long NAXIS2, long long *indpix,
		long long npix, double fhp_pix, double *dPtNP);


//! Subtract deprojected signal to the data and apply a pre-process to the data
/*!
 * Pre-process is : fill gaps, remove a polynomia, remove a baseline, apply butterworth filter and re-fill the gaps. Re-add deprojected signal to finish
 \param data Input data to be pre-processed (changed in place)
 \param flag An array specifying which "data" samples are flagged
 \param ns Number of samples for the considered scan : samples_struct.fitsvect[iframe]
 \param proc_param The param_saneProc structure
 \param fhp_pix High-pass Filter cut-off frequency (converted in samples)
 \param Ps Deprojected signal to be removed from data
 */
void MapMakePreProcessData(double *data,  int *flag, long ns, struct param_saneProc proc_param, double * bfilter, double *Ps);


void noisepectrum_estim(double *data, long ns, double *ell, int nbins, double fsamp,
		double *bfilter, double *Nell, double *Nk);


void noisecrosspectrum_estim(fftw_complex *fdata1, fftw_complex *fdata2, int ns, double *ell, int nbins, double fsamp, double *bfilter, double *Nell, double *Nk);


int readNSpectrum(std::string nameSpfile, double *bfilter, long ns, double fsamp, double *Nk);

//! Deproject map signal to a timeline format
/*!
 \param S Map signal
 \param indpix The pixels indices table
 \param samptopix Sample to pixel projection matrix
 \param ndata Number of samples for the considered scan : samples_struct.fitsvect[iframe]
 \param NAXIS1 Number of horizontal pixels (determined by sanePos)
 \param NAXIS2 Number of vertical pixels (determined by sanePos)
 \param npix Number of filled pixels (in the map)
 \param Ps Deprojected signal
 \param flgdupl True if flagged data are put in a separated map
 \param factdupl worth 2 if flagged data are put in a separated map, 1 otherwise
 \param ntotscan Total number of input scans
 \param indpsrc The masked pixels indices table
 \param npixsrc The number of pixels in the mask
 */
void deproject(double *S, long long *indpix, long long *samptopix, long long ndata, long NAXIS1, long NAXIS2,
		long long npix, double *Ps, int flgdupl, int factdupl, long ntotscan = 0,
		long long *indpsrc = NULL, long long npixsrc = 0);

//! Cast to long and compare data_compare[*a] to data_compare[*b]
/*!
 * This routine is used with standard sort function to get an array sorted values indices
 \param a A pointer to a dat_compare index to be compared
 \param b A pointer to a dat_compare index to be compared
 \return An integer which is >0 if data_compare[*array_1] > data_compare[*array_2], <0 otherwise
 */
int compare_global_array_long_long (const void *array_1, const void *array_2);

//! Cast to long and compare a variable "a" to another variable "b"
/*!
 * This routine is used with standard sort function
 \param a A pointer to a value to be compared with "b"
 \param b A pointer to a value to be compared with "a"
 \return An integer which is >0 if *a > *b, <0 otherwise
 */
int compare_long_long (const void *a, const void *b);


#endif
