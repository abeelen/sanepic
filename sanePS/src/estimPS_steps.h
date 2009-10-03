/*
 * estimPs_steps.h
 *
 *  Created on: 19 août 2009
 *      Author: matthieu
 */

#ifndef ESTIMPS_STEPS_H_
#define ESTIMPS_STEPS_H_

#include "covMatrixIO.h"
#include "psdIO.h"
#include <vector>
#include <time.h>
#include <algorithm>
#include "todprocess.h"
#include "map_making.h"
#include "inline_IO2.h"
#include <sstream>
#include <fftw3.h>

extern "C" {
#include <fitsio.h>
}


void read_mixmat_file(string MixMatfile, string dir, double **&mixmat, long ndet, long ncomp);

void common_mode_computation(double *apodwind, long ndet, long ns, long ff, int NAXIS1,int NAXIS2, long npix, bool flgdupl, int factdupl, std::vector<string> bolonames,/* string bextension, string fextension,*/
		string dirfile, int shift_data_to_point,  string dir, long iframe, double *S, long *indpix,  bool NORMLIN,
		bool NOFILLGAP, bool remove_polynomia, long napod, double **mixmat, long ncomp, double **commonm2, long *samptopix, double *Ps, double *data, double *data_lp,  /*short *flag,*/
		double *bfilter, double **Cov, double *uvec,double *p,double *ivec, double **iCov, double &factapod ,fftw_complex *fdata1, string fits_filename);


void estimate_noise_PS(std::vector<string> bolonames, string dirfile, string extentnoiseSp, string noiseSppreffile,/* string bextension, string fextension,*/ long &nbins,
		long &nbins2, long ns, long ff, long ndet, int NAXIS1,int NAXIS2, long npix,long napod, double *&ell, double **&SpN_all, double *data, /*short *flag,*/
		long *samptopix, string dir, double *S, long iframe, double *Ps, double *data_lp, double*bfilter, long *indpix, bool NORMLIN,
		bool NOFILLGAP, bool remove_polynomia,bool flgdupl, int factdupl, double *apodwind, long ncomp, double **mixmat, double **commonm2, double fsamp,
		double *Nk, double *Nell, double factapod,double **Rellth, double **N, double *commontmp, double **P, int shift_data_to_point, string outdirSpN,string fits_filename);


void estimate_CovMat_of_Rexp(long nbins, long ns, long ff, long ndet, double *ell, string dir, long ncomp, double **mixmat,double fsamp,
		double *Nk, double *Nell, double factapod,double **Rellexp, double **N, double **P, string outdirSpN, fftw_complex *fdata1, fftw_complex  *fdata2,
		double *SPref);

void expectation_maximization_algorithm(double fcut, long nbins, long ndet, long ncomp,long ns, double fsamp, long ff,
		string outdirSpN, 	double **Rellexp, double **Rellth, double **mixmat,double **P,double **N, double **Cov, double *p,
		double *uvec, double *ivec, double **iCov, double *SPref, double *ell);


double fdsf(double **Rellexp, double *w, double **A, double **P, double **N, long ndet, long ncomp, long nbins);

void rescaleAP(double **A, double **P, long ndet, long ncomp, long nbins);

void write_to_disk(string outdirSpN, long ff, std::vector<string> bolonames,
		long nbins, double *ell, double **mixmat, double **Rellth, double **Rellexp, long ncomp, long ndet,double **N, double *SPref,
		double **P);

#endif /* ESTIMPS_STEPS_H_ */
