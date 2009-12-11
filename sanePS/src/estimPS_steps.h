/*
 * estimPs_steps.h
 *
 *  Created on: 19 août 2009
 *      Author: matthieu
 */

#ifndef ESTIMPS_STEPS_H_
#define ESTIMPS_STEPS_H_

//#include "covMatrixIO.h"
//#include "psdIO.h"
//#include <vector>
//#include <time.h>
//#include <algorithm>
//#include "todprocess.h"
//#include "map_making.h"
//#include "inline_IO2.h"
//#include <sstream>
#include <string>
#include <vector>
#include "mpi_architecture_builder.h"

#include <fftw3.h>


void read_mixmat_file(std::string MixMatfile, std::string dir, double **&mixmat, long ndet, long ncomp);

void common_mode_computation(struct detectors det, struct param_process proc_param, struct param_positions pos_param,
		struct directories dir, double *apodwind,long ns, long ff, long NAXIS1, long NAXIS2, long long npix,
		long iframe, double *S, long long *indpix,double **mixmat, long ncomp, double **commonm2,
		double &factapod, std::string fits_filename);

//void estimate_noise_PS(std::vector<std::string> bolonames, std::string dirfile, std::string extentnoiseSp, std::string noiseSppreffile,/* string bextension, string fextension,*/ long &nbins,
//		long &nbins2, long ns, long ff, long ndet, long NAXIS1,long NAXIS2, long long npix,long napod, double *&ell, double *data, /*short *flag,*/
//		long long *samptopix, std::string dir, double *S, long iframe, double *Ps, double *data_lp, double*bfilter, long long *indpix, bool NORMLIN,
//		bool NOFILLGAP, bool remove_polynomia,bool flgdupl, int factdupl, double *apodwind, long ncomp, double **mixmat, double **commonm2, double fsamp,
//		double *Nk, double *Nell, double factapod,double **Rellth, double **N, double *commontmp, double **P,  std::string outdirSpN,std::string fits_filename);

void estimate_noise_PS(struct detectors det, struct param_process proc_param, struct param_positions pos_param,
		struct directories dir, long &nbins,	long &nbins2, long ns, long ff, long NAXIS1,
		long NAXIS2, long long npix, double *&ell, double *S, long iframe,long long *indpix,
		double *apodwind, long ncomp, double **mixmat, double **commonm2,
		double factapod,double **Rellth, double **N, double **P, std::string fits_filename);


//void estimate_CovMat_of_Rexp(long nbins, long ns, long ff, long ndet, double *ell, std::string dir, long ncomp, double **mixmat,double fsamp,
//		double *Nk, double *Nell, double factapod,double **Rellexp, double **N, double **P, std::string outdirSpN, fftw_complex *fdata1, fftw_complex  *fdata2,
//		double *SPref, std::vector<std::string> bolonames);
void estimate_CovMat_of_Rexp(struct directories dir, struct detectors det, long nbins, long ns, long ff, double *ell, long ncomp, double **mixmat,double fsamp,
		double factapod,double **Rellexp, double **N, double **P, double *SPref);

//void expectation_maximization_algorithm(double fcut, long nbins, long ndet, long ncomp,long ns, double fsamp, long ff,
//		std::string outdirSpN, 	double **Rellexp, double **Rellth, double **mixmat,double **P,double **N, double **Cov, double *p,
//		double *uvec, double *ivec, double **iCov, double *SPref, double *ell);

void expectation_maximization_algorithm(double fcut, long nbins, long ndet, long ncomp,long ns, double fsamp,
		long ff, std::string outdirSpN,	double **Rellexp, double **Rellth, double **mixmat,double **P,double **N,
		double *SPref, double *ell);

double fdsf(double **Rellexp, double *w, double **A, double **P, double **N, long ndet, long ncomp, long nbins);

void rescaleAP(double **A, double **P, long ndet, long ncomp, long nbins);

//void write_to_disk(std::string outdirSpN, long ff, std::vector<std::string> bolonames,
//		long nbins, double *ell, double **mixmat, double **Rellth, double **Rellexp, long ncomp, long ndet,double **N, double *SPref,
//		double **P);

void write_to_disk(std::string outdirSpN, long ff, struct detectors det, long nbins, double *ell, double **mixmat,
		double **Rellth, double **Rellexp, long ncomp,double **N, double *SPref, double **P);

#endif /* ESTIMPS_STEPS_H_ */
