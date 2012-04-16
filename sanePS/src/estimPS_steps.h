#ifndef ESTIMPS_STEPS_H_
#define ESTIMPS_STEPS_H_

#include <string>
#include <vector>
#include "mpi_architecture_builder.h"

#include <fftw3.h>


int common_mode_computation(struct samples samples_struct,  struct param_saneProc proc_param, struct param_sanePos pos_param,
		struct param_common dir, long iframe, double *apodwind, long NAXIS1, long NAXIS2, long long npix,
		double *S, long long *indpix,double **mixmat, long ncomp, double **commonm2,
		double &factapod);

int estimate_noise_PS(struct samples samples_struct, struct param_saneProc proc_param, struct param_sanePos pos_param,
		struct param_common dir, long iframe, long &nbins,	long &nbins2, long NAXIS1,
		long NAXIS2, long long npix, double *&ell, double *S, long long *indpix,
		double *apodwind, long ncomp, double **mixmat, double **commonm2,
		double factapod,double **Rellth, double **N, double **P);


int estimate_CovMat_of_Rexp(struct samples samples_struct, struct param_common dir, long iframe, long nbins, double *ell, long ncomp, double **mixmat,
		double factapod,double **Rellexp, double **N, double **P, double *SPref, int rank);

int expectation_maximization_algorithm(double fcut, long nbins, long ndet, long ncomp,long ns, double fsamp,
		std::string outdirSpN,	double **Rellexp, double **Rellth, double **mixmat,double **P,double **N,
		double *SPref, double *ell, int rank);

double fdsf(double **Rellexp, double *w, double **A, double **P, double **N, long ndet, long ncomp, long nbins);

void rescaleAP(double **A, double **P, long ndet, long ncomp, long nbins);

int write_to_disk(std::string outdirSpN, std::string fits_filename, struct param_sanePS structPS, std::vector<std::string> det, long nbins, double *ell, double **mixmat,
		double **Rellth, double **Rellexp, double **N, double *SPref, double **P);

#endif /* ESTIMPS_STEPS_H_ */
