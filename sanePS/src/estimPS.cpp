#include "estimPS.h"

#include <iostream>

#include <vector>
#include <string>

#include "todprocess.h"
#include "inputFileIO.h"
#include "estimPS_steps.h"


extern "C" {
#include <fitsio.h>
#include "nrutil.h"
}


using namespace std;



int EstimPowerSpectra(struct param_process proc_param,struct detectors det,struct common dir, struct param_positions pos_param,
		long ns, long NAXIS1, long NAXIS2, long long npix, long iframe,
		long long *indpix, double *S, string MixMatfile,string ellFile, string fits_filename, long ncomp, double fcut, int rank)
{

	long nbins = 500; // temp value
	long nbins2; // readed bins in mixmatfile
	double factapod= 0.0; // apodization factor

	double *apodwind, *SPref;
	double **commonm2, **vect; // ?
	double **P, **N, **Rellth, **Rellexp; // Rellth = theorical covariance matrix, Rellexp = experimental cov mat
	//	Nk = Noise power spectrum
	double *ell;// bins values, Nell = binned noise PS
	double **mixmat; // mixing matrix



	printf("\nEstimation procedure started : \n");

	//	data = raw data
	//	data_lp = data low passed
	//	samptopix = sample to pixel projection matrix
	//	Nk = noise PS
	//	bfilter = butterworth filter values
	//	fdata = fourier transform

	if(read_double(ellFile, ell, nbins)) // read ell in ellfile
		return 1;
	nbins = nbins-1;

	//	Nell = binned noise PS
	SPref = new double[nbins]; // first detector PS to avoid numerical problems
	P = dmatrix(0,ncomp-1,0,nbins-1); // component power spectra
	N = dmatrix(0,det.ndet-1,0,nbins-1); // uncorralated part of the noise
	Rellexp = dmatrix(0,(det.ndet)*(det.ndet)-1,0,nbins-1); // covariance estimated from signal (measured)
	Rellth = dmatrix(0,(det.ndet)*(det.ndet)-1,0,nbins-1); // covariance matrix (Theory) = pattern : invertible but not Rellexp

	// sign = sigma of the noise
	// Cov = AtN-1A
	// iCov = inverted AtN-1A

	vect = dmatrix(0,ncomp-1,0,det.ndet-1);
	commonm2 = dmatrix(0,ncomp,0,ns-1);

	// One has to initialize the two matrices for each iteration...
	init2D_double(Rellexp,0,0, (det.ndet)*(det.ndet),nbins ,0.0);
	init2D_double(Rellth,0,0, (det.ndet)*(det.ndet),nbins ,0.0);

	init2D_double(commonm2,0,0,ncomp,ns,0.0);

	//apodization
	apodwind = apodwindow(ns,int(ns*0.04));


	for (long ii=0;ii<ns;ii++)
		factapod += apodwind[ii]*apodwind[ii]/ns; // apodization factor

	//----------------------------------- READ MIXMAT PART -------------------------------//
	cout << "[ " << rank << " ] 1/6 - Reading Mixing Matrix" << endl;
	if(read_mixmat_file(MixMatfile, mixmat, det.ndet,ncomp))
		return 1;

	//----------------------------------- COMMON MODE -------------------------------//
	cout << "[ " << rank << " ] 2/6 - Common Mode Computation" << endl;
	if(common_mode_computation(det,proc_param, pos_param, dir, apodwind, ns, NAXIS1, NAXIS2, npix, iframe, S, indpix,
			mixmat, ncomp, commonm2, factapod, fits_filename)) // return commonm2
		return 1;

	//----------------------------------- ESTIMATE NOISE PS -------------------------------//
	cout << "[ " << rank << " ] 3/6 - Estimation of Noise Power Spectrum" << endl;
	if(estimate_noise_PS(det,  proc_param, pos_param, dir, nbins, nbins2, ns, NAXIS1,
			NAXIS2, npix, ell, S, iframe,indpix, apodwind, ncomp, mixmat, commonm2,
			factapod,Rellth, N, P, fits_filename))
		return 1;

	//----------------------------------- ESTIMATE COVMAT of the DATA R_exp -------------------------------//
	cout << "[ " << rank << " ] 4/6 - Estimation of Covariance Matrix" << endl;
	if(estimate_CovMat_of_Rexp(iframe, dir, det, nbins, ns, ell, ncomp, mixmat, proc_param.fsamp,
			factapod, Rellexp, N, P, SPref, fits_filename, rank))
		return 1;

	//----------------------------------- FIT COMPONENT, PS and MIXMAT -------------------------------//
	cout << "[ " << rank << " ] 5/6 - Expectation Maximization" << endl;
	if(expectation_maximization_algorithm(fcut, nbins, det.ndet, ncomp, ns, proc_param.fsamp,
			dir.output_dir,  Rellexp, Rellth, mixmat, P, N, SPref, ell, rank))
		return 1;
	//----------------------------------- WRITE TO DISK -------------------------------//
	cout << "[ " << rank << " ] 6/6 - Saving to disk" << endl;
	if(write_to_disk(dir.output_dir, fits_filename, det, nbins, ell, mixmat, Rellth,
			Rellexp, ncomp, N, SPref,P))
		return 1;
	//----------------------------------- END OF ESTIMPS -------------------------------//

	// clean up
	delete [] SPref;
	delete [] apodwind;
	free_dmatrix(Rellexp,0,(det.ndet)*(det.ndet)-1,0,nbins-1);
	free_dmatrix(Rellth,0,(det.ndet)*(det.ndet)-1,0,nbins-1);
	free_dmatrix(mixmat,0,det.ndet-1,0,ncomp-1);

	free_dmatrix(commonm2,0,ncomp,0,ns-1);
	free_dmatrix(vect,0,ncomp-1,0,det.ndet-1);
	delete [] ell;
	free_dmatrix(P,0,ncomp-1,0,nbins-1);
	free_dmatrix(N,0,det.ndet-1,0,nbins-1);

	return 0;

}






