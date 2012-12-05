#include <iostream>

#include <vector>
#include <string>

#include "todprocess.h"
#include "inputFileIO.h"
#include "estimPS_steps.h"
#include "estimPS.h"
#include "temporary_IO.h"
#include "crc.h"
#include "error_code.h"
#include "utilities.h"

#include <gsl/gsl_math.h>

extern "C" {
#include <fitsio.h>
#include "nrutil.h"
}

using namespace std;

int EstimPowerSpectra(struct param_saneProc proc_param, struct param_common dir, struct param_sanePos pos_param, struct param_sanePS structPS, struct samples samples_struct,
		long iframe, long NAXIS1, long NAXIS2, long long npix, long long *indpix, double *S, int rank)
{

	long nbins; // temp value
	long nbins2; // readed bins in mixmatfile
	double factapod= 0.0; // apodization factor

	int completed_step=0;
	int gotostep=0;

	double *apodwind, *SPref;
	double **commonm2; // ?
	double **P, **N, **Rellth, **Rellexp;	// Rellth = theorical covariance matrix, Rellexp = experimental cov mat
	//	Nk = Noise power spectrum , Nell = binned noise PS
	double *ell, *km;						// bins values
	double **mixmat; 						// mixing matrix

	std::vector<string> det = samples_struct.bolo_list[iframe];
	long ndet = (long)det.size();

	//	data = raw data
	//	data_lp = data low passed
	//	samptopix = sample to pixel projection matrix
	//	bfilter = butterworth filter values
	//	fdata = fourier transform

	// Read out the ells ...
	std::vector<double> dummy;
	std::string output;
	if ( read_file(output, samples_struct.ell_names[iframe], dummy) )
		return FILE_PROBLEM;
	vDouble2carray(dummy, &ell, &nbins);
	dummy.clear();
	nbins   = nbins-1;

	// .. transform them into km once for all
	km = new double[nbins];
	for (int ii=0;ii<nbins;ii++)
		km[ii] = exp((log(ell[ii+1])+log(ell[ii]))/2.0)*samples_struct.nsamples[iframe]/samples_struct.fsamp[iframe];

	SPref   = new double[nbins];						// first detector PS to avoid numerical problems
	P       = dmatrix(0,structPS.ncomp-1,0,nbins-1);	// component power spectra
	N       = dmatrix(0,ndet-1,0,nbins-1);				// uncorralated part of the noise
	Rellexp = dmatrix(0,(ndet)*(ndet)-1,0,nbins-1);		// covariance estimated from signal (measured)
	Rellth  = dmatrix(0,(ndet)*(ndet)-1,0,nbins-1);		// covariance matrix (Theory) = pattern : invertible but not Rellexp

	// sign = sigma of the noise
	// Cov = AtN-1A
	// iCov = inverted AtN-1A

	commonm2 = dmatrix(0,structPS.ncomp,0,samples_struct.nsamples[iframe]-1);

	// One has to initialize the two matrices for each iteration...
	init2D_double(Rellexp,0,0, (ndet)*(ndet),nbins ,0.0);
	init2D_double(Rellth,0,0, (ndet)*(ndet),nbins ,0.0);
	init2D_double(N,0,0,ndet,nbins,0.0);
	init2D_double(P,0,0,structPS.ncomp,nbins,0.0);
	init2D_double(commonm2,0,0,structPS.ncomp,samples_struct.nsamples[iframe],0.0);

	//apodization
	apodwind = apodwindow(samples_struct.nsamples[iframe],proc_param.napod);

	for (long ii=0;ii<samples_struct.nsamples[iframe];ii++)
		factapod += gsl_pow_2(apodwind[ii]);
	factapod /= samples_struct.nsamples[iframe]; // apodization factor


	if(structPS.restore)
		if(restore_session(dir.tmp_dir, FitsBasename(samples_struct.fitsvect[iframe]), gotostep, commonm2, N, P, Rellexp, Rellth, SPref, ndet, structPS.ncomp, nbins, samples_struct.nsamples[iframe])){
			cout << "ERROR restoring data... Exiting ...\n";
			return 1;
		}

	//----------------------------------- READ MIXMAT PART -------------------------------//
#ifdef DEBUG
	cout << "[ " << rank << " ] 1/6 - Reading Mixing Matrix" << endl;
#endif

	if(readMixmatTxt(samples_struct.mix_names[iframe], ndet, structPS.ncomp, mixmat))
		return 1;

	if(gotostep>=2){
#ifdef DEBUG
		cout << "skipping step 2\n";
#endif
	}else{
		//----------------------------------- COMMON MODE -------------------------------//
#ifdef DEBUG
		cout << "[ " << rank << " ] 2/6 - Common Mode Computation" << endl;
#endif
		if(common_mode_computation(samples_struct, proc_param, pos_param, dir, iframe, \
				apodwind,  NAXIS1, NAXIS2, npix, S, indpix, \
				mixmat, structPS.ncomp, commonm2, factapod)) // return commonm2
			return 1;


		completed_step=2;
		if(structPS.save_data){
			save_session(dir.tmp_dir, FitsBasename(samples_struct.fitsvect[iframe]), completed_step, commonm2, N, P, Rellexp, Rellth, SPref, ndet, structPS.ncomp, nbins, samples_struct.nsamples[iframe]);
		}
	}



	if(gotostep>=3){
#ifdef DEBUG
		cout << "skipping step 3\n";
#endif
	}else{
		//----------------------------------- ESTIMATE NOISE PS -------------------------------//
#ifdef DEBUG
		cout << "[ " << rank << " ] 3/6 - Estimation of Noise Power Spectrum" << endl;
#endif
		if(estimate_noise_PS(samples_struct, proc_param, pos_param, dir, iframe, \
				nbins, nbins2,  NAXIS1, NAXIS2, npix, km,
				S, indpix, apodwind, structPS.ncomp, mixmat, commonm2, \
				factapod,Rellth, N, P))
			return 1;

		completed_step=3;
		if(structPS.save_data)
			save_session(dir.tmp_dir, FitsBasename(samples_struct.fitsvect[iframe]), completed_step, commonm2, N, P, Rellexp, Rellth, SPref, ndet, structPS.ncomp, nbins, samples_struct.nsamples[iframe]);
	}

	if(gotostep>=4){
#ifdef DEBUG
		cout << "skipping step 4\n";
#endif
	}else{
		//----------------------------------- ESTIMATE COVMAT of the DATA R_exp -------------------------------//
#ifdef DEBUG
		cout << "[ " << rank << " ] 4/6 - Estimation of Covariance Matrix" << endl;
#endif
		if(estimate_CovMat_of_Rexp(samples_struct, dir, iframe, nbins, km, structPS.ncomp, mixmat,
				factapod, Rellexp, N, P, SPref, rank))
			return 1;

		completed_step=4;
		if(structPS.save_data)
			save_session(dir.tmp_dir, FitsBasename(samples_struct.fitsvect[iframe]), completed_step, commonm2, N, P, Rellexp, Rellth, SPref, ndet, structPS.ncomp, nbins, samples_struct.nsamples[iframe]);
	}

	if(gotostep>=5){
#ifdef DEBUG
		cout << "skipping step 5\n";
#endif
	}else{
		//----------------------------------- FIT COMPONENT, PS and MIXMAT -------------------------------//
		//TODO: dummy fcut for the moment...
#ifdef DEBUG
		cout << "[ " << rank << " ] 5/6 - Expectation Maximization" << endl;
#endif
		if(expectation_maximization_algorithm(samples_struct.fsamp[iframe], nbins, ndet, structPS.ncomp, samples_struct.nsamples[iframe], samples_struct.fsamp[iframe],
				dir.output_dir,  Rellexp, Rellth, mixmat, P, N, SPref, ell, rank))
			return 1;

		completed_step=5;
		if(structPS.save_data)
			save_session(dir.tmp_dir, FitsBasename(samples_struct.fitsvect[iframe]), completed_step, commonm2, N, P, Rellexp, Rellth, SPref, ndet, structPS.ncomp, nbins, samples_struct.nsamples[iframe]);
	}

	if(gotostep>=6){
#ifdef DEBUG
		cout << "skipping step 6\n";
#endif
	}else{
		//----------------------------------- WRITE TO DISK -------------------------------//
#ifdef DEBUG
		cout << "[ " << rank << " ] 6/6 - Saving to disk" << endl;
#endif
		if(write_to_disk(dir.output_dir, samples_struct.fitsvect[iframe], structPS, det, nbins, ell, mixmat, Rellth,
				Rellexp, N, SPref,P))
			return 1;

		completed_step=6;
		if(structPS.save_data)
			save_session(dir.tmp_dir, FitsBasename(samples_struct.fitsvect[iframe]), completed_step, commonm2, N, P, Rellexp, Rellth, SPref, ndet, structPS.ncomp, nbins, samples_struct.nsamples[iframe]);
	}
	//----------------------------------- END OF ESTIMPS -------------------------------//

	// clean up
	delete [] SPref;
	delete [] apodwind;
	delete [] ell;
	delete [] km;

	free_dmatrix(Rellexp,0,(ndet)*(ndet)-1,0,nbins-1);
	free_dmatrix(Rellth,0,(ndet)*(ndet)-1,0,nbins-1);
	free_dmatrix(mixmat,0,ndet-1,0,structPS.ncomp-1);

	free_dmatrix(commonm2,0,structPS.ncomp,0,samples_struct.nsamples[iframe]-1);
	free_dmatrix(P,0,structPS.ncomp-1,0,nbins-1);
	free_dmatrix(N,0,ndet-1,0,nbins-1);

	return 0;

}






