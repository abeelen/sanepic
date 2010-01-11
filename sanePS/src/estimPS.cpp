/*
 * estimPS.cpp
 *
 *  Created on: 28 mai 2009
 *      Author: matthieu
 */

#include "estimPS.h"

// temp
#include <iostream>

#include <vector>
#include <string>

#include "todprocess.h"
#include "inputFileIO.h"
#include "estimPS_steps.h"


extern "C" {
#include <fitsio.h>
#include "nrutil.h"
#include "nrcode.h"
}


using namespace std;



void EstimPowerSpectra(struct param_process proc_param,struct detectors det,struct directories dir, struct param_positions pos_param,
		long ns, long ff, long NAXIS1, long NAXIS2, long long npix, long iframe,
		long long *indpix, double *S, string MixMatfile,string ellFile, string extentnoiseSp,string fits_filename, long ncomp, double fcut)
{


	long nbins = 500;
	long nbins2; // readed bins in mixmatfile
	double factapod= 0.0;
	//double sign0=1;


	// data, data low passed, Prewhitened signal ? ,apodization window, ?, ?, butter filter values, first detector PS (in the cross correlation estimation part)
	double /**data, *data_lp, *Ps, *calp, */*apodwind,/**commonm_f, *//**bfilter, */*SPref;
	//	long long *samptopix; // sample to pixel projection matrix
	//unsigned char *flag; // flag array
	//	short *flag;
	double /***commonm,*/ **commonm2,/* **common_f,*/ **vect; // ?
	double **P, **N, **Rellth, **Rellexp; // Rellth = covariance matrix
	//	double /***Cov, **iCov, **iCov2,*/; //mixing matrix, AtN-1A, inverted AtN-1A, ?, SpN_all = PS for the bolo (ndet, nbins)
//	double *Nk; // Noise power spectrum
	double *ell; /**SpN,*/ /**Nell;*/ // ell = bins values, Nell = binned noise PS
	double **mixmat; // initialized to NULL to avoid warning, "mixmat may be used uninitialized" ...
	//	double *p, *uvec, *ivec;

	//fftw_plan fftplan;
//	fftw_complex *fdata1, *fdata2; // fourier transform of the data
	//fftw_complex *fdata_buffer ;
	//fdata_buffer= new fftw_complex[(ns/2+1)*ndet];

	time_t t1;


	printf("Inside EstimPowerSpectra just at the beginning\n");
	t1=time(NULL);

	//	data = new double[ns]; // raw data
	//	data_lp = new double[ns]; // data low passed
	//	Ps = new double[ns]; // prewhitened noise ?
	//	commontmp = new double[ns]; // useless
//	commonm_f = new double[ns]; // useless
	//flag = new unsigned char[ns]; // flag array
	//flag = new short[ns];
	//	samptopix = new long long[ns]; // sample to pixel proj matrix
//	Nk = new double[ns/2+1]; // noise PS
	//	bfilter = new double[ns/2+1]; // buttter filter values
	//	fdata1 = new fftw_complex[ns/2+1]; // fourier transform
//	fdata2 = new fftw_complex[ns/2+1];

	read_double(ellFile, ell, nbins);
	nbins = nbins-1;

//	Nell = new double[nbins]; // binned noise PS
	SPref = new double[nbins]; // first detector PS to avoid numerical problems
	P = dmatrix(0,ncomp-1,0,nbins-1); // component power spectra
	N = dmatrix(0,det.ndet-1,0,nbins-1); // uncorralated part of the noise
	Rellexp = dmatrix(0,(det.ndet)*(det.ndet)-1,0,nbins-1); // covariance du signal (mesurée)
	Rellth = dmatrix(0,(det.ndet)*(det.ndet)-1,0,nbins-1); // covariance matrix (theorique) = modele : inversible et pas Rellexp
	//mixmat = dmatrix(0,ndet-1,0,20);

	//sign = new double[ndet]; // sigma of the noise

	//	Cov = dmatrix(0,ncomp-1,0,ncomp-1); // AtN-1A
	//	iCov = dmatrix(0,ncomp-1,0,ncomp-1);  // inverted AtN-1A
	//	//iCov2 = dmatrix(0,ncomp-1,0,ncomp-1); // inutilisé
	//	p = new double[ncomp];
	//	uvec = new double[ncomp];
	//	ivec = new double[ncomp];

	vect = dmatrix(0,ncomp-1,0,det.ndet-1);


	commonm2 = dmatrix(0,ncomp,0,ns-1);



	init2D_double(commonm2,0,0,ncomp,ns,0.0);
	//	fill(commontmp,commontmp+ns,0.0);
//	fill(commonm_f,commonm_f+ns,0.0);
	//	init2D_double(Cov,0,0,ncomp,ncomp,0.0);
	//	init2D_double(iCov,0,0,ncomp,ncomp,0.0);


	//	for(long ii=0;ii<ns/2+1;ii++) // TODO : remettre dans les 2 prem fonctions
	//		bfilter[ii] = 1.0;

	//apodization
	apodwind = apodwindow(ns,int(ns*0.04));


	for (long ii=0;ii<ns;ii++)
		factapod += apodwind[ii]*apodwind[ii]/ns; // factapod ?? apodization factor ?

	//TODO: factapod is computed once, but used twice, it can be computed here (mat 28/10)
	//----------------------------------- READ MIXMAT PART -------------------------------//

	read_mixmat_file(MixMatfile, dir.dirfile, mixmat, det.ndet,ncomp); // TODO : la mixmat doit etre dans dirfile ou dans tmp_dir ??

	//----------------------------------- READ MIXMAT PART -------------------------------//


	// compute common mode commonm2

	common_mode_computation(det,proc_param, pos_param, dir, apodwind, ns, ff, NAXIS1, NAXIS2, npix, iframe, S, indpix,
			mixmat, ncomp, commonm2, factapod, fits_filename);

	//----------------------------------- ESTIMATE NOISE PS -------------------------------//

	estimate_noise_PS(det,  proc_param, pos_param, dir, nbins, nbins2, ns, ff, NAXIS1,
			NAXIS2, npix, ell, S, iframe,indpix, apodwind, ncomp, mixmat, commonm2,
			factapod,Rellth, N, P, fits_filename);

	//----------------------------------- ESTIMATE COVMAT of the DATA R_exp -------------------------------//

	estimate_CovMat_of_Rexp(dir, det, nbins, ns, ff, ell, ncomp, mixmat, proc_param.fsamp,
			factapod, Rellexp, N, P, SPref);

	//----------------------------------- FIT COMPONENT, PS and MIXMAT -------------------------------//


	expectation_maximization_algorithm(fcut, nbins, det.ndet, ncomp, ns, proc_param.fsamp, ff,
			dir.outdir,  Rellexp, Rellth, mixmat, P, N, SPref, ell);

	//----------------------------------- WRITE TO DISK -------------------------------//

	write_to_disk(dir.outdir, ff, det,nbins, ell, mixmat, Rellth,
			Rellexp, ncomp, N, SPref,P);
	//----------------------------------- END OF ESTIMPS -------------------------------//


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

}






