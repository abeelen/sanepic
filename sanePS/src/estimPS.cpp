/*
 * estimPS.cpp
 *
 *  Created on: 28 mai 2009
 *      Author: matthieu
 */

#include "estimPS.h"


using namespace std;


void EstimPowerSpectra(double fsamp, long ns, long ff, long ndet, int NAXIS1, int NAXIS2, long npix, long napod,
		long iframe, bool flgdupl, int factdupl, long *indpix,
		double *S, string MixMatfile, std::vector<string> bolonames, string dirfile,/* string bextension,
		string fextension,*/ int shift_data_to_point, string dir,
		bool NORMLIN, bool NOFILLGAP, bool remove_polynomia, string noiseSppreffile,
		string extentnoiseSp, string outdirSpN, string fits_filename)
{



	double fcut = 12; //fixed here to be sure that the code code will not focus on high freq
	// add it in the ini file ??
	long ncomp = 1;
	//long ncomp2 = 0;
	long nbins = 500;
	long nbins2; // readed bins in mixmatfile
	double factapod;
	//double sign0=1;


	// data, data low passed, Prewhitened signal ? ,apodization window, ?, ?, butter filter values, first detector PS (in the cross correlation estimation part)
	double *data, *data_lp, *Ps, /**calp, */*apodwind, *commontmp, *commonm_f, *bfilter, *SPref;
	long *samptopix; // sample to pixel projection matrix
	//unsigned char *flag; // flag array
//	short *flag;
	double /***commonm,*/ **commonm2,/* **common_f,*/ **vect; // ?
	double **P, **N, **Rellth, **Rellexp; // Rellth = covariance matrix
	double **Cov, **iCov,/* **iCov2,*/ **SpN_all; //mixing matrix, AtN-1A, inverted AtN-1A, ?, SpN_all = PS for the bolo (ndet, nbins)
	double *Nk; // Noise power spectrum
	double *ell, /**SpN,*/ *Nell; // ell = bins values, Nell = binned noise PS
	double **mixmat=NULL; // initialized to NULL to avoid warning, "mixmat may be used uninitialized" ...
	double *p, *uvec, *ivec;

	//fftw_plan fftplan;
	fftw_complex *fdata1, *fdata2; // fourier transform of the data
	//fftw_complex *fdata_buffer ;
	//fdata_buffer= new fftw_complex[(ns/2+1)*ndet];



	time_t t1;

	//double *data1d; // buffer used to write down 1d array


	//string testfile; // string used to remove sprintf horror
	//string nameSpfile; // noise Ps filename

	//std::ostringstream temp_stream; // used to remove sprintf horror

	//string field; // detector name in the loop
	//string tempstr1; // boloname : used to write down Rellth
	//string tempstr2; // boloname : used to write down Rellth


	printf("Inside EstimPowerSpectra just at the beginning\n");
	t1=time(NULL);

	data = new double[ns]; // raw data
	data_lp = new double[ns]; // data low passed
	Ps = new double[ns]; // prewhitened noise ?
	commontmp = new double[ns]; // useless
	commonm_f = new double[ns]; // useless
	//calp = new double[ns];
	//flag = new unsigned char[ns]; // flag array
	//flag = new short[ns];
	samptopix = new long[ns]; // sample to pixel proj matrix
	Nk = new double[ns/2+1]; // noise PS
	bfilter = new double[ns/2+1]; // buttter filter values
	fdata1 = new fftw_complex[ns/2+1]; // fourier transform
	fdata2 = new fftw_complex[ns/2+1];


	Nell = new double[nbins]; // binned noise PS
	SPref = new double[nbins]; // first detector PS to avoid numerical problems
	P = dmatrix(0,ncomp-1,0,nbins-1); // component power spectra
	N = dmatrix(0,ndet-1,0,nbins-1); // uncorralated part of the noise
	Rellexp = dmatrix(0,ndet*ndet-1,0,nbins-1); // covariance du signal (mesurée)
	Rellth = dmatrix(0,ndet*ndet-1,0,nbins-1); // covariance matrix (theorique) = modele : inversible et pas Rellexp
	//mixmat = dmatrix(0,ndet-1,0,20);

	//sign = new double[ndet]; // sigma of the noise

	Cov = dmatrix(0,ncomp-1,0,ncomp-1); // AtN-1A
	iCov = dmatrix(0,ncomp-1,0,ncomp-1);  // inverted AtN-1A
	//iCov2 = dmatrix(0,ncomp-1,0,ncomp-1); // inutilisé
	p = new double[ncomp];
	uvec = new double[ncomp];
	ivec = new double[ncomp];

	vect = dmatrix(0,ncomp-1,0,ndet-1);

	//commonm = dmatrix(0,ncomp,0,ns-1); // common mode
	commonm2 = dmatrix(0,ncomp,0,ns-1);
	//common_f = dmatrix(0,ncomp,0,ns-1); // useless

	//init2D_double(commonm,0,0,ncomp,ns,0.0);
	init2D_double(commonm2,0,0,ncomp,ns,0.0);
	//init2D_double(common_f,0,0,ncomp,ns,0.0);
	//init1D_double(commontmp,0,ns,0.0);
	//init1D_double(commonm_f,0,ns,0.0);
	fill(commontmp,commontmp+ns,0.0);
	fill(commonm_f,commonm_f+ns,0.0);
	init2D_double(Cov,0,0,ncomp,ncomp,0.0);
	init2D_double(iCov,0,0,ncomp,ncomp,0.0);
	//init2D_double(iCov2,0,0,ncomp,ncomp,0.0);

	for(long ii=0;ii<ns/2+1;ii++)
		bfilter[ii] = 1.0;

	//apodization
	apodwind = apodwindow(ns,int(ns*0.04));

	//----------------------------------- READ MIXMAT PART -------------------------------//

	read_mixmat_file(MixMatfile, dir, mixmat, ndet,ncomp);


	//----------------------------------- READ MIXMAT PART -------------------------------//


	// compute common mode commonm2
	common_mode_computation(apodwind, ndet, ns, ff, NAXIS1,NAXIS2, npix, flgdupl, factdupl, bolonames, /*bextension, fextension,*/
			dirfile,  shift_data_to_point,    dir, iframe, S, indpix,  NORMLIN,
			NOFILLGAP, remove_polynomia, napod, mixmat, ncomp, commonm2, samptopix, Ps, data, data_lp,
			bfilter, Cov, uvec, p, ivec, iCov, factapod, fdata1,fits_filename);


	//----------------------------------- ESTIMATE NOISE PS -------------------------------//

	estimate_noise_PS(bolonames, dirfile, extentnoiseSp, noiseSppreffile,/* bextension, fextension,*/ nbins,
			nbins2, ns, ff, ndet, NAXIS1, NAXIS2, npix,napod, ell, SpN_all, data,
			samptopix, dir, S, iframe,  Ps, data_lp, bfilter, indpix, NORMLIN,
			NOFILLGAP, remove_polynomia,flgdupl, factdupl, apodwind, ncomp, mixmat, commonm2, fsamp,
			Nk, Nell, factapod,Rellth, N, commontmp, P, shift_data_to_point,  outdirSpN,fits_filename);

	//----------------------------------- ESTIMATE COVMAT of the DATA R_exp -------------------------------//


	estimate_CovMat_of_Rexp(nbins, ns, ff, ndet, ell, dir,  ncomp, mixmat,fsamp,
			Nk, Nell, factapod, Rellexp, N, P,  outdirSpN, fdata1, fdata2, SPref);


	//----------------------------------- FIT COMPONENT, PS and MIXMAT -------------------------------//


	expectation_maximization_algorithm(fcut, nbins, ndet, ncomp,ns, fsamp, ff,
			outdirSpN,  Rellexp, Rellth, mixmat,P,N, Cov, p,	uvec, ivec, iCov, SPref, ell);



	//----------------------------------- WRITE TO DISK -------------------------------//
	//*****************  Write power spectra to disk  ********************//

	write_to_disk(outdirSpN, ff,  bolonames,nbins, ell, mixmat, Rellth,
			Rellexp, ncomp, ndet, N, SPref,P);
	//----------------------------------- END OF ESTIMPS -------------------------------//

	//delete [] data1d;

	delete [] data;
	delete [] data_lp;
	delete [] fdata1;
	delete [] fdata2;
	delete [] Ps;
	delete [] samptopix;
	delete [] commontmp;
	delete [] commonm_f;
	//delete [] calp;
	//delete [] flag;
	delete [] bfilter;
	delete [] apodwind;
	delete [] Nell;
	free_dmatrix(Rellexp,0,ndet*ndet-1,0,nbins-1);
	free_dmatrix(Rellth,0,ndet*ndet-1,0,nbins-1);
	free_dmatrix(mixmat,0,ndet-1,0,20);
	//delete [] sign;
	free_dmatrix(Cov,0,ncomp-1,0,ncomp-1);
	free_dmatrix(iCov,0,ncomp-1,0,ncomp-1);
	//free_dmatrix(iCov2,0,ncomp-1,0,ncomp-1);
	delete [] p;
	delete [] uvec;
	delete [] ivec;
	//free_dmatrix(commonm,0,ncomp,0,ns-1);
	free_dmatrix(commonm2,0,ncomp,0,ns-1);
	//free_dmatrix(common_f,0,ncomp,0,ns-1);
	free_dmatrix(vect,0,ncomp-1,0,ndet-1);
	delete [] Nk;
	delete [] ell;
	free_dmatrix(SpN_all,0,ndet-1,0,nbins-1);
	free_dmatrix(P,0,ncomp-1,0,nbins-1);
	free_dmatrix(N,0,ndet-1,0,nbins-1);



}






