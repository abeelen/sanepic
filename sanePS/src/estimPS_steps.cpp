#include <vector>
#include <sstream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sysexits.h>

#include "todprocess.h"
#include "map_making.h"
#include "temporary_IO.h"
#include "inputFileIO.h"
#include "dataIO.h"
#include "covMatrix_IO.h"
#include "estimPS_steps.h"
#include "cholesky.h"

extern "C" {
#include "nrutil.h"
#include <fitsio.h>
}


using namespace std;

int common_mode_computation(struct samples samples_struct, std::vector<std::string> det, struct param_saneProc proc_param, struct param_sanePos pos_param,
		struct param_common dir, double *apodwind,long ns, long NAXIS1, long NAXIS2, long long npix,
		double *S, long long *indpix,double **mixmat, long ncomp, double **commonm2,
		double &factapod, string fits_filename)
{
	//*************************** Read data and compute components

	string field; // detector name in the loop

	double *sign;
	double tmpsign, mm;
	double sign0=1;
	double **commonm;
	int *flag;

	fftw_plan fftplan;


	double *data, *data_lp, *Ps/*, *bfilter*/;
	long long *samptopix; // sample to pixel projection matrix
	double **iCov, **Cov, *ivec, **l;
	double *uvec;

	long ndet = (long)det.size();

	fftw_complex *fdata1;

	int factdupl = 1;
	if(pos_param.flgdupl==1) factdupl = 2;


	data    = new double[ns];
	flag    = new int[ns];
	data_lp = new double[ns]; // data low passed
	Ps = new double[ns]; // deprojected signal
	samptopix = new long long[ns]; // sample to pixel proj matrix
	fdata1 = new fftw_complex[ns/2+1]; // fourier transform


	Cov = dmatrix(0,ncomp-1,0,ncomp-1); // AtN-1A
	iCov = dmatrix(0,ncomp-1,0,ncomp-1);  // inverted AtN-1A
	uvec = new double[ncomp];
	ivec = new double[ncomp];
	l    = new double*[ncomp];
	for(int i=0; i<ncomp; i++)
		l[i]=new double [ncomp];


	init2D_double(Cov,0,0,ncomp,ncomp,0.0);
	init2D_double(iCov,0,0,ncomp,ncomp,0.0);

	sign = new double[ndet];
	commonm = dmatrix(0,ncomp,0,ns-1); // common mode
	init2D_double(commonm,0,0,ncomp,ns,0.0);

	// loop over detectors
	for (long idet=0;idet<ndet;idet++){

		field = det[idet];

		if(read_data_from_dirfile(samples_struct.dirfile_pointer, fits_filename, field, data, ns))
			return 1;
		if(read_flag_from_dirfile(samples_struct.dirfile_pointer, fits_filename, field, flag, ns))
			return 1;

		if (S != NULL){
			//Read pointing data
			if(read_samptopix(samples_struct.dirfile_pointer, ns, &samptopix, fits_filename, field))
				return EX_NOINPUT;

			deproject(S,indpix,samptopix,ns,NAXIS1,NAXIS2,npix,Ps,pos_param.flgdupl,factdupl);

			for(long ii=0;ii<ns;ii++)
				data[ii] = data[ii] - Ps[ii];
		}

		//TODO : f_lp_pix is hard fixed to 1.0 -> avoid errors and wasted time ?
		MapMakePreProcessData(data,flag,ns,proc_param ,1.0,data_lp, NULL);


		// should apodisation be part of MapMakePreProcess ? : no, not in Corr_preprocess !
		for (long ii=0;ii<ns;ii++)
			data[ii] = data_lp[ii]*apodwind[ii];

		//       fdata are used in cross power spectrum estimation...
		//       BUT it is done differently than power spectrum estimation WHY ?

		// compute fft and save data to disk for later
		fftplan = fftw_plan_dft_r2c_1d(ns, data, fdata1, FFTW_ESTIMATE);
		fftw_execute(fftplan);
		fftw_destroy_plan(fftplan);


		if(write_fdata(samples_struct.dirfile_pointer, ns, fdata1,  "fdata_", idet, fits_filename, det))
			return EX_DATAERR;

		/// compute sigma of the noise
		mm = 0.0;
		for (long ii=ns/2;ii<ns/2+500;ii++) mm += data[ii]; // sur 500 samples seulement ??
		mm = mm/501.0;
		tmpsign = 0.0;
		for (long ii=ns/2;ii<ns/2+500;ii++) tmpsign += (data[ii]-mm)*(data[ii]-mm);
		sign[idet] = sqrt(tmpsign/500.0);
		// normalize to the first detector
		if (idet == 0) sign0 = sign[0];
		sign[idet] = sign[idet]/sign0; // normalize to first detector sigma of the noise


		// common mode computation
		for (long jj=0;jj<ncomp;jj++)
			for (long ii=0;ii<ns;ii++)
				commonm[jj][ii] += mixmat[idet][jj]/(sign[idet]*sign[idet])*data[ii];

	}


	delete [] data;
	delete [] flag;

	for (long jj=0; jj<ncomp; jj++)
		for (long ii= 0 ;ii<ns;ii++)
			if ( isnan(commonm[jj][ii]) || isinf(commonm[jj][ii]) ){
				cout << fits_filename << " something went wrong in the computation of the common mode " << ncomp << endl;
				exit(EXIT_FAILURE);
			}

	//***************************************************************************


	/////////// AtN-1A
	for (long jj=0;jj<ncomp;jj++){
		for (long kk=0;kk<ncomp;kk++){
			for (long ii=0;ii<ndet;ii++){
				Cov[jj][kk] += mixmat[ii][jj] * mixmat[ii][kk]/sign[ii]/sign[ii];
			}
		}
	}


	// invert AtN-1A

	cholesky(ncomp,Cov,l);

	//	printf("ncomp:%ld", ncomp);

	for (long ii=0;ii<ncomp;ii++){
		for (long jj=0;jj<ncomp;jj++)
			uvec[jj] = 0.0;
		uvec[ii] = 1.0;
		solve_cholesky(Cov, uvec, l, ivec, ncomp);
		for (long jj=0;jj<ncomp;jj++)
			iCov[ii][jj] = ivec[jj];
	}

	//	printf("noise var det 0 =  %10.15g\n",sign0*sign0);


	for (long ii=0;ii<ns;ii++)
		for (long jj=0;jj<ncomp;jj++)
			for (long kk=0;kk<ncomp;kk++)
				commonm2[jj][ii] += iCov[jj][kk] * commonm[kk][ii]; //common mode * (AtN-1A)-1

	// clean up
	delete [] sign;
	delete [] data_lp;
	delete [] Ps;
	delete [] samptopix;
	delete [] fdata1;
	delete [] uvec;
	delete [] ivec;

	for (int i = 0; i <ncomp; i++)
		delete[] l[i];
	delete [] l;

	free_dmatrix(Cov,0,ncomp-1,0,ncomp-1);
	free_dmatrix(iCov,0,ncomp-1,0,ncomp-1);
	free_dmatrix(commonm,0,ncomp,0,ns-1);

	return EX_OK;
}



int estimate_noise_PS(struct samples samples_struct, std::vector<std::string> det, struct param_saneProc proc_param,struct param_sanePos pos_param,
		struct param_common dir, long &nbins,	long &nbins2, long ns, long NAXIS1,
		long NAXIS2, long long npix, double *&ell, double *S,long long *indpix,
		double *apodwind, long ncomp, double **mixmat, double **commonm2,
		double factapod,double **Rellth, double **N, double **P, string fits_filename)
{


	string nameSpfile, field;
	string testfile;
	string basename;
	std::ostringstream temp_stream;
	FILE *fp;

	int *flag;

	double *data, *data_lp, *Ps=NULL;
	double  *commontmp;
	double *Nell, *Nk;
	long long *samptopix; // sample to pixel projection matrix

	long ndet = (long)det.size();

	int factdupl = 1;
	if(pos_param.flgdupl==1) factdupl = 2; // map duplication factor

	data      = new double[ns];
	flag      = new int[ns];
	data_lp   = new double[ns]; // data low passed
	commontmp = new double[ns]; //
	Nell      = new double[nbins]; // binned noise PS
	Nk        = new double[ns/2+1]; // noise PS

	if (S != NULL){
		samptopix = new long long[ns]; // sample to pixel proj matrix
		Ps = new double[ns];
	}


	fill(commontmp,commontmp+ns,0.0);

	//----------------------------------- ESTIMATE NOISE PS -------------------------------//


	//************************************************************************//
	// second part: -- data - common mode
	//              -- estimate noise power spectra


	/////////////////////////////////////// loop again over detectors
	for (long idet=0;idet<ndet;idet++){

		field = det[idet];

		if(read_data_from_dirfile(samples_struct.dirfile_pointer, fits_filename, field, data, ns))
			return 1;
		if(read_flag_from_dirfile(samples_struct.dirfile_pointer, fits_filename, field, flag, ns))
			return 1;

		//TODO : This computation is already done when computing the common mode ->
		//       reuse the fdata if possible?
		//******************************* subtract signal

		if (S != NULL){
			//Read pointing data
			if(read_samptopix(samples_struct.dirfile_pointer, ns, &samptopix, fits_filename, field))
				return EX_NOINPUT;

			deproject(S,indpix,samptopix,ns,NAXIS1,NAXIS2,npix,Ps,pos_param.flgdupl,factdupl);

			for(long ii=0;ii<ns;ii++)
				data[ii] = data[ii] - Ps[ii];

		}

		//TODO : why f_lppix set to 1.0 ?
		MapMakePreProcessData(data,  flag, ns, proc_param, 1.0, data_lp, NULL);

		for (long ii=0;ii<ns;ii++)
			data[ii] = data_lp[ii] * apodwind[ii];


		// Subtract components
		for (long ii=0;ii<ns;ii++)
			for (long jj=0;jj<ncomp;jj++)
				data[ii] -= mixmat[idet][jj]*commonm2[jj][ii];


		// Noise power spectra
		/// measure power spectrum of the uncorrelated part of the noise
		noisepectrum_estim(data,ns,ell,(int)nbins,proc_param.fsamp,NULL,Nell,Nk);

		//TODO : normalization by factapod is also done in noisespectrum_estim ?? DONE TWICE ??

		for (long ii=0;ii<nbins;ii++){
			Rellth[idet*ndet+idet][ii] += Nell[ii]/factapod; // uncorrelated part added in covariance matrix ??
			N[idet][ii] = Nell[ii]/factapod; // uncorrelated part
		}

	}

	delete [] data;
	delete [] flag;

	////*********************** Component power spectra

	for (long ii=0;ii<ncomp;ii++){
		for (long jj=0;jj<ns;jj++)
			commontmp[jj]=commonm2[ii][jj];
		noisepectrum_estim(commontmp,ns,ell,(int)nbins,proc_param.fsamp,NULL,Nell,Nk);
		for (long jj=0;jj<nbins;jj++)
			P[ii][jj] = Nell[jj]/factapod;

		basename = FitsBasename(fits_filename);

		temp_stream << dir.output_dir + "Nellc_" << ii << "_" << basename << ".bi";
		// Get the string
		testfile= temp_stream.str();
		// Clear ostringstream buffer
		temp_stream.str("");

		if((fp = fopen(testfile.c_str(),"w"))){
			fwrite(Nell,sizeof(double), nbins, fp);
			fclose(fp);
		}else{
			cerr << "Error. Can't open " << testfile << ".Exiting.\n";
			return 1;
		}

	}

	for (long ii=0;ii<ndet;ii++)
		for (long kk=0;kk<ndet;kk++)
			for (long ll=0;ll<ncomp;ll++)
				for (long jj=0;jj<nbins;jj++)
					Rellth[ii*ndet+kk][jj] += mixmat[ii][ll] * mixmat[kk][ll] * P[ll][jj]; // add correlated part to covariance matrix


	// clean up
	if (S != NULL){
		delete [] Ps ;
		delete [] samptopix;
	}

	delete [] data_lp ;
	delete [] commontmp;
	delete [] Nell;
	delete [] Nk;

	return EX_OK;
}


int estimate_CovMat_of_Rexp(struct samples samples_struct, struct param_common dir, std::vector<std::string> det, long nbins, long ns, double *ell, long ncomp, double **mixmat,double fsamp,
		double factapod,double **Rellexp, double **N, double **P, double *SPref, string fits_filename, int rank)
{

	long ndet = (long)det.size();

	std::ostringstream temp_stream; // used to remove sprintf horror

	double *data1d; // buffer used to write down 1d array
	data1d = new double[ndet*nbins];
	string testfile;
	string basename;
	fftw_complex *fdata1, *fdata2;
	double * Nell, *Nk;


	Nell = new double[nbins]; // binned noise PS
	Nk = new double[ns/2+1]; // noise PS

	fdata1 = new fftw_complex[ns/2+1]; // fourier transform
	fdata2 = new fftw_complex[ns/2+1];

	FILE *fp;

	/////////////////////////////////////// loop again over detectors
	for (long idet1=0;idet1<ndet;idet1++){

		// read data from disk
		if(read_fdata(samples_struct.dirfile_pointer, ns, &fdata1, "fdata_", idet1, fits_filename, det))
			return EX_NOINPUT;

		for (long idet2=0;idet2<ndet;idet2++) {

			// read data from disk
			if(read_fdata(samples_struct.dirfile_pointer, ns, &fdata2, "fdata_", idet2, fits_filename, det))
				return EX_NOINPUT;

			noisecrosspectrum_estim(fdata1,fdata2,ns,ell,(int)nbins,fsamp,NULL,Nell,Nk);


			for (long ii=0;ii<nbins;ii++)
				Rellexp[idet1*ndet+idet2][ii] += Nell[ii]/factapod; // noise cross PS ?

		}
#ifdef DEBUG
		for(int rk=0;rk<rank;rk++)
			cout << "\t\t\t"; // try to deal with MPI screen outputs
		cout << "[ " << rank << " ]" << " Rellexp :" << setprecision(2) << idet1*100./ndet << "%\r" << flush ;
#endif
	}

	////// normalize to the first detector power spectrum in order to avoid numerical problems
	for (long ii=0;ii<nbins;ii++)
		SPref[ii] = Rellexp[0][ii]; // first detector PS

	for (long ii=0;ii<nbins;ii++)
		if ( isnan(SPref[ii]) || isinf(SPref[ii])){
			cout << " Problem in the first detector power spectrum\n";
			return 1;
		}
	for (long idet1=0;idet1<ndet;idet1++)
		for (long idet2=0;idet2<ndet;idet2++)
			for (long ii=0;ii<nbins;ii++)
				Rellexp[idet1*ndet+idet2][ii] = Rellexp[idet1*ndet+idet2][ii]/SPref[ii]; // normalize to first detector
	for (long jj=0;jj<ncomp;jj++)
		for (long ii=0;ii<nbins;ii++)
			P[jj][ii] = P[jj][ii]/SPref[ii]; // normalize common mode part
	for (long jj=0;jj<ndet;jj++)
		for (long ii=0;ii<nbins;ii++)
			N[jj][ii] = N[jj][ii]/SPref[ii]; //normalize uncorrelated part


	basename = FitsBasename(fits_filename);

	// write Rellexp to disk and also first guess of parameters
	temp_stream << dir.output_dir + "Rellexp_" << basename << ".txt";

	// get filename
	testfile= temp_stream.str();
	temp_stream.str("");

	fp = fopen(testfile.c_str(),"w");
	for (long jj=0;jj<nbins;jj++)
		for (long ii=0;ii<ndet;ii++)
			for (long kk=0;kk<ndet;kk++)
				fprintf(fp,"%10.15g \t",Rellexp[ii*ndet+kk][jj]); // cross power spectrum
	fprintf(fp,"\n");
	fclose(fp);

	temp_stream << dir.output_dir + "Ninit_" << basename << ".txt";

	// get filename
	testfile= temp_stream.str();
	temp_stream.str("");

	fp = fopen(testfile.c_str(),"w");
	for (long ii=0;ii<ndet;ii++){
		for (long jj=0;jj<nbins;jj++)
			fprintf(fp,"%10.15g \t",N[ii][jj]); // uncorralated part of the noise
		fprintf(fp,"\n");
	}
	fclose(fp);


	for (long i=0; i< ndet; i++)
		for (long j=0; j<nbins; j++)
			data1d[i*nbins+j] = N[i][j];

	temp_stream << "!" + dir.output_dir + "Ninit_" << basename << ".fits";

	// get filename
	testfile= temp_stream.str();
	temp_stream.str("");

	write_psd_tofits(testfile.c_str(),ndet,nbins,'d', data1d); //resized uncorralated part

	temp_stream << dir.output_dir + "Pinit_" << basename << ".txt";

	// get filename
	testfile= temp_stream.str();
	temp_stream.str("");

	fp = fopen(testfile.c_str(),"w");
	for (long ii=0;ii<ncomp;ii++)
		for (long jj=0;jj<nbins;jj++)
			fprintf(fp,"%10.15g \t",P[ii][jj]); // common mode part of the noise
	fprintf(fp,"\n");
	fclose(fp);


	temp_stream << dir.output_dir + "Ainit_" << basename << ".txt";

	// get filename
	testfile= temp_stream.str();
	temp_stream.str("");

	fp = fopen(testfile.c_str(),"w");
	for (long ii=0;ii<ndet;ii++)
		for (long jj=0;jj<ncomp;jj++){
			fprintf(fp,"%10.15g \t",mixmat[ii][jj]); // mixing matrix
			fprintf(fp,"\n");
		}
	fclose(fp);



	delete [] data1d;
	delete [] fdata1;
	delete [] fdata2;
	delete [] Nell;
	delete [] Nk;

	return EX_OK;
}


int expectation_maximization_algorithm(double fcut, long nbins, long ndet, long ncomp,long ns, double fsamp,
		string outdirSpN,	double **Rellexp, double **Rellth, double **mixmat,double **P,double **N,
		double *SPref, double *ell, int rank)
{


	//***** Fourth part
	//*********************** fit component and noise power spectra, and mixing matrix *************//
	//********* Using Expectation/Maximization algorithm
	long nbiter = 500;
	long ib=0; // used as an iterator

	double tottest=0.0;

	double f;

	double *iN, *Pr, *w;
	double **Rxs, **Rxsq, **RnRxsb, **Rxx, **Rxxq, **Rss, **Rssq, **RnRssb;
	double **Pr2, **AiNA, **Mattmp, **ImDR, **ACq, **Cq, **Wq;
	long nbins2=0;

	while ((ell[ib] < fcut) && (ib < nbins)){
		nbins2 = ib+1;
		ib++;
	}



	double **Cov, **iCov, **l;
	double *uvec, *ivec;


	Cov = dmatrix(0,ncomp-1,0,ncomp-1); // AtN-1A
	iCov = dmatrix(0,ncomp-1,0,ncomp-1);  // inverted AtN-1A
	uvec = new double[ncomp];
	ivec = new double[ncomp];
	l = new double*[ncomp];
	for(int i=0; i<ncomp;i++)
		l[i]=new double[ncomp];


	init2D_double(Cov,0,0,ncomp,ncomp,0.0);
	init2D_double(iCov,0,0,ncomp,ncomp,0.0);

//	printf("\nnbins2 = %ld                \n",nbins2);



	iN = new double[ndet];
	Pr = new double[ncomp];
	w = new double[nbins2];

	Rxs =    dmatrix(0,ndet-1 ,0,ncomp-1);
	Rxsq =   dmatrix(0,ndet-1 ,0,ncomp-1);
	RnRxsb = dmatrix(0,ndet-1 ,0,ncomp-1);
	Rxx =    dmatrix(0,ndet-1 ,0,ndet-1) ;
	Rxxq =   dmatrix(0,ndet-1 ,0,ndet-1) ;
	Rss =    dmatrix(0,ncomp-1,0,ncomp-1);
	Rssq =   dmatrix(0,ncomp-1,0,ncomp-1);
	RnRssb = dmatrix(0,ncomp-1,0,ncomp*ndet-1);
	Pr2 =    dmatrix(0,ncomp-1,0,ncomp-1);
	AiNA =   dmatrix(0,ncomp-1,0,ncomp-1);
	Mattmp = dmatrix(0,ndet-1 ,0,ndet-1) ;
	ImDR =   dmatrix(0,ndet-1, 0,ndet-1) ;
	ACq =    dmatrix(0,ndet-1, 0,ncomp-1);
	Cq =     dmatrix(0,ncomp-1,0,ncomp-1);
	Wq =     dmatrix(0,ndet-1, 0,ncomp-1);


	//// Compute weights
	for (long ii=0;ii<nbins2;ii++)
		w[ii] = (ell[ii+1] - ell[ii])*ns/fsamp;


	f = fdsf(Rellexp,w,mixmat,P,N,ndet,ncomp,nbins2) ;

#ifdef DEBUG
	printf("Pre em:   obj: %10.15g\n", f) ;

	cout << "nbins2 " << nbins2 << endl;
	cout << "ndet   " << ndet << endl;
	cout << "ncomp  " << ncomp << endl << endl;
#endif

	for (long iter=1;iter<=nbiter;iter++){

		fill(iN,iN+ndet,0.0);
		fill(Pr,Pr+ncomp,0.0);
		init2D_double(Rxs,0,0,ndet,ncomp,0.0);
		init2D_double(Rxsq,0,0,ndet,ncomp,0.0);
		init2D_double(RnRxsb,0,0,ndet,ncomp,0.0);
		init2D_double(Rxx,0,0,ndet,ndet,0.0);
		init2D_double(Rxxq,0,0,ndet,ndet,0.0);
		init2D_double(Rss,0,0,ncomp,ncomp,0.0);
		init2D_double(Rssq,0,0,ncomp,ncomp,0.0);
		init2D_double(RnRssb,0,0,ncomp,ncomp*ndet,0.0);
		init2D_double(Mattmp,0,0,ndet,ndet,0.0);
		init2D_double(ImDR,0,0,ndet,ndet,0.0);
		init2D_double(ACq,0,0,ndet,ncomp,0.0);
		init2D_double(Cq,0,0,ncomp,ncomp,0.0);
		init2D_double(Wq,0,0,ndet,ncomp,0.0);

		for (long ib=0;ib<nbins2;ib++){

			for (long idet=0;idet<ndet;idet++)
				iN[idet] = 1.0/N[idet][ib];

			for (long idet1=0;idet1<ndet;idet1++)
				for (long idet2=0;idet2<ndet;idet2++)
					Rxxq[idet1][idet2] = Rellexp[idet1*ndet + idet2][ib];

			// Robust wrt Pq=0
			for (long ii=0;ii<ncomp;ii++)
				Pr[ii] = sqrt(P[ii][ib]);

			for (long ii=0;ii<ncomp;ii++)
				for (long jj=0;jj<ncomp;jj++)
					Pr2[ii][jj] = Pr[ii]*Pr[jj];

			for (long ii=0;ii<ncomp;ii++)
				for (long jj=0;jj<ncomp;jj++){
					AiNA[ii][jj] = 0.0;
					for (long idet=0;idet<ndet;idet++)
						AiNA[ii][jj] += mixmat[idet][jj] * mixmat[idet][ii] * iN[idet];
				}

			for (long ii=0;ii<ncomp;ii++){
				for (long jj=0;jj<ncomp;jj++){
					Cov[ii][jj] = Pr2[ii][jj] * AiNA[ii][jj];
				}
				Cov[ii][ii] += 1.0;
			}

			// invert matrix
			cholesky(ncomp, Cov,l);

			for (long ii=0;ii<ncomp;ii++){
				for (long jj=0;jj<ncomp;jj++)
					uvec[jj] = 0.0;
				uvec[ii] = 1.0;
				solve_cholesky(Cov, uvec, l, ivec, ncomp);

				for (long jj=0;jj<ncomp;jj++){
					iCov[ii][jj] = ivec[jj];
				}

			}


			for (long ii=0;ii<ncomp;ii++)
				for (long jj=0;jj<ncomp;jj++)
					Cq[ii][jj] = Pr2[ii][jj] * iCov[ii][jj];

			for (long idet=0;idet<ndet;idet++)
				for (long ii=0;ii<ncomp;ii++){
					Wq[idet][ii] = 0.0;
					for (long jj=0;jj<ncomp;jj++)
						Wq[idet][ii] += mixmat[idet][jj] * iN[idet] *Cq[jj][ii] ;
				}
			for (long idet=0;idet<ndet;idet++)
				for (long ii=0;ii<ncomp;ii++){
					Rxsq[idet][ii] = 0.0;
					for (long jj=0;jj<ndet;jj++)
						Rxsq[idet][ii] += Rxxq[idet][jj] * Wq[jj][ii];
				}

			for (long kk=0;kk<ncomp;kk++)
				for (long ii=0;ii<ncomp;ii++){
					Rssq[kk][ii] = Cq[kk][ii];
					for (long jj=0;jj<ndet;jj++)
						Rssq[kk][ii] += Rxsq[jj][kk] * Wq[jj][ii];
				}

			for (long kk=1;kk<ncomp;kk++)
				for (long ii=0;ii<kk;ii++){
					Rssq[ii][kk] = 0.5*(Rssq[ii][kk]+Rssq[kk][ii]);
					Rssq[kk][ii] = Rssq[ii][kk];
				}


			// update power spectra
			for (long ii=0;ii<ncomp;ii++)
				P[ii][ib] = abs(Rssq[ii][ii]);



			for (long ii=0;ii<ndet;ii++)
				for (long jj=0;jj<ncomp;jj++)
					RnRxsb[ii][jj] += w[ib] * iN[ii]*Rxsq[ii][jj];



			for (long kk=0;kk<ncomp;kk++)
				for (long ii=0;ii<ndet;ii++)
					for (long jj=0;jj<ncomp;jj++)
						RnRssb[kk][jj+ii*ncomp] = RnRssb[kk][jj+ii*ncomp] + w[ib] * iN[ii] * Rssq[kk][jj] ;

		}



		// update mixing matrix
		for (long idet=0;idet<ndet;idet++){

			for (long ii=0;ii<ncomp;ii++){
				uvec[ii] = RnRxsb[idet][ii];
				for (long jj=0;jj<ncomp;jj++){
					Cov[ii][jj] = RnRssb[ii][jj+idet*ncomp];
				}
			}

			// solving the linear system

			cholesky(ncomp, Cov, l);
			solve_cholesky(Cov, uvec, l, ivec, ncomp);

			for (long ii=0;ii<ncomp;ii++){
				mixmat[idet][ii] = ivec[ii];
			}
		}

		// EM Step with respect to N, with the new values of A and P
		for (long ib=0;ib<nbins2;ib++){

			for (long idet = 0;idet<ndet;idet++)
				iN[idet] = 1.0/N[idet][ib];

			for (long idet1=0;idet1<ndet;idet1++)
				for (long idet2=0;idet2<ndet;idet2++)
					Rxxq[idet1][idet2] = Rellexp[idet1*ndet + idet2][ib];

			// Robust wrt Pq=0
			for (long ii=0;ii<ncomp;ii++)
				Pr[ii] = sqrt(P[ii][ib]);

			for (long ii=0;ii<ncomp;ii++)
				for (long jj=0;jj<ncomp;jj++)
					Pr2[ii][jj] = Pr[ii]*Pr[jj];

			for (long ii=0;ii<ncomp;ii++)
				for (long jj=0;jj<ncomp;jj++){
					AiNA[ii][jj] = 0.0;
					for (long idet=0;idet<ndet;idet++)
						AiNA[ii][jj] += mixmat[idet][jj] * mixmat[idet][ii] * iN[idet];
				}

			for (long ii=0;ii<ncomp;ii++){
				for (long jj=0;jj<ncomp;jj++){
					Cov[ii][jj] = Pr2[ii][jj] * AiNA[ii][jj];
				}
				Cov[ii][ii] += 1.0;
			}



			// invert matrix

			cholesky(ncomp, Cov, l);

			for (long ii=0;ii<ncomp;ii++){
				for (long jj=0;jj<ncomp;jj++)
					uvec[jj] = 0.0;
				uvec[ii] = 1.0;
				solve_cholesky(Cov, uvec, l, ivec, ncomp);

				for (long jj=0;jj<ncomp;jj++)
					iCov[ii][jj] = ivec[jj];
			}


			for (long ii=0;ii<ncomp;ii++)
				for (long jj=0;jj<ncomp;jj++)
					Cq[ii][jj] = Pr2[ii][jj] * iCov[ii][jj];

			for (long idet=0;idet<ndet;idet++)
				for (long ii=0;ii<ncomp;ii++){
					ACq[idet][ii] = 0.0;
					for (long jj=0;jj<ncomp;jj++)
						ACq[idet][ii] += mixmat[idet][jj]*Cq[jj][ii];
				}

			for (long idet=0;idet<ndet;idet++)
				for (long jj=0;jj<ndet;jj++){
					ImDR[jj][idet] = 0.0;
					if (jj == idet)
						ImDR[jj][idet] = 1.0;
					for (long ii=0;ii<ncomp;ii++)
						ImDR[jj][idet] -= ACq[jj][ii]*iN[idet]*mixmat[idet][ii] ;
				}


			for (long idet=0;idet<ndet;idet++)
				for (long ii=0;ii<ndet;ii++){
					Mattmp[idet][ii] = 0.0;
					for (long jj=0;jj<ndet;jj++)
						Mattmp[idet][ii] += Rellexp[idet+jj*ndet][ib] * ImDR[ii][jj];
				}

			for (long idet=0;idet<ndet;idet++){
				N[idet][ib] = 0.0;
				for (long ii=0;ii<ndet;ii++)
					N[idet][ib] += ImDR[idet][ii] * Mattmp[ii][idet];
			}


			for (long idet=0;idet<ndet;idet++){
				for (long ii=0;ii<ncomp;ii++){
					N[idet][ib] += ACq[idet][ii]*mixmat[idet][ii];
				}
				N[idet][ib] = abs(N[idet][ib]);
			}


		}


		tottest = 0.0;
		for (long ib=0;ib<nbins2;ib++){
			for (long idet=0;idet<ndet;idet++)
				tottest += N[idet][ib]/ndet;
			for (long idet=0;idet<ndet;idet++)
				if (N[idet][ib] < tottest*1e-8)
					N[idet][ib] = tottest*1e-8;
		}


		///// here is the problem

		f = fdsf(Rellexp,w,mixmat,P,N,ndet,ncomp,nbins2) ;

#ifdef DEBUG
		for(int rk=0;rk<rank;rk++)
			cout << "\t\t\t\t"; // try to deal with MPI screen outputs
		cout << "[ " << rank << " ]" << " em->iter: " << setprecision(2) << iter*100./nbiter << " %\r" << flush;
#endif
		if (isnan(f) || isinf(f)) {
			cout << "Nan........." << endl;
			return EX_SOFTWARE;
		}
	}


	// Fixing the indeterminacies.  Is it useful here?
	rescaleAP(mixmat, P, ndet, ncomp, nbins2) ;

	//****************************** Compute covariance matrix from the fitted model

#ifdef DEBUG
	cout << "\n[ " << rank << " ] EM step completed\n";
#endif

	for (long jj=0;jj<nbins;jj++)
		for (long ii=0;ii<ndet*ndet;ii++)
			Rellth[ii][jj] = 0.0;

	for (long jj=0;jj<nbins2;jj++){
		for (long idet=0;idet<ndet;idet++){
			Rellth[idet*ndet+idet][jj] += N[idet][jj]*SPref[jj];
			for (long ii=0;ii<ndet;ii++)
				for (long ll=0;ll<ncomp;ll++)
					Rellth[idet*ndet+ii][jj] += mixmat[idet][ll] * mixmat[ii][ll] * P[ll][jj]*SPref[jj];
		}
	}


	if (nbins2 < nbins)
		for (long jj=nbins2;jj<nbins;jj++)
			for (long idet=0;idet<ndet;idet++)
				Rellth[idet*ndet+idet][jj] = Rellexp[idet*ndet+idet][jj]*SPref[jj];


	//cleaning up

	delete [] uvec ;
	delete [] ivec;
	delete [] w;
	for(int i=0; i<ncomp; i++)
		delete [] l[i];
	delete [] l;

	free_dmatrix(Cov,0,ncomp-1,0,ncomp-1);
	free_dmatrix(iCov,0,ncomp-1,0,ncomp-1);

	delete [] iN;
	delete [] Pr;
	free_dmatrix(Rxs,0,ndet-1 ,0,ncomp-1);
	free_dmatrix(Rxsq,0,ndet-1 ,0,ncomp-1);
	free_dmatrix(RnRxsb,0,ndet-1 ,0,ncomp-1);
	free_dmatrix(Rxx,0,ndet-1 ,0,ndet-1);
	free_dmatrix(Rxxq,0,ndet-1 ,0,ndet-1);
	free_dmatrix(Rss,0,ncomp-1,0,ncomp-1);
	free_dmatrix(Rssq,0,ncomp-1,0,ncomp-1);
	free_dmatrix(RnRssb,0,ncomp-1,0,ncomp*ndet-1);
	free_dmatrix(Pr2,0,ncomp-1,0,ncomp-1);
	free_dmatrix(AiNA,0,ncomp-1,0,ncomp-1);
	free_dmatrix(Mattmp,0,ndet-1,0,ndet-1);
	free_dmatrix(ImDR,0,ndet-1,0,ndet-1);
	free_dmatrix(ACq,0,ndet-1,0,ncomp-1);
	free_dmatrix(Cq,0,ncomp-1,0,ncomp-1);
	free_dmatrix(Wq,0,ndet-1,0,ncomp-1);

	return EX_OK;
}


double fdsf(double **Rellexp, double *w, double **A, double **P, double **N, long ndet, long ncomp, long nbins){

	double f;

	double triRhR, logdetiR;

	double *uvec, *ivec;
	double *Pl, *Pnl;
	double **R, **hR, **eR, **iR, **iRhR, **l;

	uvec = new double[ndet];
	ivec = new double[ndet];
	l= new double*[ndet];
	for(int i=0; i<ndet; i++)
		l[i]=new double[ndet];

	Pl   = new double[ncomp] ;
	Pnl  = new double[ndet] ;
	R  = dmatrix(0,ndet-1,0,ndet-1);
	hR = dmatrix(0,ndet-1,0,ndet-1);
	eR = dmatrix(0,ndet-1,0,ndet-1);
	iR = dmatrix(0,ndet-1,0,ndet-1);
	iRhR = dmatrix(0,ndet-1,0,ndet-1);


	init2D_double(R,0,0,ndet,ndet,0.0);
	init2D_double(hR,0,0,ndet,ndet,0.0);
	init2D_double(eR,0,0,ndet,ndet,0.0);
	init2D_double(iR,0,0,ndet,ndet,0.0);
	fill(Pl,Pl+ncomp,0.0);
	fill(Pnl,Pnl+ndet,0.0);
	init2D_double(iRhR,0,0,ndet,ndet,0.0);


	// init
	f   = 0.0 ;

	for (long ib=0;ib<nbins;ib++){

		//// reading Rexp
		for (long ii=0;ii<ncomp;ii++)
			Pl[ii] = P[ii][ib];
		for (long ii=0;ii<ndet;ii++){
			Pnl[ii] = N[ii][ib];
			for (long jj=0;jj<ndet;jj++)
				hR[ii][jj] = Rellexp[ii+jj*ndet][ib];
		}

		/// building Rth from A, P, N
		for (long ii=0;ii<ndet;ii++)
			for (long jj=0;jj<ndet;jj++){
				R[ii][jj] = 0.0;
				if (ii == jj)
					R[ii][jj] = Pnl[ii];
				for (long kk=0;kk<ncomp;kk++)
					R[ii][jj] += A[ii][kk] * Pl[kk] * A[jj][kk];
			}


		for (long ii=0;ii<ndet;ii++){
			for (long jj=0;jj<ndet;jj++)
				eR[ii][jj] = R[ii][jj];
		}
		cholesky(ndet,eR,l);

		for (long ii=0;ii<ndet;ii++){
			for (long jj=0;jj<ndet;jj++)
				uvec[jj] = 0.0;
			uvec[ii] = 1.0;
			solve_cholesky(eR, uvec,l, ivec, ndet);

			for (long jj=0;jj<ndet;jj++)
				iR[ii][jj] = ivec[jj];

		}

		/// computing mismatch from Rexp and Rth
		for (long ii=0;ii<ndet;ii++)
			for (long jj=0;jj<ndet;jj++){
				iRhR[ii][jj] = 0.0;
				for (long kk=0;kk<ndet;kk++)
					iRhR[ii][jj] += iR[ii][kk]*hR[kk][jj] ;
			}


		triRhR = 0.0;
		logdetiR = 0;
		for (long ii=0;ii<ndet;ii++){
			triRhR += iRhR[ii][ii];
			logdetiR -= log(l[ii][ii]*l[ii][ii]);
		}


		f   +=  w[ib] * (triRhR - logdetiR - ndet ) ; //pb when hR non inversible

	}


	delete [] uvec;
	delete [] ivec;
	delete [] Pl;
	delete [] Pnl;
	for(int i=0; i<ndet; i++)
		delete [] l[i];
	delete [] l;

	free_dmatrix(R,0,ndet-1,0,ndet-1);
	free_dmatrix(hR,0,ndet-1,0,ndet-1);
	free_dmatrix(eR,0,ndet-1,0,ndet-1);
	free_dmatrix(iR,0,ndet-1,0,ndet-1);
	free_dmatrix(iRhR,0,ndet-1,0,ndet-1);


	return f;

}


void rescaleAP(double **A, double **P, long ndet, long ncomp, long nbins){

	double *norm2ratio;

	norm2ratio = new double[ncomp];

	fill(norm2ratio,norm2ratio+ncomp,0.0);

	for (long ii=0;ii<ncomp;ii++){
		for (long jj=0;jj<ndet;jj++)
			norm2ratio[ii] += A[jj][ii] * A[jj][ii];
		norm2ratio[ii] = 1.0/norm2ratio[ii];
	}

	for (long ii=0;ii<ncomp;ii++){
		for (long jj=0;jj<ndet;jj++)
			A[jj][ii] = A[jj][ii] * sqrt(norm2ratio[ii]) ;
		for (long ib=0;ib<nbins;ib++)
			P[ii][ib] = P[ii][ib] / norm2ratio[ii] ;
	}

	delete [] norm2ratio;

}


int write_to_disk(string outdirSpN, string fits_filename, struct param_sanePS structPS, std::vector<std::string> det,  long nbins, double *ell, double **mixmat,
		double **Rellth, double **Rellexp, double **N, double *SPref, double **P)
{


	std::ostringstream temp_stream;
	string testfile;
	string tempstr1, tempstr2;
	string nameSpfile;
	string basename;

	long ndet = (long)det.size();

	FILE *fp;
	double *data1d;

	basename = FitsBasename(fits_filename);

	temp_stream << outdirSpN << basename << structPS.ell_suffix;
	nameSpfile= temp_stream.str();
	temp_stream.str("");

	fp = fopen(nameSpfile.c_str(),"w");
	for (long ii=0;ii<nbins;ii++)
		fprintf(fp,"%g\n",ell[ii]);
	fclose(fp);

	temp_stream << "!" << outdirSpN << basename << structPS.cov_matrix_suffix;
	nameSpfile= temp_stream.str();
	temp_stream.str("");

	if(write_CovMatrix(nameSpfile, det, nbins, ell, Rellth))
		return EX_IOERR;


	temp_stream << "!" << outdirSpN << basename << "_exp" << structPS.cov_matrix_suffix;
	nameSpfile= temp_stream.str();
	temp_stream.str("");

	if(write_CovMatrix(nameSpfile, det, nbins, ell, Rellexp))
		return EX_IOERR;


#ifdef DEBUG

	temp_stream << outdirSpN << basename << "_exp.psd";
	nameSpfile= temp_stream.str();
	temp_stream.str("");
	fp = fopen(nameSpfile.c_str(),"w");

	for (long idet1=0;idet1<ndet;idet1++){
		for (long idet2=0;idet2<ndet;idet2++){

			///// write power spectrum to disk
			tempstr1 = det[idet1];
			tempstr2 = det[idet2];
			fprintf(fp,"%s%s%s\n",tempstr1.c_str(),"-",tempstr2.c_str());
			fprintf(fp,"%d\n",(int)nbins);
			for (long ii=0;ii<nbins;ii++){
				fprintf(fp,"%g\t",ell[ii]);
				fprintf(fp,"%10.15g\n",(Rellexp[idet1*(ndet)+idet2][ii]+Rellexp[idet2*(ndet)+idet1][ii])/2.0*SPref[ii]);
			}
			fprintf(fp,"%g\n",ell[nbins]);
		}
	}
	fclose(fp);

#endif

	temp_stream << outdirSpN  << basename << structPS.mix_suffix;
	testfile= temp_stream.str();
	temp_stream.str("");

	fp = fopen(testfile.c_str(),"w");
	fprintf(fp,"%i\n",structPS.ncomp);
	for (long ii=0;ii<ndet;ii++)
		for (long jj=0;jj<structPS.ncomp;jj++)
			fprintf(fp,"%10.15g \n",mixmat[ii][jj]);
	fprintf(fp,"\n");
	fclose(fp);


#ifdef DEBUG

	//**************** Write component power spectra to disk
	for (long idet1=0;idet1<ndet;idet1++){

		tempstr1 = det[idet1];
		temp_stream << outdirSpN + tempstr1 + "_uncnoise_" << basename << ".psd";
		nameSpfile= temp_stream.str();
		temp_stream.str("");

		fp = fopen(nameSpfile.c_str(),"w");
		for (long ii=0;ii<nbins;ii++){
			fprintf(fp,"%10.15g\n",N[idet1][ii]*SPref[ii]);
		}
		fprintf(fp,"\n");
		fclose(fp);
	}

	data1d = new double[ndet*nbins];
	for (long i=0; i< ndet; i++)
		for (long j=0; j<nbins; j++)
			data1d[i*nbins+j] = N[i][j]*SPref[j];

	temp_stream << "!" + outdirSpN + "Nfinal_" << basename << "_uncnoise.fits";
	testfile= temp_stream.str();
	temp_stream.str("");
	write_psd_tofits(testfile,nbins,ndet,'d',data1d);
	delete [] data1d;

	for (long jj=0;jj<structPS.ncomp;jj++){

		temp_stream << outdirSpN + "Comp_" << jj << "_uncnoise_" << basename << ".psd";
		nameSpfile= temp_stream.str();
		temp_stream.str("");

		fp = fopen(nameSpfile.c_str(),"w");
		for (long ii=0;ii<nbins;ii++){
			fprintf(fp,"%10.15g\n",P[jj][ii]*SPref[ii]);
		}
		fprintf(fp,"\n");
		fclose(fp);
	}

#endif


	data1d = new double[structPS.ncomp*nbins];
	for (long i=0; i< structPS.ncomp; i++)
		for (long j=0; j<nbins; j++)
			data1d[i*nbins+j] = P[i][j]*SPref[j];

	temp_stream << "!" + outdirSpN + "Nfinal_" << basename << "_cnoise.fits";

	// get filename
	testfile= temp_stream.str();
	temp_stream.str("");
	write_psd_tofits(testfile,nbins,structPS.ncomp,'d',data1d);



	delete [] data1d;

	return EX_OK;

}
