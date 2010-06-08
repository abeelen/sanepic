/*
 * estimPS_steps.cpp
 *
 *  Created on: 19 août 2009
 *      Author: matthieu
 */

#include "estimPS_steps.h"


#include "psdIO.h"
#include <vector>
#include "todprocess.h"
#include "map_making.h"
#include "inline_IO2.h"
#include "inputFileIO.h"
#include "dataIO.h"
#include "covMatrix_IO.h"

#include <sstream>
#include <cmath>
#include <cstring>



extern "C" {
#include <fitsio.h>
}


using namespace std;

void read_mixmat_file(string MixMatfile, string dir, double **&mixmat, long ndet, long ncomp){

	read_mixmat_txt(MixMatfile, ndet, ncomp,mixmat);

}




void common_mode_computation(struct detectors det, struct param_process proc_param, struct param_positions pos_param,
		struct common dir, double *apodwind,long ns, long ff, long NAXIS1, long NAXIS2, long long npix,
		long iframe, double *S, long long *indpix,double **mixmat, long ncomp, double **commonm2,
		double &factapod, string fits_filename)
{
	//*************************** Read data and compute components

	// samptopix, Ps, data, data_lp, fdata1, bfilter, cov, uvec,p,ivec, icov

	string field; // detector name in the loop

	double *sign;
	double tmpsign, mm;
	double sign0=1;
	double **commonm;
	int *flag;

	fftw_plan fftplan;


	double *data, *data_lp, *Ps, *bfilter;
	long long *samptopix; // sample to pixel projection matrix
	double **Cov, **iCov;
	double *p, *uvec, *ivec;


	fftw_complex *fdata1;


	int factdupl = 1;
	if(pos_param.flgdupl==1) factdupl = 2;


	data_lp = new double[ns]; // data low passed
	Ps = new double[ns]; // deprojected signal
	samptopix = new long long[ns]; // sample to pixel proj matrix
	bfilter = new double[ns/2+1]; // buttter filter values
	fdata1 = new fftw_complex[ns/2+1]; // fourier transform


	Cov = dmatrix(0,ncomp-1,0,ncomp-1); // AtN-1A
	iCov = dmatrix(0,ncomp-1,0,ncomp-1);  // inverted AtN-1A
	p = new double[ncomp];
	uvec = new double[ncomp];
	ivec = new double[ncomp];


	init2D_double(Cov,0,0,ncomp,ncomp,0.0);
	init2D_double(iCov,0,0,ncomp,ncomp,0.0);


	for(long ii=0;ii<ns/2+1;ii++)
		bfilter[ii] = 1.0;


	sign = new double[det.ndet];
	commonm = dmatrix(0,ncomp,0,ns-1); // common mode
	init2D_double(commonm,0,0,ncomp,ns,0.0);

	// test
	//	fftplan = fftw_plan_dft_r2c_1d(ns, data, fdata1, FFTW_ESTIMATE);

	// loop over detectors
	for (long idet=0;idet<det.ndet;idet++){

		field = det.boloname[idet];


		long test_ns;
		read_signal_from_fits(fits_filename, field, data, test_ns);
		if (test_ns != ns) {
			cerr << "Read signal does not correspond to frame size : Check !!" << endl;
			exit(-1);
		}

		read_flag_from_fits(fits_filename , field, flag, test_ns);
		if (test_ns != ns) {
			cerr << "Read flag does not correspond to frame size : Check !!" << endl;
			exit(-1);
		}



		//TODO: subtract the signal only if needed

		if (S != NULL){
			//Read pointing data
			read_samptopix(ns, samptopix,  dir.tmp_dir, iframe, det.boloname[idet]);


			//TODO : Check this function on what it does/should do
			deproject(S,indpix,samptopix,ns,NAXIS1,NAXIS2,npix,Ps,pos_param.flgdupl,factdupl);

#ifdef DEBUG
			FILE * fp;
			fp = fopen("testDeproject.txt","w");
			for (int i =0;i<ns;i++)
				fprintf(fp,"%lf %lf\n",data[i],Ps[i]);
			fclose(fp);
			exit(-1);
#endif


			for(long ii=0;ii<ns;ii++)
				data[ii] = data[ii] - Ps[ii];
		}

		//TODO : f_lp_pix is hard fixed to 1.0 ??????????
		MapMakPreProcessData(data,flag,ns,proc_param.napod,proc_param.poly_order,1.0,data_lp,bfilter,
				proc_param.NORMLIN,proc_param.NOFILLGAP,proc_param.remove_polynomia);

		// TODO: should apodisation be part of MapMakePreProcess ?
		for (long ii=0;ii<ns;ii++)
			data[ii] = data_lp[ii]*apodwind[ii];


		// TODO: Do we need to compute and save ffts here ???
		//       fdata are used in cross power spectrum estimation...
		//       BUT it is done differently than power spectrum estimation WHY ?
		// TODO: Check how to speed up ffts :  why destroy all plan ?
		// compute fft and save data to disk for later
		fftplan = fftw_plan_dft_r2c_1d(ns, data, fdata1, FFTW_ESTIMATE);
		fftw_execute(fftplan);
		fftw_destroy_plan(fftplan);


		write_fdata(ns, fdata1,  "fdata_", dir.tmp_dir, idet, ff, det.boloname);
		/*for(long ii=0;ii<ns/2+1;ii++){
				fdata_buffer[idet*(ns/2+1)+ii][0]=fdata1[ii][0];
				fdata_buffer[idet*(ns/2+1)+ii][1]=fdata1[ii][1];
			}*/



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

		delete [] flag; // ajout mat 26/11

	}

	//test
	//	fftw_destroy_plan(fftplan);

	for (long jj=0; jj<ncomp; jj++)
		for (long ii= 0 ;ii<ns;ii++)
			if ( isnan(commonm[jj][ii]) || isinf(commonm[jj][ii]) )
				cout << "Something went wrong in the computation of the common mode " << ncomp << endl;


	//***************************************************************************


	/////////// AtN-1A
	for (long jj=0;jj<ncomp;jj++)
		for (long kk=0;kk<ncomp;kk++)
			for (long ii=0;ii<det.ndet;ii++)
				Cov[jj][kk] += mixmat[ii][jj] * mixmat[ii][kk]/sign[ii]/sign[ii];


	// invert AtN-1A
	dcholdc(Cov,ncomp,p);
	for (long ii=0;ii<ncomp;ii++){
		for (long jj=0;jj<ncomp;jj++)
			uvec[jj] = 0.0;
		uvec[ii] = 1.0;
		dcholsl(Cov,ncomp,p,uvec,ivec);
		for (long jj=0;jj<ncomp;jj++)
			iCov[ii][jj] = ivec[jj];
	}



	printf("noise var det 0 =  %10.15g\n",sign0*sign0);


	for (long ii=0;ii<ns;ii++)
		for (long jj=0;jj<ncomp;jj++)
			for (long kk=0;kk<ncomp;kk++)
				commonm2[jj][ii] += iCov[jj][kk] * commonm[kk][ii]; //common mode * (AtN-1A)-1



	//	factapod = 0.0;
	//	for (long ii=0;ii<ns;ii++)
	//		factapod += apodwind[ii]*apodwind[ii]/ns; // factapod ?? apodization factor ?



	/// filter common mode to fcut Hz
	//for (long ii=0;ii<ncomp;ii++){
	//for (long jj=0;jj<ns;jj++)
	//commontmp[jj] = commonm2[ii][jj];
	//butterworth(commontmp,ns,fcut/fsamp*double(ns),8,commonm_f,bfilter,1,0,0);
	//for (long jj=0;jj<ns;jj++)
	//common_f[ii][jj] = commontmp[jj];                         //// - commonm_f[jj];
	//}

	// clean up
	delete [] sign;
	delete [] data ;
	delete [] data_lp ;
	delete [] Ps ;
	delete [] samptopix;
	delete [] bfilter ;
	delete [] fdata1 ;
	delete [] p ;
	delete [] uvec ;
	delete [] ivec;

	free_dmatrix(Cov,0,ncomp-1,0,ncomp-1);
	free_dmatrix(iCov,0,ncomp-1,0,ncomp-1);
	free_dmatrix(commonm,0,ncomp,0,ns-1);


	//----------------------------------- END -------------------------------//

}



void estimate_noise_PS(struct detectors det, struct param_process proc_param,struct param_positions pos_param,
		struct common dir, long &nbins,	long &nbins2, long ns, long ff, long NAXIS1,
		long NAXIS2, long long npix, double *&ell, double *S, long iframe,long long *indpix,
		double *apodwind, long ncomp, double **mixmat, double **commonm2,
		double factapod,double **Rellth, double **N, double **P, string fits_filename)
{


	string nameSpfile, field;
	string testfile;
	std::ostringstream temp_stream; // used to remove sprintf horror
	FILE *fp;

	int *flag;

	double *data, *data_lp, *Ps=NULL, *bfilter;
	double  *commontmp;
	double *Nell, *Nk;
	long long *samptopix; // sample to pixel projection matrix


	int factdupl = 1;
	if(pos_param.flgdupl==1) factdupl = 2; // map duplication factor


	data_lp = new double[ns]; // data low passed
	bfilter = new double[ns/2+1]; // buttter filter values
	commontmp = new double[ns]; //
	Nell = new double[nbins]; // binned noise PS
	Nk = new double[ns/2+1]; // noise PS

	if (S != NULL){
		samptopix = new long long[ns]; // sample to pixel proj matrix
		Ps = new double[ns];
	}

	for(long ii=0;ii<ns/2+1;ii++)
		bfilter[ii] = 1.0;


	fill(commontmp,commontmp+ns,0.0);

	//----------------------------------- ESTIMATE NOISE PS -------------------------------//


	//************************************************************************//
	// second part: -- data - common mode
	//              -- estimate noise power spectra


	/////////////////////////////////////// loop again over detectors
	for (long idet=0;idet<det.ndet;idet++){

		field = det.boloname[idet];

		long test_ns;
		read_signal_from_fits(fits_filename, field, data, test_ns);
		if (test_ns != ns) {
			cerr << "Read signal does not correspond to frame size : Check !!" << endl;
			exit(-1);
		}

		read_flag_from_fits(fits_filename , field, flag, test_ns);
		if (test_ns != ns) {
			cerr << "Read flag does not correspond to frame size : Check !!" << endl;
			exit(-1);
		}


		//TODO : This computation is already done when computing the common mode
		//       reuse the fdata if possible?
		//******************************* subtract signal

		if (S != NULL){
			//Read pointing data
			read_samptopix(ns, samptopix,  dir.tmp_dir, iframe, field);

			deproject(S,indpix,samptopix,ns,NAXIS1,NAXIS2,npix,Ps,pos_param.flgdupl,factdupl);

			for(long ii=0;ii<ns;ii++)
				data[ii] = data[ii] - Ps[ii];

		}


		MapMakPreProcessData(data,flag,ns,proc_param.napod,proc_param.poly_order,1.0,data_lp,bfilter,
				proc_param.NORMLIN,proc_param.NOFILLGAP,proc_param.remove_polynomia);

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
			Rellth[idet*det.ndet+idet][ii] += Nell[ii]/factapod; // uncorrelated part added in covariance matrix ??
			N[idet][ii] = Nell[ii]/factapod; // uncorrelated part
		}

		delete [] flag;
	}




	////*********************** Component power spectra

	for (long ii=0;ii<ncomp;ii++){
		for (long jj=0;jj<ns;jj++)
			commontmp[jj]=commonm2[ii][jj];
		noisepectrum_estim(commontmp,ns,ell,(int)nbins,proc_param.fsamp,NULL,Nell,Nk);
		for (long jj=0;jj<nbins;jj++)
			P[ii][jj] = Nell[jj]/factapod;

		// subtract a factor to correct from noise
		//for (jj=0;jj<nbins;jj++)
		//  P[ii][jj] -= iCov[ii][ii]*sign0*sign0;
		//for (jj=0;jj<nbins;jj++)
		//  if (P[ii][jj] < 0)
		//	P[ii][jj] = 0.0;

		//printf("iCov[%ld] = %10.15g\n",ii,iCov[ii][ii]*sign0*sign0);

		temp_stream << dir.output_dir + "Nellc_" << ii << "_" << ff << ".bi";
		// Get the string
		testfile= temp_stream.str();
		// Clear ostringstream buffer
		temp_stream.str("");

		if((fp = fopen(testfile.c_str(),"w"))){
			fwrite(Nell,sizeof(double), nbins, fp);
			fclose(fp);
		}else{
			cerr << "Error. Can't open " << testfile << ".Exiting.\n";
			exit(0);
		}

	}

	for (long ii=0;ii<det.ndet;ii++)
		for (long kk=0;kk<det.ndet;kk++)
			for (long ll=0;ll<ncomp;ll++)
				for (long jj=0;jj<nbins;jj++)
					Rellth[ii*det.ndet+kk][jj] += mixmat[ii][ll] * mixmat[kk][ll] * P[ll][jj]; // add correlated part to covariance matrix

	// clean up
	if (S != NULL){
		delete [] Ps ;
		delete [] samptopix;
	}

	delete [] data_lp ;
	delete [] data ;
	delete [] bfilter ;
	delete [] commontmp;
	delete [] Nell;
	delete [] Nk;

	//----------------------------------- END -------------------------------//
}


void estimate_CovMat_of_Rexp(struct common dir, struct detectors det, long nbins, long ns, long ff, double *ell, long ncomp, double **mixmat,double fsamp,
		double factapod,double **Rellexp, double **N, double **P, double *SPref)
{

	std::ostringstream temp_stream; // used to remove sprintf horror

	double *data1d; // buffer used to write down 1d array
	data1d = new double[det.ndet*nbins];
	string testfile;
	fftw_complex *fdata1, *fdata2;
	double * Nell, *Nk;

	Nell = new double[nbins]; // binned noise PS
	Nk = new double[ns/2+1]; // noise PS

	fdata1 = new fftw_complex[ns/2+1]; // fourier transform
	fdata2 = new fftw_complex[ns/2+1];

	FILE *fp;

	/////////////////////////////////////// loop again over detectors
	for (long idet1=0;idet1<det.ndet;idet1++){

		// read data from disk
		read_fdata(ns, fdata1, "fdata_",  dir.tmp_dir, idet1, ff, det.boloname);

		/*for(long ii=0;ii<ns/2+1;ii++){
			fdata1[ii][0]=fdata_buffer[(ns/2+1)*idet1+ii][0];
			fdata1[ii][1]=fdata_buffer[(ns/2+1)*idet1+ii][1];
		}*/


		for (long idet2=0;idet2<det.ndet;idet2++) {

			// read data from disk
			read_fdata(ns, fdata2, "fdata_",  dir.tmp_dir, idet2, ff, det.boloname);

			/*for(long ii=0;ii<ns/2+1;ii++){
						fdata2[ii][0]=fdata_buffer[(ns/2+1)*idet2+ii][0];
						fdata2[ii][1]=fdata_buffer[(ns/2+1)*idet2+ii][1];
	}*/

			noisecrosspectrum_estim(fdata1,fdata2,ns,ell,(int)nbins,fsamp,NULL,Nell,Nk);


			for (long ii=0;ii<nbins;ii++)
				Rellexp[idet1*det.ndet+idet2][ii] += Nell[ii]/factapod; // noise cross PS ?

		}

		cout << "Computing Rellexp :" << idet1*100./det.ndet << "%\r" << flush ;

	}

	////// normalize to the first detector power spectrum in order to avoid numerical problems
	for (long ii=0;ii<nbins;ii++)
		SPref[ii] = Rellexp[0][ii]; // first detector PS

	for (long ii=0;ii<nbins;ii++)
		if ( isnan(SPref[ii]) || isinf(SPref[ii]))
			cout << " Problem in the first detector power spectrum\n";

	for (long idet1=0;idet1<det.ndet;idet1++)
		for (long idet2=0;idet2<det.ndet;idet2++)
			for (long ii=0;ii<nbins;ii++)
				Rellexp[idet1*det.ndet+idet2][ii] = Rellexp[idet1*det.ndet+idet2][ii]/SPref[ii]; // normalize to first detector
	for (long jj=0;jj<ncomp;jj++)
		for (long ii=0;ii<nbins;ii++)
			P[jj][ii] = P[jj][ii]/SPref[ii]; // normalize common mode part
	for (long jj=0;jj<det.ndet;jj++)
		for (long ii=0;ii<nbins;ii++)
			N[jj][ii] = N[jj][ii]/SPref[ii]; //normalize uncorrelated part


	// write Rellexp to disk and also first guess of parameters
	temp_stream << dir.output_dir + "Rellexp_" << ff << ".txt";

	// get filename
	testfile= temp_stream.str();
	temp_stream.str("");

	fp = fopen(testfile.c_str(),"w");
	for (long jj=0;jj<nbins;jj++)
		for (long ii=0;ii<det.ndet;ii++)
			for (long kk=0;kk<det.ndet;kk++)
				fprintf(fp,"%10.15g \t",Rellexp[ii*det.ndet+kk][jj]); // cross power spectrum
	fprintf(fp,"\n");
	fclose(fp);

	temp_stream << dir.output_dir + "Ninit_" << ff << ".txt";

	// get filename
	testfile= temp_stream.str();
	temp_stream.str("");

	fp = fopen(testfile.c_str(),"w");
	for (long ii=0;ii<det.ndet;ii++){
		for (long jj=0;jj<nbins;jj++)
			fprintf(fp,"%10.15g \t",N[ii][jj]); // uncorralated part of the noise
		fprintf(fp,"\n");
	}
	fclose(fp);


	for (long i=0; i< det.ndet; i++)
		for (long j=0; j<nbins; j++)
			data1d[i*nbins+j] = N[i][j];

	temp_stream << "!" + dir.output_dir + "Ninit_" << ff << ".fits";

	// get filename
	testfile= temp_stream.str();
	temp_stream.str("");

	write_psd_tofits(testfile.c_str(),det.ndet,nbins,'d', data1d); //resized uncorralated part

	temp_stream << dir.output_dir + "Pinit_" << ff << ".txt";

	// get filename
	testfile= temp_stream.str();
	temp_stream.str("");

	fp = fopen(testfile.c_str(),"w");
	for (long ii=0;ii<ncomp;ii++)
		for (long jj=0;jj<nbins;jj++)
			fprintf(fp,"%10.15g \t",P[ii][jj]); // common mode part of the noise
	fprintf(fp,"\n");
	fclose(fp);


	temp_stream << dir.output_dir + "Ainit_" << ff << ".txt";

	// get filename
	testfile= temp_stream.str();
	temp_stream.str("");

	fp = fopen(testfile.c_str(),"w");
	for (long ii=0;ii<det.ndet;ii++)
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
	//----------------------------------- END -------------------------------//
}


void expectation_maximization_algorithm(double fcut, long nbins, long ndet, long ncomp,long ns, double fsamp,
		long ff, string outdirSpN,	double **Rellexp, double **Rellth, double **mixmat,double **P,double **N,
		double *SPref, double *ell)
{


	//***** Fourth part
	//*********************** fit component and noise power spectra, and mixing matrix *************//
	//********* Using Expectation/Maximization algorithm
	long nbiter = 500;
	long ib=0; // used as an iterator
	//long iter;

	double tottest=0.0;

	double f;
	//	FILE *fp;

	double *iN, *Pr, *w;
	double **Rxs, **Rxsq, **RnRxsb, **Rxx, **Rxxq, **Rss, **Rssq, **RnRssb;
	double **Pr2, **AiNA, **Mattmp, **ImDR, **ACq, **Cq, **Wq;
	long nbins2=0;

	//ib = 0;
	while ((ell[ib] < fcut) && (ib < nbins)){
		nbins2 = ib+1;
		ib++;
	}



	double **Cov, **iCov;
	double *p, *uvec, *ivec;



	Cov = dmatrix(0,ncomp-1,0,ncomp-1); // AtN-1A
	iCov = dmatrix(0,ncomp-1,0,ncomp-1);  // inverted AtN-1A
	p = new double[ncomp];
	uvec = new double[ncomp];
	ivec = new double[ncomp];


	init2D_double(Cov,0,0,ncomp,ncomp,0.0);
	init2D_double(iCov,0,0,ncomp,ncomp,0.0);

	printf("\nnbins2 = %ld                \n",nbins2);



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



	//	std::ostringstream temp_stream;
	//
	//	//sprintf(testfile,"%s%s%d%s%s",outdirSpN.c_str(),"weights_",(int)ff,termin.c_str(),".txt");
	//	temp_stream << outdirSpN + "weights_" << ff << ".txt";
	//
	//	string testfile;
	//// récupérer une chaîne de caractères
	//	testfile= temp_stream.str();
	//	temp_stream.str("");
	//
	//	fp = fopen(testfile.c_str(),"w");
	//	for (long jj=0;jj<nbins2;jj++){
	//		fprintf(fp,"%10.15g \t",w[jj]);
	//	}
	//	fprintf(fp,"\n");
	//	fclose(fp);


	f = fdsf(Rellexp,w,mixmat,P,N,ndet,ncomp,nbins2) ;
	printf("Pre em:   obj: %10.15g\n", f) ;

	cout << "nbins2 " << nbins2 << endl;
	cout << "ndet   " << ndet << endl;
	cout << "ncomp  " << ncomp << endl << endl;

	for (long iter=1;iter<=nbiter;iter++){

		//init1D_double(iN,0,ndet,0.0);
		//init1D_double(Pr,0,ncomp,0.0);
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
				for (long jj=0;jj<ncomp;jj++)
					Cov[ii][jj] = Pr2[ii][jj] * AiNA[ii][jj];
				Cov[ii][ii] += 1.0;
			}

			// invert matrix
			dcholdc(Cov,ncomp,p);

			for (long ii=0;ii<ncomp;ii++){
				for (long jj=0;jj<ncomp;jj++)
					uvec[jj] = 0.0;
				uvec[ii] = 1.0;
				dcholsl(Cov,ncomp,p,uvec,ivec);
				for (long jj=0;jj<ncomp;jj++)
					iCov[ii][jj] = ivec[jj];
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
				for (long jj=0;jj<ncomp;jj++)
					Cov[ii][jj] = RnRssb[ii][jj+idet*ncomp];
			}

			// solving the linear system
			dcholdc(Cov,ncomp,p);
			dcholsl(Cov,ncomp,p,uvec,ivec);
			for (long ii=0;ii<ncomp;ii++)
				mixmat[idet][ii] = ivec[ii];
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
				for (long jj=0;jj<ncomp;jj++)
					Cov[ii][jj] = Pr2[ii][jj] * AiNA[ii][jj];
				Cov[ii][ii] += 1.0;
			}



			// invert matrix
			dcholdc(Cov,ncomp,p);
			for (long ii=0;ii<ncomp;ii++){
				for (long jj=0;jj<ncomp;jj++)
					uvec[jj] = 0.0;
				uvec[ii] = 1.0;
				dcholsl(Cov,ncomp,p,uvec,ivec);
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




		//printf("A[1][2] =  %10.15g\n", mixmat[1][2]) ;
		//printf("N[2][3] =  %10.15g\n", N[2][3]) ;Rellth
		//printf("P[2][3] =  %10.15g\n", P[2][3]) ;


		// Fixing the indeterminacies.  Is it useful here?
		//    rescaleAP(mixmat, P, ndet, ncomp, nbins2) ;



		///// here is the problem

		f = fdsf(Rellexp,w,mixmat,P,N,ndet,ncomp,nbins2) ;
		cout << "em->iter: " << iter*100./nbiter << " %\r" << flush;
		if (isnan(f) || isinf(f)) {
			cout << "Nan........." << endl;
			exit(1);
		}


	}


	// Fixing the indeterminacies.  Is it useful here?
	rescaleAP(mixmat, P, ndet, ncomp, nbins2) ;




	//****************************** Compute covariance matrix from the fitted model


	printf("EM step completed\n");


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

	//----------------------------------- END -------------------------------//


	//cleaning up

	delete [] p ;
	delete [] uvec ;
	delete [] ivec;
	delete [] w;

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
}


double fdsf(double **Rellexp, double *w, double **A, double **P, double **N, long ndet, long ncomp, long nbins){

	double f;

	//long ib, ii, jj, kk;
	double triRhR, logdetiR;

	double *p, *uvec, *ivec;
	double *Pl, *Pnl;
	double **R, **hR, **eR, **iR, **iRhR;

	p = new double[ndet];
	uvec = new double[ndet];
	ivec = new double[ndet];

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
	//init1D_double(Pl,0,ncomp,0.0);
	//init1D_double(Pnl,0,ndet,0.0);
	fill(Pl,Pl+ncomp,0.0);
	fill(Pnl,Pnl+ndet,0.0);
	init2D_double(iRhR,0,0,ndet,ndet,0.0);


	// init
	f   = 0. ;

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


		//printf("ib=%d\n",ib);
		/// inverting Rth
		for (long ii=0;ii<ndet;ii++)
			for (long jj=0;jj<ndet;jj++)
				eR[ii][jj] = R[ii][jj];
		dcholdc(eR,ndet,p);
		for (long ii=0;ii<ndet;ii++){
			for (long jj=0;jj<ndet;jj++)
				uvec[jj] = 0.0;
			uvec[ii] = 1.0;
			dcholsl(eR,ndet,p,uvec,ivec);
			for (long jj=0;jj<ndet;jj++)
				iR[ii][jj] = ivec[jj];
		}
		// printf("ib=%d\n",ib);


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
			//      cout << ii << " " << iRhR[ii][ii] << " , " << p[ii] << " : " << triRhR << " , " << detiR << endl;
			triRhR += iRhR[ii][ii];
			logdetiR -= log(p[ii]*p[ii]);
		}


		//f   +=  w[ib] * ( triRhR - log(det(iRhR)) - ndet ) ;
		f   +=  w[ib] * (triRhR - logdetiR - ndet ) ; //pb when hR non inversible

	}


	delete [] p;
	delete [] uvec;
	delete [] ivec;
	delete [] Pl;
	delete [] Pnl;

	free_dmatrix(R,0,ndet-1,0,ndet-1);
	free_dmatrix(hR,0,ndet-1,0,ndet-1);
	free_dmatrix(eR,0,ndet-1,0,ndet-1);
	free_dmatrix(iR,0,ndet-1,0,ndet-1);
	free_dmatrix(iRhR,0,ndet-1,0,ndet-1);




	return f;

}


void rescaleAP(double **A, double **P, long ndet, long ncomp, long nbins){

	//long ii, jj, ib;
	double *norm2ratio;

	norm2ratio = new double[ncomp];

	//init1D_double(norm2ratio,0,ncomp,0.0);
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


void write_to_disk(string outdirSpN,struct samples samples_struct, long ff, struct detectors det,	long nbins, double *ell, double **mixmat,
		double **Rellth, double **Rellexp, long ncomp,double **N, double *SPref, double **P)
{


	std::ostringstream temp_stream;
	string testfile;
	string tempstr1, tempstr2;
	string nameSpfile;
	string base_n;

	FILE *fp;
	double *data1d;

	base_n = FitsBasename(samples_struct.fitsvect[ff]);

	temp_stream << "!" + outdirSpN + "BoloPS_" << base_n << "_psd.fits";

	// récupérer une chaîne de caractères
	nameSpfile= temp_stream.str();
	temp_stream.str("");

	write_CovMatrix(nameSpfile, det.boloname, nbins, ell, Rellth);


	temp_stream << outdirSpN + "Ell_" << base_n << ".psd";

	// get filename
	nameSpfile= temp_stream.str();
	temp_stream.str("");

	fp = fopen(nameSpfile.c_str(),"w");
	for (long ii=0;ii<nbins;ii++){
		fprintf(fp,"%g\n",ell[ii]);
	}
	fclose(fp);

	temp_stream << outdirSpN + "BoloPS" << base_n << "_exp.psd";

	// get filename
	nameSpfile= temp_stream.str();
	temp_stream.str("");

	fp = fopen(nameSpfile.c_str(),"w");

	for (long idet1=0;idet1<det.ndet;idet1++){
		for (long idet2=0;idet2<det.ndet;idet2++){

			///// write power spectrum to disk
			tempstr1 = det.boloname[idet1];
			tempstr2 = det.boloname[idet2];
			fprintf(fp,"%s%s%s\n",tempstr1.c_str(),"-",tempstr2.c_str());
			fprintf(fp,"%d\n",(int)nbins);
			for (long ii=0;ii<nbins;ii++){
				fprintf(fp,"%g\t",ell[ii]);
				fprintf(fp,"%10.15g\n",(Rellexp[idet1*(det.ndet)+idet2][ii]+Rellexp[idet2*(det.ndet)+idet1][ii])/2.0*SPref[ii]);
			}
			fprintf(fp,"%g\n",ell[nbins]);
		}
	}
	fclose(fp);



	temp_stream << outdirSpN + "Afinal_" << base_n << ".txt";

	// get filename
	testfile= temp_stream.str();
	temp_stream.str("");

	fp = fopen(testfile.c_str(),"w");
	for (long ii=0;ii<det.ndet;ii++)
		for (long jj=0;jj<ncomp;jj++)
			fprintf(fp,"%10.15g \n",mixmat[ii][jj]);
	fprintf(fp,"\n");
	fclose(fp);



	//TODO : One should define a fits format for that
	//**************** Write component power spectra to disk
	for (long idet1=0;idet1<det.ndet;idet1++){

		tempstr1 = det.boloname[idet1];
		temp_stream << outdirSpN + tempstr1 + "_uncnoise" << base_n << ".psd";

		// récupérer une chaîne de caractères
		nameSpfile= temp_stream.str();
		temp_stream.str("");

		fp = fopen(nameSpfile.c_str(),"w");
		for (long ii=0;ii<nbins;ii++){
			fprintf(fp,"%10.15g\n",N[idet1][ii]*SPref[ii]);
		}
		fprintf(fp,"\n");
		fclose(fp);
	}

	data1d = new double[det.ndet*nbins];
	for (long i=0; i< det.ndet; i++)
		for (long j=0; j<nbins; j++)
			data1d[i*nbins+j] = N[i][j]*SPref[j];

	temp_stream << "!" + outdirSpN + "Nfinal_" << base_n << "_uncnoise.fits";

	// get filename
	testfile= temp_stream.str();
	temp_stream.str("");
	write_psd_tofits(testfile,nbins,det.ndet,'d',data1d);
	delete [] data1d;



	//TODO: One should define a fits format for that
	//TODO: Does not appear in the output????
	for (long jj=0;jj<ncomp;jj++){

		temp_stream << outdirSpN + "Comp_" << jj << "_uncnoise" << base_n << ".psd";

		// get filename
		nameSpfile= temp_stream.str();
		temp_stream.str("");

		fp = fopen(nameSpfile.c_str(),"w");
		for (long ii=0;ii<nbins;ii++){
			fprintf(fp,"%10.15g\n",P[jj][ii]*SPref[ii]);
		}
		fprintf(fp,"\n");
		fclose(fp);
	}

	data1d = new double[ncomp*nbins];
	for (long i=0; i< ncomp; i++)
		for (long j=0; j<nbins; j++)
			data1d[i*nbins+j] = P[i][j]*SPref[j];

	temp_stream << "!" + outdirSpN + "Nfinal_" << base_n << "_cnoise.fits";

	// get filename
	testfile= temp_stream.str();
	temp_stream.str("");
	write_psd_tofits(testfile,nbins,ncomp,'d',data1d);


	delete [] data1d;

}
