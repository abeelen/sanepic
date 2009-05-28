/*
 * Corr_preprocess.cpp
 *
 *  Created on: 28 mai 2009
 *      Author: matthieu
 */

#include "Corr_preprocess.h"

using namespace std;


void write_tfAS(double *S, long *indpix, int nn, long npix, bool flgdupl, int factdupl, string dir, string termin, long ff, long ns, long marge, long ndet, long iframe){


	long idet1;
	long ndata = ns+2*marge;

	FILE *fp;
	char testfile[100];

	double *Ps;
	long *samptopix;

	fftw_plan fftplan;
	fftw_complex *fdata;

	samptopix = new long[ns];
	Ps = new double[ndata];
	fdata = new fftw_complex[ndata/2+1];


	//for (idet1=rank*ndet/size;idet1<(rank+1)*ndet/size;idet1++){
	for (idet1=0;idet1<ndet;idet1++){

		//Read pointing data
		sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"samptopix_",iframe,"_",idet1,"_",termin.c_str(),".bi");
		fp = fopen(testfile,"r");
		fread(samptopix,sizeof(long),ns,fp);
		fclose(fp);

		deproject(S,indpix,samptopix,ndata,marge,nn,npix,Ps,flgdupl,factdupl);

		//Fourier transform of the data
		fftplan = fftw_plan_dft_r2c_1d(ndata, Ps, fdata, FFTW_ESTIMATE);
		fftw_execute(fftplan);
		fftw_destroy_plan(fftplan);

		sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"fPs_",iframe,"_",idet1,"_",termin.c_str(),".bi");
		fp = fopen(testfile,"w");
		fwrite(fdata,sizeof(double), (ndata/2+1)*2, fp);
		fclose(fp);

	}

	delete[] samptopix;
	delete[] Ps;
	delete[] fdata;

}






void write_ftrProcesdata(double *S, long *indpix, long *indpsrc, int nn, long npix,
		long npixsrc, long ntotscan, long addnpix, bool flgdupl, int factdupl,
		int fillg, string dir, string termin, double errarcsec, string dirfile,
		string scerr_field, string flpoint_field, string *bolonames,
		string bextension, string fextension, string cextension,
		int shift_data_to_point, double f_lppix, long ff, long ns,
		long marge, long napod, long ndet, bool NORMLIN, bool NOFILLGAP,
		long iframe){


	long ii, idet1;
	long ndata = ns+2*marge;

	double *scerr, *data, *calp, *bfilter, *data_lp, *Ps;
	unsigned char *flpoint, *flag, *rejectsamp;
	long *samptopix;

	fftw_plan fftplan;
	fftw_complex *fdata;

	string field1;

	char testfile[100];

	FILE *fp;

	scerr = new double[ndata];
	data =  new double[ndata];
	data_lp = new double[ndata];
	calp =  new double[ndata];
	flag =  new unsigned char[ndata];
	flpoint = new unsigned char[ndata];
	rejectsamp = new unsigned char[ndata];

	samptopix = new long[ns];
	Ps = new double[ndata];
	bfilter = new double[ndata/2+1];
	fdata = new fftw_complex[ndata/2+1];



	for (idet1=0;idet1<ndet;idet1++){

		field1 = bolonames[idet1];
		//     cout << "Inside write_ftrProcessdata " << endl;
		//     cout << " field 1 : " << field1 << endl;

		if (S != NULL){
			read_data_std(dirfile, ff, 0, ns, scerr, scerr_field, 'd');
			read_data_std(dirfile, ff, 0, ns, flpoint, flpoint_field, 'c');
		}

		read_data_std(dirfile, ff, shift_data_to_point, ns, data, field1+bextension, 'd');

		if (fextension != "NOFLAG"){
			read_data_std(dirfile, ff, shift_data_to_point, ns, flag, field1+fextension,  'c');
		} else {
			//      printf("NOFLAG\n");
			for (ii=0;ii<ns;ii++)
				flag[ii] = 0;
		}

		if (cextension != "NOCALP"){
			read_data_std(dirfile, ff, 0, ns/20, calp, field1+cextension, 'd');
		} else {
			//      printf("NOCALP\n");
			for (ii=0;ii<ns/20;ii++)
				calp[ii] = 1.0;
		}



		if (S != NULL){
			//// Read pointing
			sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"samptopix_",iframe,"_",idet1,"_",termin.c_str(),".bi");
			fp = fopen(testfile,"r");
			fread(samptopix,sizeof(long),ns,fp);
			fclose(fp);

			if (addnpix){
				deproject(S,indpix,samptopix,ns+2*marge,marge,nn,npix,Ps,fillg,factdupl,ntotscan,indpsrc,npixsrc);
			} else {
				deproject(S,indpix,samptopix,ns+2*marge,marge,nn,npix,Ps,fillg,factdupl);
			}

			for (ii=0;ii<ns;ii++) rejectsamp[ii] = 0;
			for (ii=0;ii<ns;ii++)
				if ((flag[ii] & 1) != 0 || (scerr[ii] > errarcsec) || (flpoint[ii] & 1) != 0)
					rejectsamp[ii] = 1;
		}


		if (S != NULL){
			//********************  pre-processing of data ********************//
			MapMakPreProcessData(data,flag,calp,ns,marge,napod,4,f_lppix,data_lp,bfilter,
					NORMLIN,NOFILLGAP,Ps);
		}
		else {
			MapMakPreProcessData(data,flag,calp,ns,marge,napod,4,f_lppix,data_lp,bfilter,
					NORMLIN,NOFILLGAP);
		}

		//Fourier transform of the data
		fftplan = fftw_plan_dft_r2c_1d(ndata, data_lp, fdata, FFTW_ESTIMATE);
		fftw_execute(fftplan);
		fftw_destroy_plan(fftplan);


		sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"fdata_",iframe,"_",idet1,"_",termin.c_str(),".bi");
		fp = fopen(testfile,"w");
		fwrite(fdata,sizeof(double), (ndata/2+1)*2, fp);
		fclose(fp);

	}


	delete[] scerr;
	delete[] data;
	delete[] data_lp;
	delete[] calp;
	delete[] flag;
	delete[] flpoint;
	delete[] samptopix;
	delete[] Ps;
	delete[] bfilter;
	delete[] fdata;
	delete[] rejectsamp;


}




void do_PtNd(double *PNd, string *extentnoiseSp_all, string noiseSppreffile,
		string dir, string prefixe, string termin, string *bolonames,
		double f_lppix, double fsamp, long ff, long ns, long marge, long ndet, int size,
		int rank, long *indpix, long nn, long npix, long iframe, double *Mp, long *hits){


	long ii, jj, idet1, idet2, nbins;
	double dnbins;
	long ndata = ns+2*marge;
	string field1, field2;
	string extentnoiseSp;

	char nameSpfile[100];
	char testfile[100];

	long *samptopix;
	double *ell, *SpN, *bfilter, *bfilter_, *Nk, *Nd;

	fftw_plan fftplan;
	fftw_complex *fdata, *Ndf;


	samptopix = new long[ns];
	Nd = new double[ndata];
	bfilter = new double[ndata/2+1];
	bfilter_ = new double[ndata/2+1];
	Nk = new double[ndata/2+1];
	fdata = new fftw_complex[ndata/2+1];
	Ndf = new fftw_complex[ndata/2+1];

	double **SpN_all;

	FILE *fp;




	for (idet1=rank*ndet/size;idet1<(rank+1)*ndet/size;idet1++){
		field1 = bolonames[idet1];

		//Read pointing data
		sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"samptopix_",iframe,"_",idet1,"_",termin.c_str(),".bi");
		fp = fopen(testfile,"r");
		fread(samptopix,sizeof(long),ns,fp);
		fclose(fp);


		//**************************************** Noise power spectrum
		extentnoiseSp = extentnoiseSp_all[iframe];
		sprintf(nameSpfile,"%s%s%s%s",noiseSppreffile.c_str(),field1.c_str(),"-all",extentnoiseSp.c_str());
		if ((fp = fopen(nameSpfile,"r")) == NULL){
			cerr << "ERROR: Can't find noise power spectra file" << nameSpfile << " , check -k or -K in command line. Exiting. \n";
			exit(1);
		}
		fread(&dnbins,sizeof(double), 1, fp);
		nbins = (long)dnbins;
		SpN_all = dmatrix(0,ndet-1,0,nbins-1);
		ell = new double[nbins+1];
		SpN = new double[nbins];
		fread(ell,sizeof(double), nbins+1, fp);
		fread(*SpN_all,sizeof(double), nbins*ndet, fp);
		fclose(fp);
		//*****************************************

		for (ii=0;ii<ndata/2+1;ii++)
			bfilter[ii] = pow(double(ii)/f_lppix, 16) /(1.0+pow(double(ii)/f_lppix, 16));
		for (ii=0;ii<ndata/2+1;ii++)
			bfilter_[ii] = 1.0/(bfilter[ii]+0.000001);


		//Init N-1d
		for (ii=0;ii<ndata/2+1;ii++){
			Ndf[ii][0] = 0;
			Ndf[ii][1] = 0;
		}



		for (idet2=0;idet2<ndet;idet2++){
			field2 = bolonames[idet2];

			//read Fourier transform of the data
			sprintf(testfile,"%s%s%s%ld%s%ld%s%s%s",dir.c_str(),prefixe.c_str(),"_",iframe,"_",idet2,"_",termin.c_str(),".bi");
			fp = fopen(testfile,"r");
			fread(fdata,sizeof(double), (ndata/2+1)*2, fp);
			fclose(fp);


			//****************** Cross power spectrum of the noise  ***************//
			for (ii=0;ii<nbins;ii++)
				SpN[ii] = SpN_all[idet2][ii];


			// interpolate logarithmically the noise power spectrum
			InvbinnedSpectrum2log_interpol(ell,SpN,bfilter_,nbins,ndata,fsamp,Nk);


			for (jj=0;jj<ndata/2+1;jj++)
				if (isnan(Nk[jj])) {
					printf("Ca ne va pas fr %ld, det1 %ld, det2 %ld\n",iframe, idet1, idet2);
					exit(1);
				}

			//********************************* compute N^-1 d  ***********************//
			for (ii=0;ii<ndata/2+1;ii++){
				Ndf[ii][0] += fdata[ii][0]*Nk[ii];
				Ndf[ii][1] += fdata[ii][1]*Nk[ii];
			}



			//Compute weight map for preconditioner
			if ((Mp != NULL) && (idet2 == idet1))
				compute_diagPtNPCorr(Nk,samptopix,ndata,marge,nn,indpix,npix,f_lppix,Mp);


		}// end of idet2 loop


		fftplan = fftw_plan_dft_c2r_1d(ndata, Ndf, Nd, FFTW_ESTIMATE);
		fftw_execute(fftplan);
		fftw_destroy_plan(fftplan);



		for (ii=-marge;ii<ndata-marge;ii++){
			if ((ii < 0) || (ii >= ndata-2*marge)){
				PNd[npix-2] += Nd[ii+marge];
			} else {
				PNd[indpix[samptopix[ii]]] += Nd[ii+marge];
			}
		}

		//compute hit counts
		if (hits != NULL){
			for (ii=0;ii<ndata-2*marge;ii++){
				hits[indpix[samptopix[ii]]] += 1;
			}
		}


		delete[] ell;
		delete[] SpN;
		free_dmatrix(SpN_all,0,ndet-1,0,nbins-1);

	}// end of idet1 loop



	delete[] samptopix;
	delete[] Nd;
	delete[] bfilter;
	delete[] bfilter_;
	delete[] Nk;
	delete[] fdata;
	delete[] Ndf;


}

