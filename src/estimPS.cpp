/*
 * estimPS.cpp
 *
 *  Created on: 28 mai 2009
 *      Author: matthieu
 */

#include "estimPS.h"

using namespace std;


void EstimPowerSpectra(double fsamp, long ns, long ff, long ndet, int nn, long npix, long napod,
		long iframe, bool flgdupl, int factdupl, long *indpix,
		double *S, string MixMatfile, string *bolonames, string dirfile, string bextension,
		string fextension, string cextension, int shift_data_to_point, string dir,
		string termin, bool NORMLIN, bool NOFILLGAP, string noiseSppreffile,
		string extentnoiseSp, string outdirSpN)
{



	double fcut = 12;
	long ncomp = 1;

	long ii, jj, kk, ll, idet, idet1, idet2, ib;
	long ncomp2 = 0;
	int nbins = 500;
	int nbins2;
	double tmpsign, mm, sign0, factapod;

	double dnbins, dummy1;

	double *data, *data_lp, *Ps, *calp, *apodwind, *commontmp, *commonm_f, *bfilter, *SPref;
	long *samptopix;
	unsigned char *flag;
	double **commonm, **commonm2, **common_f, **vect;
	double **P, **N, **Rellth, **Rellexp;
	double **aa, **Cov, **iCov, **iCov2, **SpN_all;
	double *Nk;
	double *ell, *SpN, *Nell;
	double *sign;
	double *p, *uvec, *ivec;

	fftw_plan fftplan;
	fftw_complex *fdata1, *fdata2;

	FILE *fp;

	double *data1d;


	char testfile[100];
	char nameSpfile[100];

	string field;
	string tempstr1;
	string tempstr2;


	printf("Inside EstimPowerSpectra just at the beginning\n");


	data = new double[ns];
	data_lp = new double[ns];
	Ps = new double[ns];
	commontmp = new double[ns];
	commonm_f = new double[ns];
	calp = new double[ns];
	flag = new unsigned char[ns];
	samptopix = new long[ns];
	Nk = new double[ns/2+1];
	bfilter = new double[ns/2+1];
	fdata1 = new fftw_complex[ns/2+1];
	fdata2 = new fftw_complex[ns/2+1];


	Nell = new double[nbins];
	SPref = new double[nbins];
	P = dmatrix(0,ncomp-1,0,nbins-1);
	N = dmatrix(0,ndet-1,0,nbins-1);
	Rellexp = dmatrix(0,ndet*ndet-1,0,nbins-1);
	Rellth = dmatrix(0,ndet*ndet-1,0,nbins-1);
	aa = dmatrix(0,ndet-1,0,20);

	sign = new double[ndet];

	Cov = dmatrix(0,ncomp-1,0,ncomp-1);
	iCov = dmatrix(0,ncomp-1,0,ncomp-1);
	iCov2 = dmatrix(0,ncomp-1,0,ncomp-1);
	p = new double[ncomp];
	uvec = new double[ncomp];
	ivec = new double[ncomp];

	vect = dmatrix(0,ncomp-1,0,ndet-1);

	commonm = dmatrix(0,ncomp,0,ns-1);
	commonm2 = dmatrix(0,ncomp,0,ns-1);
	common_f = dmatrix(0,ncomp,0,ns-1);

	init2D_double(commonm,0,0,ncomp,ns,0.0);
	init2D_double(commonm2,0,0,ncomp,ns,0.0);
	init2D_double(common_f,0,0,ncomp,ns,0.0);
	init1D_double(commontmp,0,ns,0.0);
	init1D_double(commonm_f,0,ns,0.0);
	init2D_double(Cov,0,0,ncomp,ncomp,0.0);
	init2D_double(iCov,0,0,ncomp,ncomp,0.0);
	init2D_double(iCov2,0,0,ncomp,ncomp,0.0);

	for (ii=0;ii<ns/2+1;ii++)
		bfilter[ii] = 1.0;




	//////////////////////////////////////
	// read mixing parameters
	if ((fp = fopen(MixMatfile.c_str(),"r")) == NULL){
		cerr << "ERROR: Can't find Mixing Matrix file. Exiting. \n";
		exit(1);
	}
	fscanf(fp,"%ld",&ncomp2); // modified d => ld to avoid warning, mat-27/05

	for (ii=0;ii<ndet;ii++){
		for (jj=0;jj<ncomp2;jj++){
			fscanf(fp,"%lf",&dummy1);
			aa[ii][jj] = dummy1;
		}
	}
	fclose(fp);

	if (ncomp2 < ncomp) ncomp = ncomp2;




	//**************************************************************************************
	//**************************** Read data and compute components

	//apodization
	apodwind = apodwindow(ns,int(ns*0.04));


	// loop over detectors
	for (idet=0;idet<ndet;idet++){

		field = bolonames[idet];


		read_data_std(dirfile, ff, shift_data_to_point, ns, data, field+bextension, 'd');

		if (fextension != "NOFLAG"){
			read_data_std(dirfile, ff, shift_data_to_point, ns, flag, field+fextension,  'c');
		} else {
			//      printf("NOFLAG\n");
			for (ii=0;ii<ns;ii++)
				flag[ii] = 0;
		}

		if (cextension != "NOCALP"){
			read_data_std(dirfile, ff, 0, ns/20, calp, field+cextension, 'd');
		} else {
			//      printf("NOCALP\n");
			for (ii=0;ii<ns/20;ii++)
				calp[ii] = 1.0;
		}


		//******************************* subtract signal

		//Read pointing data
		sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"samptopix_",iframe,"_",idet,"_",termin.c_str(),".bi");
		fp = fopen(testfile,"r");
		fread(samptopix,sizeof(long),ns,fp);
		fclose(fp);


		deproject(S,indpix,samptopix,ns,nn,npix,Ps,flgdupl,factdupl);


		data[ii] = data[ii] - Ps[ii];

		MapMakPreProcessData(data,flag,calp,ns,napod,4,1.0,data_lp,bfilter,
				NORMLIN,NOFILLGAP);

		for (ii=0;ii<ns;ii++)
			data[ii] = data_lp[ii]*apodwind[ii];


		// compute fft and save data to disk for later
		fftplan = fftw_plan_dft_r2c_1d(ns, data, fdata1, FFTW_ESTIMATE);
		fftw_execute(fftplan);
		fftw_destroy_plan(fftplan);


		sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"fdata_",ff,"_",idet,"_",termin.c_str(),".bi");
		fp = fopen(testfile,"w");
		fwrite(fdata1,sizeof(double), (ns/2+1)*2, fp);
		fclose(fp);



		/// compute sigma of the noise
		mm = 0.0;
		for (ii=ns/2;ii<ns/2+500;ii++) mm += data[ii];
		mm = mm/501.0;
		tmpsign = 0.0;
		for (ii=ns/2;ii<ns/2+500;ii++) tmpsign += (data[ii]-mm)*(data[ii]-mm);
		sign[idet] = sqrt(tmpsign/500.0);
		// normalize to the first detector
		if (idet == 0) sign0 = sign[0];
		sign[idet] = sign[idet]/sign0;



		for (jj=0;jj<ncomp;jj++)
			for (ii=0;ii<ns;ii++)
				commonm[jj][ii] += aa[idet][jj]/(sign[idet]*sign[idet])*data[ii];


	}

	for (jj=0; jj<ncomp; jj++)
		for (ii= 0 ;ii<ns;ii++)
			if ( isnan(commonm[jj][ii]) || isinf(commonm[jj][ii]) )
				cout << "Something went wrong in the computation of the common mode " << ncomp << endl;


	//***************************************************************************


	/////////// AtN-1A
	for (jj=0;jj<ncomp;jj++)
		for (kk=0;kk<ncomp;kk++)
			for (ii=0;ii<ndet;ii++)
				Cov[jj][kk] += aa[ii][jj] * aa[ii][kk]/sign[ii]/sign[ii];


	// invert AtN-1A
	dcholdc(Cov,ncomp,p);
	for (ii=0;ii<ncomp;ii++){
		for (jj=0;jj<ncomp;jj++)
			uvec[jj] = 0.0;
		uvec[ii] = 1.0;
		dcholsl(Cov,ncomp,p,uvec,ivec);
		for (jj=0;jj<ncomp;jj++)
			iCov[ii][jj] = ivec[jj];
	}


	for (ii=0;ii<ns;ii++)
		for (jj=0;jj<ncomp;jj++)
			for (kk=0;kk<ncomp;kk++)
				commonm2[jj][ii] += iCov[jj][kk] * commonm[kk][ii];



	factapod = 0.0;
	for (ii=0;ii<ns;ii++)
		factapod += apodwind[ii]*apodwind[ii]/ns;



	/// filter common mode to fcut Hz
	for (ii=0;ii<ncomp;ii++){
		for (jj=0;jj<ns;jj++)
			commontmp[jj] = commonm2[ii][jj];
		//butterworth(commontmp,ns,fcut/fsamp*double(ns),8,commonm_f,bfilter,1,0,0);
		for (jj=0;jj<ns;jj++)
			common_f[ii][jj] = commontmp[jj];                         //// - commonm_f[jj];
	}



	//**************************************** Read pre-estimated power spectrum for reference
	sprintf(nameSpfile,"%s%s%s%s",noiseSppreffile.c_str(),field.c_str(),"-all",extentnoiseSp.c_str());
	if ((fp = fopen(nameSpfile,"r")) == NULL){
		cerr << "ERROR: Can't find noise power spectra file" << nameSpfile << " , check -k or -K in command line. Exiting. \n";
		exit(1);
	}
	fread(&dnbins,sizeof(double), 1, fp);
	nbins = (long)dnbins;
	nbins2 = nbins;
	SpN_all = dmatrix(0,ndet-1,0,nbins-1);
	ell = new double[nbins+1];
	SpN = new double[nbins];
	fread(ell,sizeof(double), nbins+1, fp);
	fread(*SpN_all,sizeof(double), nbins*ndet, fp);
	fclose(fp);
	delete [] SpN;
	//*****************************************

	//    printf("Inside EstimPowerSpectra after 1st step\n");



	//   ///////////////////////////////////// 2nd step: reestimate mixing matrix

	//   for (ii=0;ii<ncomp;ii++)
	//     for (jj=0;jj<ncomp;jj++)
	//       Cov[ii][jj] = 0.0;

	//   for (ii=0;ii<ncomp;ii++)
	//     for (jj=0;jj<ndet;jj++)
	//       vect[ii][jj] = 0.0;

	//   for (idet=0;idet<ndet;idet++){

	//     // read data from disk
	//     field = bolonames[idet];
	//     read_data_std(dirfile, ff, shift_data_to_point, ns, data, field+bextension, 'd');

	//     for (ii=0;ii<ncomp;ii++)
	//       for (jj=0;jj<ns;jj++)
	// 	vect[ii][idet] += common_f[ii][jj] * data[jj]*apodwind[jj];


	//     ///// measure power spectrum
	//     noisepectrum_estim(data,ns,ell,nbins,fsamp,NULL,Nell,Nk);


	//     ///// write power spectrum to disk
	//     tempstr2 = bolonames[idet];
	//     sprintf(nameSpfile,"%s%s%s%s%s%ld%s",outdirSpN.c_str(),tempstr2.c_str(),"-",tempstr2.c_str(),"_",ff,"exp.psd");
	//     fp = fopen(nameSpfile,"w");
	//     fprintf(fp,"%d\n",nbins);
	//     for (ii=0;ii<nbins;ii++){
	//       fprintf(fp,"%g\t",ell[ii]);
	//       fprintf(fp,"%10.15g\n",Nell[ii]);
	//     }
	//     fprintf(fp,"%g\n",ell[nbins]);
	//     fprintf(fp,"\n");
	//     fclose(fp);


	//   }

	//   for (ii=0;ii<ncomp;ii++)
	//     for (jj=0;jj<ncomp;jj++)
	//       for (kk=0;kk<ns;kk++)
	// 	Cov[ii][jj] += common_f[ii][kk] * common_f[jj][kk];


	//   // invert AtN-1A
	//   dcholdc(Cov,ncomp,p);
	//   for (ii=0;ii<ncomp;ii++){
	//     for (jj=0;jj<ncomp;jj++)
	//       uvec[jj] = 0.0;
	//     uvec[ii] = 1.0;
	//     dcholsl(Cov,ncomp,p,uvec,ivec);
	//     for (jj=0;jj<ncomp;jj++)
	//       iCov2[ii][jj] = ivec[jj];
	//   }


	//   //update mixing matrix

	//   for (ii=0;ii<ndet;ii++){
	//     for (jj=0;jj<ncomp;jj++){
	//       aa[ii][jj] = 0.0;
	//       for (kk=0;kk<ncomp;kk++)
	// 	aa[ii][jj] += iCov2[jj][kk] * vect[kk][ii];
	//     }
	//   }



	//   //////////////////////////////////////
	//   // write new mixing parameters to disk
	//   sprintf(testfile,"%s%s%d%s%s",outdirSpN.c_str(),"MixMat_",ff,termin.c_str(),".txt");
	//   fp = fopen(testfile,"w");
	//   fprintf(fp,"%d \n",ncomp);
	//   for (ii=0;ii<ndet;ii++){
	//     for (jj=0;jj<ncomp;jj++){
	//       fprintf(fp,"%lf \t",aa[ii][jj]);
	//     }
	//     fprintf(fp,"\n");
	//   }









	//************************************************************************//
	// second part: -- data - common mode
	//              -- estimate noise power spectra


	/////////////////////////////////////// loop again over detectors
	for (idet=0;idet<ndet;idet++){

		field = bolonames[idet];


		read_data_std(dirfile, ff, shift_data_to_point, ns, data, field+bextension, 'd');

		if (fextension != "NOFLAG"){
			read_data_std(dirfile, ff, shift_data_to_point, ns, flag, field+fextension,  'c');
		} else {
			//      printf("NOFLAG\n");
			for (ii=0;ii<ns;ii++)
				flag[ii] = 0;
		}

		if (cextension != "NOCALP"){
			read_data_std(dirfile, ff, 0, ns/20, calp, field+cextension, 'd');
		} else {
			//      printf("NOCALP\n");
			for (ii=0;ii<ns/20;ii++)
				calp[ii] = 1.0;
		}



		//******************************* subtract signal

		//Read pointing data
		sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"samptopix_",iframe,"_",idet,"_",termin.c_str(),".bi");
		fp = fopen(testfile,"r");
		fread(samptopix,sizeof(long),ns,fp);
		fclose(fp);

		deproject(S,indpix,samptopix,ns,nn,npix,Ps,flgdupl,factdupl);


		data[ii] = data[ii] - Ps[ii];




		MapMakPreProcessData(data,flag,calp,ns,napod,4,1.0,data_lp,bfilter,
				NORMLIN,NOFILLGAP);



		for (ii=0;ii<ns;ii++)
			data[ii] = data_lp[ii] * apodwind[ii];


		// Subtract components
		for (ii=0;ii<ns;ii++)
			for (jj=0;jj<ncomp;jj++)
				data[ii] -= aa[idet][jj]*common_f[jj][ii];


		//Noise power spectra

		///// measure power spectrum of the uncorrelated part of the noise
		noisepectrum_estim(data,ns,ell,nbins,fsamp,NULL,Nell,Nk);


		for (ii=0;ii<nbins;ii++){
			Rellth[idet*ndet+idet][ii] += Nell[ii]/factapod;
			N[idet][ii] = Nell[ii]/factapod;
		}
	}




	////*********************** Component power spectra

	for (ii=0;ii<ncomp;ii++){
		for (jj=0;jj<ns;jj++)
			commontmp[jj]=common_f[ii][jj];
		noisepectrum_estim(commontmp,ns,ell,nbins,fsamp,NULL,Nell,Nk);
		for (jj=0;jj<nbins;jj++)
			P[ii][jj] = Nell[jj]/factapod;

		// subtract a factor to correct from noise
		//for (jj=0;jj<nbins;jj++)
		//  P[ii][jj] -= iCov[ii][ii]*sign0*sign0;
		//for (jj=0;jj<nbins;jj++)
		//  if (P[ii][jj] < 0)
		//	P[ii][jj] = 0.0;


		printf("iCov[%ld] = %10.15g\n",ii,iCov[ii][ii]*sign0*sign0);


		sprintf(testfile,"%s%s%ld%s%ld%s%s%s",outdirSpN.c_str(),"Nellc_",ii,"_",ff,"_",termin.c_str(),".bi");
		fp = fopen(testfile,"w");
		fwrite(Nell,sizeof(double), nbins, fp);
		fclose(fp);


	}
	printf("noise var det 0 =  %10.15g\n",sign0*sign0);




	for (ii=0;ii<ndet;ii++)
		for (kk=0;kk<ndet;kk++)
			for (ll=0;ll<ncomp;ll++)
				for (jj=0;jj<nbins;jj++)
					Rellth[ii*ndet+kk][jj] += aa[ii][ll] * aa[kk][ll] * P[ll][jj];







	//************************************************************************//
	// third part: -- estimate the covariance matrix of the data R_exp
	//



	/////////////////////////////////////// loop again over detectors
	for (idet1=0;idet1<ndet;idet1++){

		// read data from disk
		sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"fdata_",ff,"_",idet1,"_",termin.c_str(),".bi");
		fp = fopen(testfile,"r");
		fread(fdata1,sizeof(double), (ns/2+1)*2, fp);
		fclose(fp);


		for (idet2=0;idet2<ndet;idet2++) {

			// read data from disk
			sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"fdata_",ff,"_",idet2,"_",termin.c_str(),".bi");
			fp = fopen(testfile,"r");
			fread(fdata2,sizeof(double), (ns/2+1)*2, fp);
			fclose(fp);


			noisecrosspectrum_estim(fdata1,fdata2,ns,ell,nbins,fsamp,NULL,Nell,Nk);


			for (ii=0;ii<nbins;ii++)
				Rellexp[idet1*ndet+idet2][ii] += Nell[ii]/factapod;

		}

		cout << "Computing Rellexp :" << idet1*100./ndet << "%\r" << flush ;

	}

	////// normalize to the first detector power spectrum in order to avoid numerical problems
	for (ii=0;ii<nbins;ii++)
		SPref[ii] = Rellexp[0][ii];

	for (ii=0;ii<nbins;ii++)
		if ( isnan(SPref[ii]) || isinf(SPref[ii]))
			cout << " Problem in the first detector power spectrum\n";

	for (idet1=0;idet1<ndet;idet1++)
		for (idet2=0;idet2<ndet;idet2++)
			for (ii=0;ii<nbins;ii++)
				Rellexp[idet1*ndet+idet2][ii] = Rellexp[idet1*ndet+idet2][ii]/SPref[ii];
	for (jj=0;jj<ncomp;jj++)
		for (ii=0;ii<nbins;ii++)
			P[jj][ii] = P[jj][ii]/SPref[ii];
	for (jj=0;jj<ndet;jj++)
		for (ii=0;ii<nbins;ii++)
			N[jj][ii] = N[jj][ii]/SPref[ii];




	printf("Before saving data to disk\n");



	//// write Rellexp to disk and also first guess of parameters
	sprintf(testfile,"%s%s%d%s%s",outdirSpN.c_str(),"Rellexp_",(int)ff,termin.c_str(),".txt");
	fp = fopen(testfile,"w");
	for (jj=0;jj<nbins;jj++)
		for (ii=0;ii<ndet;ii++)
			for (kk=0;kk<ndet;kk++)
				fprintf(fp,"%10.15g \t",Rellexp[ii*ndet+kk][jj]);
	fprintf(fp,"\n");
	fclose(fp);

	sprintf(testfile,"%sNinit_%d_%s.txt",outdirSpN.c_str(),(int)ff,termin.c_str());
	fp = fopen(testfile,"w");
	for (ii=0;ii<ndet;ii++){
		for (jj=0;jj<nbins;jj++)
			fprintf(fp,"%10.15g \t",N[ii][jj]);
		fprintf(fp,"\n");
	}
	fclose(fp);

	data1d = new double[ndet*nbins];
	for (long i=0; i< ndet; i++)
		for (long j=0; j<nbins; j++)
			data1d[i*nbins+j] = N[i][j];

	sprintf(testfile,"!%sNinit_%d_%s.fits",outdirSpN.c_str(),(int)ff,termin.c_str());
	write_psd_tofits(testfile,ndet,nbins,'d', data1d);
	delete [] data1d;


	sprintf(testfile,"%s%s%d%s%s",outdirSpN.c_str(),"Pinit_",(int)ff,termin.c_str(),".txt");
	fp = fopen(testfile,"w");
	for (ii=0;ii<ncomp;ii++)
		for (jj=0;jj<nbins;jj++)
			fprintf(fp,"%10.15g \t",P[ii][jj]);
	fprintf(fp,"\n");
	fclose(fp);


	sprintf(testfile,"%s%s%d%s%s",outdirSpN.c_str(),"Ainit_",(int)ff,termin.c_str(),".txt");
	fp = fopen(testfile,"w");
	for (ii=0;ii<ndet;ii++)
		for (jj=0;jj<ncomp;jj++){
			fprintf(fp,"%10.15g \t",aa[ii][jj]);
			fprintf(fp,"\n");
		}
	fclose(fp);









	//***** Fourth part
	//*********************** fit component and noise power spectra, and mixing matrix *************//
	//********* Using Expectation/Maximization algorithm
	long nbiter = 500;
	long iter;

	double tottest=0.0;

	double f;

	double *iN, *Pr, *w;
	double **Rxs, **Rxsq, **RnRxsb, **Rxx, **Rxxq, **Rss, **Rssq, **RnRssb;
	double **Pr2, **AiNA, **Mattmp, **ImDR, **ACq, **Cq, **Wq;


	ib = 0;
	while ((ell[ib] < fcut) && (ib < nbins)){
		nbins2 = ib+1;
		ib++;
	}

	/////// to be removed:
	//nbins2 = nbins;

	printf("nbins2 = %d                \n",nbins2);



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
	for (ib=0;ib<nbins2;ib++)
		w[ib] = (ell[ib+1] - ell[ib])*ns/fsamp;



	sprintf(testfile,"%s%s%d%s%s",outdirSpN.c_str(),"weights_",(int)ff,termin.c_str(),".txt");
	fp = fopen(testfile,"w");
	for (jj=0;jj<nbins2;jj++){
		fprintf(fp,"%10.15g \t",w[jj]);
	}
	fprintf(fp,"\n");
	fclose(fp);


	f = fdsf(Rellexp,w,aa,P,N,ndet,ncomp,nbins2) ;
	printf("Pre em:   obj: %10.15g\n", f) ;


	for (iter=1;iter<=nbiter;iter++){

		init1D_double(iN,0,ndet,0.0);
		init1D_double(Pr,0,ncomp,0.0);
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

		for (ib=0;ib<nbins2;ib++){

			for (idet=0;idet<ndet;idet++)
				iN[idet] = 1.0/N[idet][ib];

			for (idet1=0;idet1<ndet;idet1++)
				for (idet2=0;idet2<ndet;idet2++)
					Rxxq[idet1][idet2] = Rellexp[idet1*ndet + idet2][ib];

			// Robust wrt Pq=0
			for (ii=0;ii<ncomp;ii++)
				Pr[ii] = sqrt(P[ii][ib]);

			for (ii=0;ii<ncomp;ii++)
				for (jj=0;jj<ncomp;jj++)
					Pr2[ii][jj] = Pr[ii]*Pr[jj];

			for (ii=0;ii<ncomp;ii++)
				for (jj=0;jj<ncomp;jj++){
					AiNA[ii][jj] = 0.0;
					for (idet=0;idet<ndet;idet++)
						AiNA[ii][jj] += aa[idet][jj] * aa[idet][ii] * iN[idet];
				}

			for (ii=0;ii<ncomp;ii++){
				for (jj=0;jj<ncomp;jj++)
					Cov[ii][jj] = Pr2[ii][jj] * AiNA[ii][jj];
				Cov[ii][ii] += 1.0;
			}

			// invert matrix
			dcholdc(Cov,ncomp,p);

			for (ii=0;ii<ncomp;ii++){
				for (jj=0;jj<ncomp;jj++)
					uvec[jj] = 0.0;
				uvec[ii] = 1.0;
				dcholsl(Cov,ncomp,p,uvec,ivec);
				for (jj=0;jj<ncomp;jj++)
					iCov[ii][jj] = ivec[jj];
			}


			for (ii=0;ii<ncomp;ii++)
				for (jj=0;jj<ncomp;jj++)
					Cq[ii][jj] = Pr2[ii][jj] * iCov[ii][jj];

			for (idet=0;idet<ndet;idet++)
				for (ii=0;ii<ncomp;ii++){
					Wq[idet][ii] = 0.0;
					for (jj=0;jj<ncomp;jj++)
						Wq[idet][ii] += aa[idet][jj] * iN[idet] *Cq[jj][ii] ;
				}
			for (idet=0;idet<ndet;idet++)
				for (ii=0;ii<ncomp;ii++){
					Rxsq[idet][ii] = 0.0;
					for (jj=0;jj<ndet;jj++)
						Rxsq[idet][ii] += Rxxq[idet][jj] * Wq[jj][ii];
				}

			for (kk=0;kk<ncomp;kk++)
				for (ii=0;ii<ncomp;ii++){
					Rssq[kk][ii] = Cq[kk][ii];
					for (jj=0;jj<ndet;jj++)
						Rssq[kk][ii] += Rxsq[jj][kk] * Wq[jj][ii];
				}

			for (kk=1;kk<ncomp;kk++)
				for (ii=0;ii<kk;ii++){
					Rssq[ii][kk] = 0.5*(Rssq[ii][kk]+Rssq[kk][ii]);
					Rssq[kk][ii] = Rssq[ii][kk];
				}


			// update power spectra
			for (ii=0;ii<ncomp;ii++)
				P[ii][ib] = abs(Rssq[ii][ii]);



			for (ii=0;ii<ndet;ii++)
				for (jj=0;jj<ncomp;jj++)
					RnRxsb[ii][jj] += w[ib] * iN[ii]*Rxsq[ii][jj];



			for (kk=0;kk<ncomp;kk++)
				for (ii=0;ii<ndet;ii++)
					for (jj=0;jj<ncomp;jj++)
						RnRssb[kk][jj+ii*ncomp] = RnRssb[kk][jj+ii*ncomp] + w[ib] * iN[ii] * Rssq[kk][jj] ;

		}



		// update mixing matrix
		for (idet=0;idet<ndet;idet++){

			for (ii=0;ii<ncomp;ii++){
				uvec[ii] = RnRxsb[idet][ii];
				for (jj=0;jj<ncomp;jj++)
					Cov[ii][jj] = RnRssb[ii][jj+idet*ncomp];
			}

			// solving the linear system
			dcholdc(Cov,ncomp,p);
			dcholsl(Cov,ncomp,p,uvec,ivec);
			for (ii=0;ii<ncomp;ii++)
				aa[idet][ii] = ivec[ii];
		}



		// EM Step with respect to N, with the new values of A and P
		for (ib=0;ib<nbins2;ib++){

			for (idet = 0;idet<ndet;idet++)
				iN[idet] = 1.0/N[idet][ib];

			for (idet1=0;idet1<ndet;idet1++)
				for (idet2=0;idet2<ndet;idet2++)
					Rxxq[idet1][idet2] = Rellexp[idet1*ndet + idet2][ib];

			// Robust wrt Pq=0
			for (ii=0;ii<ncomp;ii++)
				Pr[ii] = sqrt(P[ii][ib]);

			for (ii=0;ii<ncomp;ii++)
				for (jj=0;jj<ncomp;jj++)
					Pr2[ii][jj] = Pr[ii]*Pr[jj];

			for (ii=0;ii<ncomp;ii++)
				for (jj=0;jj<ncomp;jj++){
					AiNA[ii][jj] = 0.0;
					for (idet=0;idet<ndet;idet++)
						AiNA[ii][jj] += aa[idet][jj] * aa[idet][ii] * iN[idet];
				}

			for (ii=0;ii<ncomp;ii++){
				for (jj=0;jj<ncomp;jj++)
					Cov[ii][jj] = Pr2[ii][jj] * AiNA[ii][jj];
				Cov[ii][ii] += 1.0;
			}



			// invert matrix
			dcholdc(Cov,ncomp,p);
			for (ii=0;ii<ncomp;ii++){
				for (jj=0;jj<ncomp;jj++)
					uvec[jj] = 0.0;
				uvec[ii] = 1.0;
				dcholsl(Cov,ncomp,p,uvec,ivec);
				for (jj=0;jj<ncomp;jj++)
					iCov[ii][jj] = ivec[jj];
			}


			for (ii=0;ii<ncomp;ii++)
				for (jj=0;jj<ncomp;jj++)
					Cq[ii][jj] = Pr2[ii][jj] * iCov[ii][jj];

			for (idet=0;idet<ndet;idet++)
				for (ii=0;ii<ncomp;ii++){
					ACq[idet][ii] = 0.0;
					for (jj=0;jj<ncomp;jj++)
						ACq[idet][ii] += aa[idet][jj]*Cq[jj][ii];
				}

			for (idet=0;idet<ndet;idet++)
				for (jj=0;jj<ndet;jj++){
					ImDR[jj][idet] = 0.0;
					if (jj == idet)
						ImDR[jj][idet] = 1.0;
					for (ii=0;ii<ncomp;ii++)
						ImDR[jj][idet] -= ACq[jj][ii]*iN[idet]*aa[idet][ii] ;
				}


			for (idet=0;idet<ndet;idet++)
				for (ii=0;ii<ndet;ii++){
					Mattmp[idet][ii] = 0.0;
					for (jj=0;jj<ndet;jj++)
						Mattmp[idet][ii] += Rellexp[idet+jj*ndet][ib] * ImDR[ii][jj];
				}

			for (idet=0;idet<ndet;idet++){
				N[idet][ib] = 0.0;
				for (ii=0;ii<ndet;ii++)
					N[idet][ib] += ImDR[idet][ii] * Mattmp[ii][idet];
			}


			for (idet=0;idet<ndet;idet++){
				for (ii=0;ii<ncomp;ii++){
					N[idet][ib] += ACq[idet][ii]*aa[idet][ii];
				}
				N[idet][ib] = abs(N[idet][ib]);
			}


		}


		tottest = 0.0;
		for (ib=0;ib<nbins2;ib++){
			for (idet=0;idet<ndet;idet++)
				tottest += N[idet][ib]/ndet;
			for (idet=0;idet<ndet;idet++)
				if (N[idet][ib] < tottest*1e-8)
					N[idet][ib] = tottest*1e-8;
		}




		//printf("A[1][2] =  %10.15g\n", aa[1][2]) ;
		//printf("N[2][3] =  %10.15g\n", N[2][3]) ;
		//printf("P[2][3] =  %10.15g\n", P[2][3]) ;


		// Fixing the indeterminacies.  Is it useful here?
		//    rescaleAP(aa, P, ndet, ncomp, nbins2) ;



		///// here is the problem

		f = fdsf(Rellexp,w,aa,P,N,ndet,ncomp,nbins2) ;
		cout << "em->iter: " << iter*100./nbiter << " %\r" << flush;
		if (isnan(f) || isinf(f)) {
			cout << "Nan........." << endl;
			exit(1);
		}


	}


	// Fixing the indeterminacies.  Is it useful here?
	rescaleAP(aa, P, ndet, ncomp, nbins2) ;




	//****************************** Compute covariance matrix from the fitted model


	printf("EM step completed\n");


	for (jj=0;jj<nbins;jj++)
		for (ii=0;ii<ndet*ndet;ii++)
			Rellth[ii][jj] = 0.0;

	for (jj=0;jj<nbins2;jj++){
		for (idet=0;idet<ndet;idet++){
			Rellth[idet*ndet+idet][jj] += N[idet][jj]*SPref[jj];
			for (ii=0;ii<ndet;ii++)
				for (ll=0;ll<ncomp;ll++)
					Rellth[idet*ndet+ii][jj] += aa[idet][ll] * aa[ii][ll] * P[ll][jj]*SPref[jj];
		}
	}


	if (nbins2 < nbins)
		for (jj=nbins2;jj<nbins;jj++)
			for (idet=0;idet<ndet;idet++)
				Rellth[idet*ndet+idet][jj] = Rellexp[idet*ndet+idet][jj]*SPref[jj];



	//*****************  Write power spectra to disk  ********************//


	sprintf(nameSpfile,"%s%s%d%s%s",outdirSpN.c_str(),"BoloPS",(int)ff,termin.c_str(),".psd");
	fp = fopen(nameSpfile,"w");
	for (idet1=0;idet1<ndet;idet1++){
		for (idet2=0;idet2<ndet;idet2++){

			///// write power spectrum to disk
			tempstr1 = bolonames[idet1];
			tempstr2 = bolonames[idet2];
			//sprintf(nameSpfile,"%s%s%s%s%s%ld%s",outdirSpN.c_str(),tempstr1.c_str(),"-",tempstr2.c_str(),"_",ff,".psd");
			//fp = fopen(nameSpfile,"w");
			fprintf(fp,"%s%s%s\n",tempstr1.c_str(),"-",tempstr2.c_str());
			fprintf(fp,"%d\n",nbins);
			for (ii=0;ii<nbins;ii++){
				fprintf(fp,"%g\t",ell[ii]);
				fprintf(fp,"%10.15g\n",(Rellth[idet1*ndet+idet2][ii]+Rellth[idet2*ndet+idet1][ii])/2.0);
			}
			fprintf(fp,"%g\n",ell[nbins]);
			//fprintf(fp,"\n");
			//fclose(fp);
		}
	}
	fclose(fp);

	// New format for output of power spectra

	sprintf(nameSpfile,"%s%s%d%s%s",outdirSpN.c_str(),"BoloPS",(int)ff,termin.c_str(),"_binary.psd");
	fp = fopen(nameSpfile,"w");
	fwrite(&ndet,sizeof(long),1,fp);
	fwrite(&nbins,sizeof(long),1,fp);
	fwrite(ell,sizeof(double),nbins+1,fp);
	for (idet1=0;idet1<ndet;idet1++)
		for (idet2=0;idet2<ndet;idet2++)
			for (ii=0;ii<nbins; ii++)
				fwrite(&Rellth[idet1*ndet+idet2][ii],sizeof(double),1,fp);
	fclose(fp);

	// write ell
	sprintf(nameSpfile,"%s%s%d%s%s",outdirSpN.c_str(),"Ell_",(int)ff,termin.c_str(),".psd");
	fp = fopen(nameSpfile,"w");
	for (ii=0;ii<nbins;ii++){
		fprintf(fp,"%g\n",ell[ii]);
	}
	fclose(fp);


	sprintf(nameSpfile,"%s%s%d%s%s",outdirSpN.c_str(),"BoloPS",(int)ff,termin.c_str(),"_exp.psd");
	fp = fopen(nameSpfile,"w");

	for (idet1=0;idet1<ndet;idet1++){
		for (idet2=0;idet2<ndet;idet2++){

			///// write power spectrum to disk
			tempstr1 = bolonames[idet1];
			tempstr2 = bolonames[idet2];
			//sprintf(nameSpfile,"%s%s%s%s%s%ld%s",outdirSpN.c_str(),tempstr1.c_str(),"-",tempstr2.c_str(),"_",ff,"_exp.psd");
			//fp = fopen(nameSpfile,"w");
			fprintf(fp,"%s%s%s\n",tempstr1.c_str(),"-",tempstr2.c_str());
			fprintf(fp,"%d\n",nbins);
			for (ii=0;ii<nbins;ii++){
				fprintf(fp,"%g\t",ell[ii]);
				fprintf(fp,"%10.15g\n",(Rellexp[idet1*ndet+idet2][ii]+Rellexp[idet2*ndet+idet1][ii])/2.0*SPref[ii]);
			}
			fprintf(fp,"%g\n",ell[nbins]);
			//fprintf(fp,"\n");
			//fclose(fp);
		}
	}
	fclose(fp);



	sprintf(testfile,"%s%s%d%s%s",outdirSpN.c_str(),"Afinal_",(int)ff,termin.c_str(),".txt");
	fp = fopen(testfile,"w");
	for (ii=0;ii<ndet;ii++)
		for (jj=0;jj<ncomp;jj++)
			fprintf(fp,"%10.15g \t",aa[ii][jj]);
	fprintf(fp,"\n");
	fclose(fp);



	//**************** Write component power spectra to disk
	for (idet1=0;idet1<ndet;idet1++){

		tempstr1 = bolonames[idet1];
		sprintf(nameSpfile,"%s%s%s%ld%s",outdirSpN.c_str(),tempstr1.c_str(),"_uncnoise",ff,".psd");
		fp = fopen(nameSpfile,"w");
		//fprintf(fp,"%d\n",nbins);
		for (ii=0;ii<nbins;ii++){
			//fprintf(fp,"%g\t",ell[ii]);
			fprintf(fp,"%10.15g\n",N[idet1][ii]*SPref[ii]);
		}
		//fprintf(fp,"%g\n",ell[nbins]);
		fprintf(fp,"\n");
		fclose(fp);
	}

	data1d = new double[ndet*nbins];
	for (long i=0; i< ndet; i++)
		for (long j=0; j<nbins; j++)
			data1d[i*nbins+j] = N[i][j]*SPref[j];

	sprintf(testfile,"!%sNfinal_%d_%s_uncnoise.fits",outdirSpN.c_str(),(int)ff,termin.c_str());
	write_psd_tofits(testfile,nbins,ndet,'d',data1d);
	delete [] data1d;



	for (jj=0;jj<ncomp;jj++){

		sprintf(nameSpfile,"%s%s%ld%s%ld%s",outdirSpN.c_str(),"Comp_",jj,"_uncnoise",ff,".psd");
		fp = fopen(nameSpfile,"w");
		//fprintf(fp,"%d\n",nbins);
		for (ii=0;ii<nbins;ii++){
			//fprintf(fp,"%g\t",ell[ii]);
			fprintf(fp,"%10.15g\n",P[jj][ii]*SPref[ii]);
		}
		//fprintf(fp,"%g\n",ell[nbins]);
		fprintf(fp,"\n");
		fclose(fp);
	}

	data1d = new double[ncomp*nbins];
	for (long i=0; i< ncomp; i++)
		for (long j=0; j<nbins; j++)
			data1d[i*nbins+j] = P[i][j]*SPref[j];

	sprintf(testfile,"!%sNfinal_%d_%s_cnoise.fits",outdirSpN.c_str(),(int)ff,termin.c_str());
	write_psd_tofits(testfile,nbins,ncomp,'d',data1d);
	delete [] data1d;




	delete [] data;
	delete [] data_lp;
	delete [] fdata1;
	delete [] fdata2;
	delete [] Ps;
	delete [] samptopix;
	delete [] commontmp;
	delete [] commonm_f;
	delete [] calp;
	delete [] flag;
	delete [] bfilter;
	delete [] apodwind;
	delete [] Nell;
	free_dmatrix(Rellexp,0,ndet*ndet-1,0,nbins-1);
	free_dmatrix(Rellth,0,ndet*ndet-1,0,nbins-1);
	free_dmatrix(aa,0,ndet-1,0,20);
	delete [] sign;
	free_dmatrix(Cov,0,ncomp-1,0,ncomp-1);
	free_dmatrix(iCov,0,ncomp-1,0,ncomp-1);
	free_dmatrix(iCov2,0,ncomp-1,0,ncomp-1);
	delete [] p;
	delete [] uvec;
	delete [] ivec;
	free_dmatrix(commonm,0,ncomp,0,ns-1);
	free_dmatrix(commonm2,0,ncomp,0,ns-1);
	free_dmatrix(common_f,0,ncomp,0,ns-1);
	free_dmatrix(vect,0,ncomp-1,0,ndet-1);
	delete [] Nk;
	delete [] ell;
	free_dmatrix(SpN_all,0,ndet-1,0,nbins-1);
	free_dmatrix(P,0,ncomp-1,0,nbins-1);
	free_dmatrix(N,0,ndet-1,0,nbins-1);
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

	long ib, ii, jj, kk;
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
	init1D_double(Pl,0,ncomp,0.0);
	init1D_double(Pnl,0,ndet,0.0);
	init2D_double(iRhR,0,0,ndet,ndet,0.0);


	// init
	f   = 0. ;

	for (ib=0;ib<nbins;ib++){

		//// reading Rexp
		for (ii=0;ii<ncomp;ii++)
			Pl[ii] = P[ii][ib];
		for (ii=0;ii<ndet;ii++){
			Pnl[ii] = N[ii][ib];
			for (jj=0;jj<ndet;jj++)
				hR[ii][jj] = Rellexp[ii+jj*ndet][ib];
		}

		/// building Rth from A, P, N
		for (ii=0;ii<ndet;ii++)
			for (jj=0;jj<ndet;jj++){
				R[ii][jj] = 0.0;
				if (ii == jj)
					R[ii][jj] = Pnl[ii];
				for (kk=0;kk<ncomp;kk++)
					R[ii][jj] += A[ii][kk] * Pl[kk] * A[jj][kk];
			}


		//printf("ib=%d\n",ib);
		/// inverting Rth
		for (ii=0;ii<ndet;ii++)
			for (jj=0;jj<ndet;jj++)
				eR[ii][jj] = R[ii][jj];
		dcholdc(eR,ndet,p);
		for (ii=0;ii<ndet;ii++){
			for (jj=0;jj<ndet;jj++)
				uvec[jj] = 0.0;
			uvec[ii] = 1.0;
			dcholsl(eR,ndet,p,uvec,ivec);
			for (jj=0;jj<ndet;jj++)
				iR[ii][jj] = ivec[jj];
		}
		// printf("ib=%d\n",ib);


		/// computing mismatch from Rexp and Rth
		for (ii=0;ii<ndet;ii++)
			for (jj=0;jj<ndet;jj++){
				iRhR[ii][jj] = 0.0;
				for (kk=0;kk<ndet;kk++)
					iRhR[ii][jj] += iR[ii][kk]*hR[kk][jj] ;
			}


		triRhR = 0.0;
		logdetiR = 0;
		for (ii=0;ii<ndet;ii++){
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

	long ii, jj, ib;
	double *norm2ratio;

	norm2ratio = new double[ncomp];

	init1D_double(norm2ratio,0,ncomp,0.0);

	for (ii=0;ii<ncomp;ii++){
		for (jj=0;jj<ndet;jj++)
			norm2ratio[ii] += A[jj][ii] * A[jj][ii];
		norm2ratio[ii] = 1.0/norm2ratio[ii];
	}

	for (ii=0;ii<ncomp;ii++){
		for (jj=0;jj<ndet;jj++)
			A[jj][ii] = A[jj][ii] * sqrt(norm2ratio[ii]) ;
		for (ib=0;ib<nbins;ib++)
			P[ii][ib] = P[ii][ib] / norm2ratio[ii] ;
	}

	delete [] norm2ratio;

}
