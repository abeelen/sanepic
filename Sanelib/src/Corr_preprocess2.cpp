/*
 * Corr_preprocess.cpp
 *
 *  Created on: 20 juil. 2009
 *      Author: matthieu
 */


#include <string>
#include <vector>

#include "Corr_preprocess2.h"


using namespace std;

// TODO : S = signal ???
void write_tfAS2(double *S, long *indpix, int nn, long npix, bool flgdupl, int factdupl, string dir, string termin, long ff, long ns, long ndet, long iframe){


	//long idet1;
	//long ndata = ns+2*marge;

	//FILE *fp;
	//char testfile[100];

	double *Ps; // TODO PS = ?
	long *samptopix;

	fftw_plan fftplan;
	fftw_complex *fdata;

	samptopix = new long[ns];
	Ps = new double[ns];
	fdata = new fftw_complex[ns/2+1];


	//for (idet1=rank*ndet/size;idet1<(rank+1)*ndet/size;idet1++){
	for (long idet1=0;idet1<ndet;idet1++){

		//Read pointing data
		read_samptopix(ns, samptopix, termin, dir, idet1, iframe);
		/*sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"samptopix_",iframe,"_",idet1,"_",termin.c_str(),".bi");
		fp = fopen(testfile,"r");
		fread(samptopix,sizeof(long),ns,fp);
		fclose(fp);*/

		deproject(S,indpix,samptopix,ns,nn,npix,Ps,flgdupl,factdupl);

		//Fourier transform of the data
		fftplan = fftw_plan_dft_r2c_1d(ns, Ps, fdata, FFTW_ESTIMATE);
		fftw_execute(fftplan);
		fftw_destroy_plan(fftplan);


		write_fPs(ns, fdata, termin, dir, idet1, iframe);
		/*sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"fPs_",iframe,"_",idet1,"_",termin.c_str(),".bi");
		fp = fopen(testfile,"w");
		fwrite(fdata,sizeof(double), (ns/2+1)*2, fp);
		fclose(fp);*/

	}

	delete[] samptopix;
	delete[] Ps;
	delete[] fdata;

}






void write_ftrProcesdata2(double *S, long *indpix, long *indpsrc, int nn, long npix,
		long npixsrc, long ntotscan, long addnpix, bool flgdupl, int factdupl,
		int fillg, string dir, string termin, double errarcsec, string dirfile,
		string scerr_field, string flpoint_field, std::vector<string> bolonames,
		string bextension, string fextension, /*string cextension,*/
		int shift_data_to_point, double f_lppix, long ff, long ns,
		long napod, long ndet, bool NORMLIN, bool NOFILLGAP, bool remove_polynomia, long iframe,fftw_complex *&fdata_buffer/*,fftw_complex **&fdatas*/){


	//long ii, idet1;
	//long ndata = ns+2*marge;

	double *scerr, *data,/* *calp,*/ *bfilter, *data_lp, *Ps;
	unsigned char *flpoint, *flag, *rejectsamp;
	long *samptopix;

	fftw_plan fftplan;
	fftw_complex *fdata;
	//fftw_complex *fdata_buffer;


	string field1;

	//char testfile[100];
	//char testfile2[100];

	//FILE *fp,*fp2;

	scerr = new double[ns];
	data =  new double[ns];
	data_lp = new double[ns];
	//calp =  new double[ns];
	flag =  new unsigned char[ns];
	flpoint = new unsigned char[ns];
	rejectsamp = new unsigned char[ns];

	samptopix = new long[ns];
	Ps = new double[ns];
	bfilter = new double[ns/2+1];
	fdata = new fftw_complex[ns/2+1];



	//std::vector<double*> fdatas_vect[ndet];


	//sprintf(testfile,"%s%s",dir.c_str(),"test_file_fdatas.txt");
	//sprintf(testfile2,"%s%s",dir.c_str(),"test_file_fdatas_buffer.txt");
	//fp=fopen(testfile,"a");
	//fp2=fopen(testfile2,"a");

	for (long idet1=0;idet1<ndet;idet1++){

		field1 = bolonames[idet1];
		//     cout << "Inside write_ftrProcessdata " << endl;
		//    cout << " field 1 : " << field1 << endl;

		if (S != NULL){
			read_data_std(dirfile, ff, 0, ns, scerr, scerr_field, 'd');
			read_data_std(dirfile, ff, 0, ns, flpoint, flpoint_field, 'c');
		}

		read_data_std(dirfile, ff, shift_data_to_point, ns, data, field1+bextension, 'd');

		if (fextension != "NOFLAG"){
			read_data_std(dirfile, ff, shift_data_to_point, ns, flag, field1+fextension,  'c');
		} else {
			//      printf("NOFLAG\n");
			for (long ii=0;ii<ns;ii++)
				flag[ii] = 0;
		}

		//if (cextension != "NOCALP"){
		//read_data_std(dirfile, ff, 0, ns/20, calp, field1+cextension, 'd'); // attention avec le samples_per_frame !!!
		//read_data_std(dirfile, ff, 0, ns/samples_per_frame, calp, field1+cextension, 'd'); // attention avec le samples_per_frame !!!
		//} else {
		//      printf("NOCALP\n");
		//for (ii=0;ii<ns/20;ii++)
		//for (ii=0;ii<ns/samples_per_frame;ii++) // attention 20 avant
		//calp[ii] = 1.0;
		//}



		if (S != NULL){
			//// Read pointing
			read_samptopix(ns, samptopix, termin, dir, idet1, iframe);

			/*sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"samptopix_",iframe,"_",idet1,"_",termin.c_str(),".bi");
			fp = fopen(testfile,"r");
			fread(samptopix,sizeof(long),ns,fp);
			fclose(fp);*/

			if (addnpix){
				deproject(S,indpix,samptopix,ns,nn,npix,Ps,fillg,factdupl,ntotscan,indpsrc,npixsrc);
			} else {
				deproject(S,indpix,samptopix,ns,nn,npix,Ps,fillg,factdupl);
			}

			for (long ii=0;ii<ns;ii++) rejectsamp[ii] = 0;
			for (long ii=0;ii<ns;ii++)
				if ((flag[ii] & 1) != 0 || (scerr[ii] > errarcsec) || (flpoint[ii] & 1) != 0)
					rejectsamp[ii] = 1;
		}


		if (S != NULL){
			//********************  pre-processing of data ********************//
			MapMakPreProcessData(data,flag,/*calp,*/ns,napod,4,f_lppix,data_lp,bfilter,
					NORMLIN,NOFILLGAP,remove_polynomia,Ps);
		}
		else {
			MapMakPreProcessData(data,flag,/*calp,*/ns,napod,4,f_lppix,data_lp,bfilter,
					NORMLIN,NOFILLGAP,remove_polynomia);
		}

		//Fourier transform of the data
		fftplan = fftw_plan_dft_r2c_1d(ns, data_lp, fdata, FFTW_ESTIMATE);
		fftw_execute(fftplan);
		fftw_destroy_plan(fftplan);

		for(long ii=0;ii<ns/2+1;ii++){
			fdata_buffer[idet1*(ns/2+1)+ii][0]=fdata[ii][0];
			fdata_buffer[idet1*(ns/2+1)+ii][1]=fdata[ii][1];
		}

		//write fourier transform to disk
		//write_fdata(ns, fdata, termin, dir, idet1, iframe);


		//fprintf(fp,"%lf, %lf\n",fdata[0][0],fdata[0][1]);

		/*sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"fdata_",iframe,"_",idet1,"_",termin.c_str(),".bi");
		fp = fopen(testfile,"w");
		fwrite(fdata,sizeof(double), (ns/2+1)*2, fp);
		fclose(fp);*/

	}

	//fclose(fp);

	delete[] scerr;
	delete[] data;
	delete[] data_lp;
	//delete[] calp;
	delete[] flag;
	delete[] flpoint;
	delete[] samptopix;
	delete[] Ps;
	delete[] bfilter;
	delete[] fdata;
	delete[] rejectsamp;

	/*for(long idet=0;idet<ndet;idet++){
		fprintf(fp2,"%lf, %lf\n",fdata_buffer[idet*(ns/2+1)][0],fdata_buffer[idet*(ns/2+1)][1]);
	}*/


//	fclose(fp2);
	//cout << "fin" << endl;
	//getchar();

}



void do_PtNd2(double *PNd, string *extentnoiseSp_all, string noiseSppreffile,
		string dir, string prefixe, string termin, std::vector<string> bolonames,
		double f_lppix, double fsamp, long ff, long ns, long ndet, int size,
		int rank, long *indpix, long nn, long npix, long iframe,/*fftw_complex **fdatas,*/ double *Mp, long *hits, fftw_complex *fdata_buffer){


	long  nbins;
	//double dnbins;
	//long ndata = ns+2*marge;
	string field1, field2;
	string extentnoiseSp;

	string nameSpfile;
	//char testfile[100];
	//FILE *fp;


	long *samptopix;
	double *ell, *SpN, *bfilter, *bfilter_, *Nk, *Nd;

	double powered;

	fftw_plan fftplan;
	fftw_complex *fdata, *Ndf;


	samptopix = new long[ns];
	Nd = new double[ns];
	bfilter = new double[ns/2+1];
	bfilter_ = new double[ns/2+1];
	Nk = new double[ns/2+1];
	fdata = new fftw_complex[ns/2+1];
	Ndf = new fftw_complex[ns/2+1];





	double **SpN_all;

	//FILE *fp;

	//cout << ndet << endl;
	//cout << rank << " " << size << endl;

	for (long idet1=rank*ndet/size;idet1<(rank+1)*ndet/size;idet1++){
		field1 = bolonames[idet1];



		//Read pointing data
		read_samptopix(ns, samptopix, termin, dir, idet1, iframe);
		//cout << "samptopix" << endl;

		//**************************************** Noise power spectrum
		extentnoiseSp = extentnoiseSp_all[iframe];
		nameSpfile = noiseSppreffile + field1 + "-all" + extentnoiseSp;
		//cout << nameSpfile << endl;

		//		sprintf(nameSpfile,"%s%s%s%s",noiseSppreffile.c_str(),field1.c_str(),"-all",extentnoiseSp.c_str());

		//read noise PS file
		read_noise_file(nbins, ell, SpN_all, nameSpfile, ndet);
		SpN = new double[nbins];


		//cout << "apres read" << endl;

		//*****************************************

		for (long ii=0;ii<ns/2+1;ii++){
			powered=gsl_pow_int(double(ii)/f_lppix, 16);
			bfilter[ii] = powered /(1.0+powered);
		}
		for (long ii=0;ii<ns/2+1;ii++)
			bfilter_[ii] = 1.0/(bfilter[ii]+0.000001);


		//Init N-1d
		for (long ii=0;ii<ns/2+1;ii++){
			Ndf[ii][0] = 0;
			Ndf[ii][1] = 0;
		}



		for (long idet2=0;idet2<ndet;idet2++){
			//field2 = bolonames[idet2];
			//cout << field2 << endl;


			//read Fourier transform of the data
			//read_fdata(ns, fdata, prefixe, termin, dir, idet2, iframe);
			for(long ii=0;ii<ns/2+1;ii++){
				fdata[ii][0]=fdata_buffer[(ns/2+1)*idet2+ii][0];
				fdata[ii][1]=fdata_buffer[(ns/2+1)*idet2+ii][1];
			}

			//cout << "apres fdata" << endl;



			//****************** Cross power spectrum of the noise  ***************//
			for (int ii=0;ii<nbins;ii++)
				SpN[ii] = SpN_all[idet2][ii];


			// interpolate logarithmically the noise power spectrum
			InvbinnedSpectrum2log_interpol(ell,SpN,bfilter_,nbins,ns,fsamp,Nk);


			for (long jj=0;jj<ns/2+1;jj++)
				if (isnan(Nk[jj])) {
					printf("isnan has been found : iframe %ld, det1 %ld, det2 %ld\n",iframe, idet1, idet2);
					exit(1);
				}

			//********************************* compute N^-1 d  ***********************//
			for (long ii=0;ii<ns/2+1;ii++){
				Ndf[ii][0] += fdata[ii][0]*Nk[ii];
				Ndf[ii][1] += fdata[ii][1]*Nk[ii];
			}



			//Compute weight map for preconditioner
			if ((Mp != NULL) && (idet2 == idet1))
				compute_diagPtNPCorr(Nk,samptopix,ns,nn,indpix,npix,f_lppix,Mp);
			//


		}// end of idet2 loop

		//cout << "fftw" << endl;
		fftplan = fftw_plan_dft_c2r_1d(ns, Ndf, Nd, FFTW_ESTIMATE);
		fftw_execute(fftplan);
		fftw_destroy_plan(fftplan);
		//cout << "apres fftw" << endl;


		for (long ii=0;ii<ns;ii++){
			if ((ii < 0) || (ii >= ns)){
				PNd[npix-2] += Nd[ii];
			} else {
				PNd[indpix[samptopix[ii]]] += Nd[ii]; // Nd real
			}
		}

		//compute hit counts
		if (hits != NULL){
			for (long ii=0;ii<ns;ii++){
				hits[indpix[samptopix[ii]]] += 1;
			}
		}

		//cout << "avant delete" << endl;
		delete[] ell;
		delete[] SpN;
		free_dmatrix(SpN_all,0,ndet-1,0,nbins-1);
		//cout << "delete" << endl;


	}// end of idet1 loop


	delete[] samptopix;
	delete[] Nd;
	delete[] bfilter;
	delete[] bfilter_;
	delete[] Nk;
	delete[] fdata;
	delete[] Ndf;


}

