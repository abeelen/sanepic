/*
 * NoCorr_preprocess.cpp
 *
 *  Created on: 28 mai 2009
 *      Author: matthieu
 */

#include <vector>
#include <string>

#include "NoCorr_preprocess.h"
#include "inline_IO.h"

using namespace std;

void do_PtNd_nocorr(double *PNd, string *extentnoiseSp_all, string noiseSppreffile,
		string dir, string termin, double errarcsec, string dirfile,
		string scerr_field, string flpoint_field, std::vector<string> bolonames,
		string bextension, string fextension, /*string cextension,*/
		int shift_data_to_point, double f_lppix, double f_lppix_Nk,
		double fsamp, long ntotscan, long addnpix, bool flgdupl, int factdupl,
		int fillg, long ff, long ns, long napod, long ndet,
		int size, int rank, long *indpix, long *indpsrc, long nn, long npix,
		long npixsrc, bool NORMLIN, bool NOFILLGAP,bool remove_polynomia, long iframe, double *S){



	//long ii, idet;
	//long ndata = ns+2*marge;
	string field;
	string extentnoiseSp;

	char nameSpfile[100];
	//char testfile[100];

	long *samptopix;
	double *bfilter, *Nk, *data, *data_lp, *scerr,/* *calp,*/ *Ps;
	unsigned char *flag, *flpoint, *rejectsamp;


	samptopix = new long[ns];
	bfilter = new double[ns/2+1];
	Nk = new double[ns/2+1];

	scerr = new double[ns];
	data =  new double[ns];
	data_lp = new double[ns];
	//calp =  new double[ns];
	flag =  new unsigned char[ns];
	flpoint = new unsigned char[ns];
	rejectsamp = new unsigned char[ns];
	Ps = new double[ns];


	//FILE *fp;



	for (long idet=rank*ndet/size;idet<(rank+1)*ndet/size;idet++){

		field = bolonames[idet];



		if (S != NULL){
			read_data_std(dirfile, ff, 0, ns, scerr, scerr_field, 'd');
			read_data_std(dirfile, ff, 0, ns, flpoint, flpoint_field, 'c');
		}

		read_data_std(dirfile, ff, shift_data_to_point, ns, data, field+bextension, 'd');

		if (fextension != "NOFLAG"){
			read_data_std(dirfile, ff, shift_data_to_point, ns, flag, field+fextension,  'c');
		} else {
			//      printf("NOFLAG\n");
			for (long ii=0;ii<ns;ii++)
				flag[ii] = 0;
		}

		//if (cextension != "NOCALP"){
		//	read_data_std(dirfile, ff, 0, ns/20, calp, field+cextension, 'd');
		//} else {
		//      printf("NOCALP\n");
		//for (ii=0;ii<ns/20;ii++)
		//	calp[ii] = 1.0;
		//}


		//// Read pointing
		read_samptopix(ns, samptopix, termin, dir, idet, iframe);

		/*sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"samptopix_",iframe,"_",idet,"_",termin.c_str(),".bi");
		fp = fopen(testfile,"r");
		fread(samptopix,sizeof(long),ns,fp);
		fclose(fp);*/


		if (S != NULL){

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

// end of preprocess begin of fdata

		for (long ii=0;ii<ns/2+1;ii++)
			bfilter[ii] = pow(double(ii)/f_lppix_Nk, 16) /(1.0+pow(double(ii)/f_lppix_Nk, 16));



		//****************** Compute (or read) input power spectrum of the NOISE  ***************//
		extentnoiseSp = extentnoiseSp_all[iframe];
		//nameSpfile = noiseSppreffile + field + extentnoiseSp;
		sprintf(nameSpfile,"%s%s%s",noiseSppreffile.c_str(),field.c_str(),extentnoiseSp.c_str());
		readNSpectrum(nameSpfile,bfilter,ns,fsamp,Nk);


		//********************** compute P^t N-1 d ************************//
		compute_PtNmd(data_lp,Nk,ns,nn,indpix,samptopix,npix,PNd);


	}// end of idet loop


	delete[] samptopix;
	delete[] bfilter;
	delete[] Nk;
	delete[] scerr;
	delete[] data;
	delete[] data_lp;
	//delete[] calp;
	delete[] flag;
	delete[] flpoint;
	delete[] rejectsamp;
	delete[] Ps;


}





void do_PtNPS_nocorr(double *S, string *extentnoiseSp_all, string noiseSppreffile, string dir,
		string termin, string dirfile, std::vector<string> bolonames, double f_lppix,
		double fsamp, bool flgdupl, int factdupl, long ff, long ns,
		long ndet, int size, int rank, long *indpix, long nn, long npix,
		long iframe, double *PtNPmatS, double *Mp, long *hits){



	//long ii, idet;
	//long ndata = ns;
	string field;
	string extentnoiseSp;

	char nameSpfile[100];
	//char testfile[100];

	long *samptopix;
	double *bfilter, *Nk, *Ps;


	samptopix = new long[ns];
	bfilter = new double[ns/2+1];
	Nk = new double[ns/2+1];
	Ps = new double[ns];

	//FILE *fp;



	for (long idet=rank*ndet/size;idet<(rank+1)*ndet/size;idet++){

		field = bolonames[idet];

		read_samptopix(ns, samptopix, termin, dir, idet, iframe);
		/*sprintf(testfile,"%s%s%ld%s%ld%s%s%s",dir.c_str(),"samptopix_",iframe,"_",idet,"_",termin.c_str(),".bi");
		fp = fopen(testfile,"r");
		fread(samptopix,sizeof(long),ns,fp);
		fclose(fp);*/


		// AS
		deproject(S,indpix,samptopix,ns,nn,npix,Ps,flgdupl,factdupl);




		for (long ii=0;ii<(ns)/2+1;ii++)
			bfilter[ii] = pow(double(ii)/f_lppix, 16) /(1.0+pow(double(ii)/f_lppix, 16));

		extentnoiseSp = extentnoiseSp_all[iframe];
		sprintf(nameSpfile,"%s%s%s",noiseSppreffile.c_str(),field.c_str(),extentnoiseSp.c_str());
		readNSpectrum(nameSpfile,bfilter,ns,fsamp,Nk);


		//AtN-1A AS (espensive part)
		compute_PtNmd(Ps,Nk,ns,nn,indpix,samptopix,npix,PtNPmatS);


		//Compute weight map for preconditioner
		if ((Mp != NULL))
			compute_diagPtNP(Nk,samptopix,ns,nn,indpix,npix,f_lppix,Mp);


		//compute hit counts
		if (hits != NULL){
			for (long ii=0;ii<ns;ii++){
				hits[indpix[samptopix[ii]]] += 1;
			}
		}


	} // end of idet loop



	delete[] samptopix;
	delete[] bfilter;
	delete[] Nk;
	delete[] Ps;


}













