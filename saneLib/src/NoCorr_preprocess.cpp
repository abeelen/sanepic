/*
 * NoCorr_preprocess.cpp
 *
 *  Created on: 28 mai 2009
 *      Author: matthieu
 */
#include <iostream>
#include <vector>
#include <string>

#include <fftw3.h>
#include <gsl/gsl_math.h>

#include "todprocess.h"
#include "map_making.h"
#include "inline_IO2.h"
#include "dataIO.h"
#include "NoCorr_preprocess.h"


using namespace std;

//void do_PtNd_nocorr(double *PNd, string *extentnoiseSp_all, string noiseSppreffile,
//		string dir, string dirfile,
//		std::vector<string> bolonames, string *fits_table,
//		double f_lppix, double f_lppix_Nk,
//		double fsamp, long ntotscan, long addnpix, bool flgdupl, int factdupl,
//		int fillg, long ns, long napod, long ndet,
//		long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2, long long npix,
//		long long npixsrc, bool NORMLIN, bool NOFILLGAP,bool remove_polynomia, long iframe, double *S)

void do_PtNd_nocorr(double *PNd,string tmp_dir, struct user_options u_opt,
		struct samples samples_struct, struct input_commons com, struct detectors det, double f_lppix, double f_lppix_Nk,
		long addnpix, long ns, long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2, long long npix,
		long long npixsrc, long iframe, double *S)
{


	string field;
	string extentnoiseSp;

	string nameSpfile;


	long long *samptopix;
	double *bfilter, *Nk, *data, *data_lp, *Ps;
	short *flag, *flpoint;
	double powered;

	samptopix = new long long[ns];
	bfilter = new double[ns/2+1];
	Nk = new double[ns/2+1];


	data =  new double[ns];
	data_lp = new double[ns];
	flpoint = new short[ns];
	Ps = new double[ns];

	int factdupl = 1;
	if (com.flgdupl) factdupl=2;


	string fits_filename;

	fits_filename = samples_struct.fits_table[iframe];


	//for (long idet=rank*ndet/size;idet<(rank+1)*ndet/size;idet++){
	for (long idet=0;idet<det.ndet;idet++){
		field = det.boloname[idet];



		if (S != NULL){
			read_flpoint_from_fits(fits_filename, flpoint);
		}

		read_signal_from_fits(fits_filename, data, field);


		long test_ns;
		read_flag_from_fits(fits_filename , field, flag, test_ns);
		// TODO : Test sizes...



		//// Read pointing
		read_samptopix(ns, samptopix, tmp_dir, idet, iframe,det.boloname);


		if (S != NULL){

			if (addnpix){
				deproject(S,indpix,samptopix,ns,NAXIS1, NAXIS2,npix,Ps,2,factdupl,samples_struct.ntotscan,indpsrc,npixsrc);
			} else {
				deproject(S,indpix,samptopix,ns,NAXIS1, NAXIS2,npix,Ps,2,factdupl);
			}

		}


		if (S != NULL){
			//********************  pre-processing of data ********************//
			MapMakPreProcessData(data,flag,ns,com.napod,4,f_lppix,data_lp,bfilter,
					u_opt.NORMLIN,com.NOFILLGAP,u_opt.remove_polynomia,Ps);
		}
		else {
			MapMakPreProcessData(data,flag,ns,com.napod,4,f_lppix,data_lp,bfilter,
					u_opt.NORMLIN,com.NOFILLGAP,u_opt.remove_polynomia);
		}


		/************************************************************************************/

		// end of preprocess begin of fdata

		for (long ii=0;ii<ns/2+1;ii++){
			powered=pow(double(ii)/f_lppix, 16);
			bfilter[ii] = powered /(1.0+powered);
		}



		//****************** Compute (or read) input power spectrum of the NOISE  ***************//
		extentnoiseSp = samples_struct.noise_table[iframe];
		nameSpfile = tmp_dir + field + extentnoiseSp;

		readNSpectrum(nameSpfile,bfilter,ns,u_opt.fsamp,Nk);


		//********************** compute P^t N-1 d ************************//
		compute_PtNmd(data_lp,Nk,ns,NAXIS1, NAXIS2,indpix,samptopix,npix,PNd);


	}// end of idet loop


	delete[] samptopix;
	delete[] bfilter;
	delete[] Nk;

	delete[] data;
	delete[] data_lp;

	delete[] flag;
	delete[] flpoint;

	delete[] Ps;


}





//void do_PtNPS_nocorr(double *S, string *extentnoiseSp_all, string noiseSppreffile, string dir,
//		/* string termin, */ string dirfile, std::vector<string> bolonames, double f_lppix,
//		double fsamp, bool flgdupl, int factdupl, long ns,
//		long ndet, long long *indpix, long NAXIS1, long NAXIS2, long long npix,
//		long iframe, double *PtNPmatS, double *Mp, long *hits)

void do_PtNPS_nocorr(double *S, string *extentnoiseSp_all, struct directories dir,
		struct detectors det,double f_lppix,double fsamp, bool flgdupl, long ns,
		long long *indpix, long NAXIS1, long NAXIS2, long long npix,
		long iframe, double *PtNPmatS, double *Mp, long *hits)

{

	string field;
	string extentnoiseSp;
	string nameSpfile;


	long long *samptopix;
	double *bfilter, *Nk, *Ps;
	double powered;

	samptopix = new long long[ns];
	bfilter = new double[ns/2+1];
	Nk = new double[ns/2+1];
	Ps = new double[ns];


	int factdupl = 1;
	if(flgdupl==1) factdupl = 2;


	//for (long idet=rank*ndet/size;idet<(rank+1)*ndet/size;idet++){
	for (long idet=0;idet<det.ndet;idet++){
		field = det.boloname[idet];

		read_samptopix(ns, samptopix, dir.tmp_dir, idet, iframe,det.boloname);


		// AS
		deproject(S,indpix,samptopix,ns,NAXIS1, NAXIS2,npix,Ps,flgdupl,factdupl);




		for (long ii=0;ii<(ns)/2+1;ii++){
			powered=pow(double(ii)/f_lppix, 16);
			bfilter[ii] = powered /(1.0+powered);
		}

		extentnoiseSp = extentnoiseSp_all[iframe];
		//sprintf(nameSpfile,"%s%s%s",noiseSppreffile.c_str(),field.c_str(),extentnoiseSp.c_str());
		nameSpfile = dir.tmp_dir + field + extentnoiseSp;
		readNSpectrum(nameSpfile,bfilter,ns,fsamp,Nk);


		//AtN-1A AS (espensive part)
		compute_PtNmd(Ps,Nk,ns,NAXIS1, NAXIS2,indpix,samptopix,npix,PtNPmatS);


		//Compute weight map for preconditioner
		if ((Mp != NULL))
			compute_diagPtNP(Nk,samptopix,ns,NAXIS1, NAXIS2,indpix,npix,f_lppix,Mp);


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













