/*
 * Corr_preprocess.cpp
 *
 *  Created on: 20 juil. 2009
 *      Author: matthieu
 */

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

#include "Corr_preprocess.h"
#include "covMatrixIO.h"
#include "dataIO.h"
#include "inline_IO2.h"
#include "todprocess.h"
#include "map_making.h"


//temp
#include <fstream>

#include <gsl/gsl_math.h>
#include <fftw3.h>
//#include <time.h>

using namespace std;


//void write_tfAS(double *S, long long *indpix, long NAXIS1, long NAXIS2, long long npix,
//		bool flgdupl, string dir, long ns, long ndet, long iframe, std::vector<string> bolonames)

void write_tfAS(double *S, struct detectors det,long long *indpix, long NAXIS1, long NAXIS2, long long npix,
		bool flgdupl, string dir, long ns, long iframe)
{



	double *Ps;
	long long *samptopix;

	fftw_plan fftplan;
	fftw_complex *fdata;

	samptopix = new long long[ns];
	Ps = new double[ns];
	fdata = new fftw_complex[ns/2+1];

	fill(samptopix,samptopix+ns,0);
	fill(Ps,Ps+ns,0.0);
	for(long ii =0; ii<ns/2+1;ii++){
		fdata[ii][0]=0.0;
		fdata[ii][1]=0.0;
	}
	// map duplication factor
	int factdupl;
	(flgdupl) ? factdupl = 2: factdupl = 1; //  default 1 : if flagged data are put in a duplicated map

	//for (idet1=rank*ndet/size;idet1<(rank+1)*ndet/size;idet1++){
	for (long idet1=0;idet1<det.ndet;idet1++){

		//Read pointing data
		read_samptopix(ns, samptopix, dir, idet1, iframe, det.boloname);

		//		cout << "samptopix : " << endl;
		//		cout << samptopix[0] << " " << samptopix[1] << " " << samptopix[2] << endl;

		// temporary down
		deproject(S,indpix,samptopix,ns,NAXIS1, NAXIS2,npix,Ps,flgdupl,factdupl);

		//Fourier transform of the data
		fftplan = fftw_plan_dft_r2c_1d(ns, Ps, fdata, FFTW_ESTIMATE);
		fftw_execute(fftplan);
		fftw_destroy_plan(fftplan);


		write_fPs(ns, fdata, dir, idet1, iframe, det.boloname);

	}

	delete[] samptopix;
	delete[] Ps;
	delete[] fdata;

}


//void write_ftrProcesdata(double *S, long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2, long long npix,
//		long long npixsrc, long ntotscan, long long addnpix, bool flgdupl, int factdupl,
//		int fillg, string dir, string dirfile,
//		std::vector<string> bolonames,string *fits_table, double f_lppix, long ns,
//		long napod, long ndet, bool NORMLIN, bool NOFILLGAP, bool remove_polynomia,
//		long iframe)

void write_ftrProcesdata(double *S, struct user_options u_opt, struct samples samples_struct, struct input_commons com,
		string tmp_dir,	struct detectors det, long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2,
		long long npix,	long long npixsrc, long long addnpix, double f_lppix, long ns, long iframe)
{


	//long ii, idet1;
	//long ndata = ns+2*marge;

	double *data, *bfilter, *data_lp, *Ps;
	short *flag;
	//	short *flpoint;
	long long *samptopix;

	fftw_plan fftplan;
	fftw_complex *fdata;


	string field1, fits_filename;

	//char testfile[100];

	//FILE *fp;
	//	cout << "avant les alloc" << endl;

	//scerr = new double[ns];
//	data =  new double[ns];
	data_lp = new double[ns];
	//calp =  new double[ns];
	//flag =  new unsigned char[ns];
	//flpoint = new unsigned char[ns];
	//flag = new short[ns];
	//	flpoint = new short[ns];
	//rejectsamp = new unsigned char[ns];

	samptopix = new long long[ns];
	Ps = new double[ns];
	bfilter = new double[ns/2+1];
	fdata = new fftw_complex[ns/2+1];

	//	unsigned char *flag;
	//	double *data2;
	//	data2= new double[ns];
	//	flag = new unsigned char[ns];

//	cout << "avant les fill" << endl;

//	fill(flag,flag+ns,0);
//	fill(data,data+ns,0.0);
	fill(data_lp,data_lp+ns,0.0);
	fill(Ps,Ps+ns,0.0);
	fill(bfilter,bfilter+(ns/2+1),0.0);
	fill(samptopix,samptopix+ns,0);
	//	fill(flpoint,flpoint+ns,0);

	for (long ii=0;ii<ns/2+1;ii++){
		fdata[ii][0] = 0.0;
		fdata[ii][1] = 0.0;
	}

//	cout << "apres les fill" << endl;

	int factdupl = 1;
	if(com.flgdupl==1)		factdupl = 2;

	fits_filename = samples_struct.fits_table[iframe];
	cout << "fits file : " << fits_filename << endl;

	for (long idet1=0;idet1<det.ndet;idet1++){

		field1 = det.boloname[idet1];
//				cout << field1 << endl;


		//		fill(data,data+ns,0.0);
		//		fill(data_lp,data_lp+ns,0.0);
		//		fill(Ps,Ps+ns,0.0);
		//		fill(bfilter,bfilter+ns,0.0);
		//		fill(samptopix,samptopix+ns,0);
		//		fill(flpoint,flpoint+ns,0);

		for (long ii=0;ii<ns/2+1;ii++){
			fdata[ii][0] = 0.0;
			fdata[ii][1] = 0.0;
		}
		//		cout << field1 << "  apres ALLOC "  << endl;

//		if (S != NULL){
//			// TODO: What is the point of this ?? Do we really need this ?
//			//			read_flpoint_from_fits(fits_filename, flpoint);
//			//cout << "flpoint : " << flpoint[0] <<  flpoint[1] << flpoint[2] << flpoint[3] << endl;
//		}

		long test_ns;
		read_signal_from_fits(fits_filename, field1, data, test_ns);
		if (test_ns != ns) {
			cerr << "Read signal does not correspond to frame size : Check !!" << endl;
			exit(-1);
		}

		read_flag_from_fits(fits_filename , field1, flag, test_ns);
		if (test_ns != ns) {
			cerr << "Read flag does not correspond to frame size : Check !!" << endl;
			exit(-1);
		}


		//********************  pre-processing of data ********************//

		if (S != NULL){

			//TODO : Ps should not be here...  remove the signal before or make the deproject inside MapMakePreProcess
			//TODO : write fdata inside MapMakePreProcess.. or create a function same is true in sanePS

			cout << "S !=NULL\n";
			// Read pointing
			read_samptopix(ns, samptopix, tmp_dir, idet1, iframe, det.boloname);
			//TODO : Fix that... same number of argument... not the same calling as in sanePS
			// Deproject
			if (addnpix)
				deproject(S,indpix,samptopix,ns,NAXIS1, NAXIS2,npix,Ps,2,factdupl,samples_struct.ntotscan,indpsrc,npixsrc);
			else
				deproject(S,indpix,samptopix,ns,NAXIS1, NAXIS2,npix,Ps,2,factdupl);

			//********************  pre-processing of data ********************//
			MapMakPreProcessData(data,flag,ns,com.napod,u_opt.poly_order,f_lppix,data_lp,bfilter, // default poly order = 4
					u_opt.NORMLIN,com.NOFILLGAP,u_opt.remove_polynomia,Ps);
		}
		else {
			MapMakPreProcessData(data,flag,ns,com.napod,u_opt.poly_order,f_lppix,data_lp,bfilter, // default poly order = 4
					u_opt.NORMLIN,com.NOFILLGAP,u_opt.remove_polynomia);
		}

		//cout << "data apres map : " << setprecision(14)  << data_lp[0] << " " << data_lp[1] << " " << data_lp[2] << " "  << data_lp[3] << endl;

		//Fourier transform of the data
		fftplan = fftw_plan_dft_r2c_1d(ns, data_lp, fdata, FFTW_ESTIMATE);
		fftw_execute(fftplan);
		fftw_destroy_plan(fftplan);

		//		cout << "write fdata_" << iframe << "_" << idet1 << endl;

		//write fourier transform to disk
		write_fdata(ns, fdata, tmp_dir, idet1, iframe, det.boloname);

		delete [] data;
		delete [] flag;
		//		cout << "write fdata_" << iframe << "_" << idet1 << endl;
	}


	//	cout << "avant les clean" << endl;

//	delete[] data;
	delete[] data_lp;
	//	delete[] flag;
	//	delete[] flpoint;
	delete[] samptopix;
	delete[] Ps;
	delete[] bfilter;
	delete[] fdata;

//	cout << "apres les clean" << endl;
//	getchar();
}

void do_PtNd(double *PNd, string *extentnoiseSp_all, string dir, string prefixe,
		struct detectors det, double f_lppix, double fsamp, long ns,
		long long *indpix, long NAXIS1, long NAXIS2, long long npix, long iframe,
		double *Mp, long *hits)

{

	long  nbins;
	string field1, field2;
	string extentNoiseSp;

	string nameSpfile;

	long long *samptopix;
	double *ell, *SpN, *bfilter, *bfilter_, *Nk, *Nd;

	double powered;

	fftw_plan fftplan;
	fftw_complex *fdata, *Ndf;


	samptopix = new long long[ns];
	Nd = new double[ns];
	bfilter = new double[ns/2+1];
	bfilter_ = new double[ns/2+1];
	Nk = new double[ns/2+1];
	fdata = new fftw_complex[ns/2+1];
	Ndf = new fftw_complex[ns/2+1];

	double **SpN_all;

	//TODO : This is a butterworth filter.... why not use butterworth()
	for (long ii=0;ii<ns/2+1;ii++){
		powered=gsl_pow_int(double(ii)/f_lppix,16);
		bfilter[ii] = powered /(1.0+powered);
	}
	for (long ii=0;ii<ns/2+1;ii++)
		bfilter_[ii] = 1.0/(bfilter[ii]+0.000001);

	fill(Nd,Nd+ns,0.0);
	fill(Nk,Nk+(ns/2+1),0.0);
	fill(samptopix,samptopix+ns,0);

	//	for (long ii=0;ii<ns/2+1;ii++){
	//		fdata[ii][0] = 0.0;
	//		fdata[ii][1] = 0.0;
	//	}





	//cout << rank << " " << size << endl;

	//for (long idet1=rank*ndet/size;idet1<(rank+1)*ndet/size;idet1++){
	for (long idet1=0;idet1<det.ndet;idet1++){
		field1 = det.boloname[idet1];

		//Read pointing data
		read_samptopix(ns, samptopix, dir, idet1, iframe, det.boloname);


		//**************************************** Noise power spectrum
		extentNoiseSp = extentnoiseSp_all[iframe];
		nameSpfile = dir + field1 + "-all" + extentNoiseSp;

		//read noise PS file
		long ndet2;
		read_InvNoisePowerSpectra(dir, field1,  extentNoiseSp, &nbins, &ndet2, &ell, &SpN_all);
		if(det.ndet!=ndet2) cout << "Error. The number of detector in noisePower Spectra file must be egal to input bolofile number\n";

		SpN = new double[nbins];
		fill(SpN,SpN+nbins,0.0);

		//Init N-1d
		for (long ii=0;ii<ns/2+1;ii++){
			Ndf[ii][0] = 0.0;
			Ndf[ii][1] = 0.0;
		}

		//		if(idet1>0)
		//			cout << "avant double boucle" << endl;

		for (long idet2=0;idet2<det.ndet;idet2++){
			field2 = det.boloname[idet2];

			fill(Nd,Nd+ns,0.0);
			fill(Nk,Nk+(ns/2+1),0.0);

			//			for (long ii=0;ii<ns/2+1;ii++){
			//				fdata[ii][0] = 0.0;
			//				fdata[ii][1] = 0.0;
			//			}

			//read Fourier transform of the data
			read_fdata(ns, fdata, prefixe, dir, idet2, iframe, det.boloname);


			//****************** Cross power spectrum of the noise  ***************//
			for (int ii=0;ii<nbins;ii++){
				SpN[ii] = SpN_all[idet2][ii];
				//cout << SpN[ii] << " ";
			}


			// TODO : Why do we need to reinterpolate the noise power spectrum here ?
			// interpolate logarithmically the noise power spectrum
			InvbinnedSpectrum2log_interpol(ell,SpN,bfilter_,nbins,ns,fsamp,Nk);
			//InvbinnedSpectrum2bis(ell,SpN,bfilter_,nbins,ns,fsamp,Nk);

			for (long jj=0;jj<ns/2+1;jj++){
				//				filee << Nk[jj] << endl;
				if (isnan(Nk[jj])) {
					printf("isnan has been found : iframe %ld, det1 %ld, det2 %ld\n",iframe, idet1, idet2);
					exit(1);
				}
			}

			//			filee.close();
			//********************************* compute N^-1 d  ***********************//
			for (long ii=0;ii<ns/2+1;ii++){
				Ndf[ii][0] += fdata[ii][0]*Nk[ii];
				Ndf[ii][1] += fdata[ii][1]*Nk[ii];
			}



			//Compute weight map for preconditioner
			if ((Mp != NULL) && (idet2 == idet1))
				compute_diagPtNPCorr(Nk,samptopix,ns,NAXIS1, NAXIS2,indpix,npix,f_lppix,Mp);
			//


		}// end of idet2 loop

		fftplan = fftw_plan_dft_c2r_1d(ns, Ndf, Nd, FFTW_ESTIMATE);
		fftw_execute(fftplan);
		fftw_destroy_plan(fftplan);

		for (long ii=0;ii<ns;ii++){
			//			if ((ii < 0) || (ii >= ns)){
			//				PNd[npix-2] += Nd[ii];
			//			} else {
			PNd[indpix[samptopix[ii]]] += Nd[ii]; // Nd real
			//			}
		}

		//compute hit counts
		if (hits != NULL){
			for (long ii=0;ii<ns;ii++){
				hits[indpix[samptopix[ii]]] += 1;
			}
		}


		delete[] ell;
		delete[] SpN;
		free_dmatrix(SpN_all,0,det.ndet-1,0,nbins-1);


	}// end of idet1 loop

	delete[] samptopix;
	delete[] Nd;
//	delete[] bfilter_; // For some reason there is a double free on this one....
	delete[] bfilter;
	delete[] Nk;
	delete[] fdata;
	delete[] Ndf;

}
