#include <iostream>
#include <vector>
#include <string>

#include <fftw3.h>
#include <gsl/gsl_math.h>

#include "todprocess.h"
#include "map_making.h"
#include "temporary_IO.h"
#include "dataIO.h"
#include "NoCorr_preprocess.h"


using namespace std;

void do_PtNd_nocorr(double *PNd,string tmp_dir, struct param_process proc_param, struct param_positions pos_param,
		struct samples samples_struct, struct detectors det, double f_lppix, double f_lppix_Nk,
		long addnpix, long ns, long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2, long long npix,
		long long npixsrc, long iframe, double *S, int rank, int size)
{


	string field;
	string extentnoiseSp;

	string nameSpfile;


	long long *samptopix;
	double *bfilter, *Nk, *data, *data_lp, *Ps;
	int *flag=NULL;
	double powered;

	samptopix = new long long[ns];
	bfilter = new double[ns/2+1];
	Nk = new double[ns/2+1];


	data_lp = new double[ns];
	Ps = new double[ns];

	int factdupl = 1;
	if (pos_param.flgdupl) factdupl=2;


	string fits_filename;

	fits_filename = samples_struct.fits_table[iframe];

	for (long idet=rank*det.ndet/size;idet<(rank+1)*det.ndet/size;idet++){
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


		//// Read pointing
		read_samptopix(ns, samptopix, tmp_dir, fits_filename,field);


		if (S != NULL){

			if (addnpix){
				deproject(S,indpix,samptopix,ns,NAXIS1, NAXIS2,npix,Ps,2,factdupl,samples_struct.ntotscan,indpsrc,npixsrc);
			} else {
				deproject(S,indpix,samptopix,ns,NAXIS1, NAXIS2,npix,Ps,2,factdupl);
			}

		}


		if (S != NULL){
			//********************  pre-processing of data ********************//
			MapMakPreProcessData(data,flag,ns,proc_param.napod,proc_param.poly_order,f_lppix,data_lp,
					proc_param.NORMLIN,proc_param.NOFILLGAP,proc_param.remove_polynomia,Ps);
		}
		else {
			MapMakPreProcessData(data,flag,ns,proc_param.napod,proc_param.poly_order,f_lppix,data_lp,
					proc_param.NORMLIN,proc_param.NOFILLGAP,proc_param.remove_polynomia);
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

		readNSpectrum(nameSpfile,bfilter,ns,proc_param.fsamp,Nk);


		//********************** compute P^t N-1 d ************************//
		compute_PtNmd(data_lp,Nk,ns,NAXIS1, NAXIS2,indpix,samptopix,npix,PNd);


	}// end of idet loop


	delete[] samptopix;
	delete[] bfilter;
	delete[] Nk;

	delete[] data;
	delete[] data_lp;

	delete[] flag;

	delete[] Ps;


}






void do_PtNPS_nocorr(double *S, string *extentnoiseSp_all, struct common dir,
		struct detectors det,double f_lppix,double fsamp, bool flgdupl, long ns,
		long long *indpix, long NAXIS1, long NAXIS2, long long npix,
		long iframe,std::string fname, double *PtNPmatS, double *Mp, long *hits, int rank, int size)

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

	for (long idet=rank*det.ndet/size;idet<(rank+1)*det.ndet/size;idet++){

		field = det.boloname[idet];

		read_samptopix(ns, samptopix, dir.tmp_dir, fname,field);


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













