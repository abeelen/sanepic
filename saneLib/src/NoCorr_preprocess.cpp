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

int do_PtNd_nocorr(double *PNd,string tmp_dir, struct param_sanePre proc_param, struct param_sanePos pos_param,
		struct samples samples_struct, std::vector<std::string> det, long ndet, double f_lppix, double f_lppix_Nk,
		long addnpix, long ns, long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2, long long npix,
		long long npixsrc, long iframe, double *S, int para_bolo_indice, int para_bolo_size)
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


	string fits_filename, dirfile_filename;

	fits_filename = samples_struct.fitsvect[iframe];
	dirfile_filename = samples_struct.basevect[iframe];

	for (long idet=para_bolo_indice*ndet/para_bolo_size;idet<(para_bolo_indice+1)*ndet/para_bolo_size;idet++){
		field = det[idet];

		if(read_data_from_dirfile(samples_struct.dirfile_pointer, dirfile_filename, field, data, ns))
			return 1;
		if(read_flag_from_dirfile(samples_struct.dirfile_pointer, dirfile_filename, field, flag, ns))
			return 1;

		//		long test_ns;
		//		read_signal_from_fits(fits_filename, field, data, test_ns);
		//		if (test_ns != ns) {
		//			cerr << "Read signal does not correspond to frame size : Check !!" << endl;
		//			return 1;
		//		}
		//
		//		read_flag_from_fits(fits_filename , field, flag, test_ns);
		//		if (test_ns != ns) {
		//			cerr << "Read flag does not correspond to frame size : Check !!" << endl;
		//			return 1;
		//		}


		//// Read pointing
		read_samptopix(samples_struct.dirfile_pointer, ns, samptopix, dirfile_filename,field);


		if (S != NULL){

			deproject(S,indpix,samptopix,ns,NAXIS1, NAXIS2,npix,Ps,2,factdupl,samples_struct.ntotscan,indpsrc,npixsrc);

			//********************  pre-processing of data ********************//
			MapMakePreProcessData(data,  flag, ns, proc_param, f_lppix, data_lp,  Ps);

		} else {
			MapMakePreProcessData(data,  flag, ns, proc_param, f_lppix, data_lp, NULL);
		}


		/************************************************************************************/

		// end of preprocess begin of fdata

		for (long ii=0;ii<ns/2+1;ii++){
			powered=pow(double(ii)/f_lppix, 16);
			bfilter[ii] = powered /(1.0+powered);
		}



		//****************** Compute (or read) input power spectrum of the NOISE  ***************//
		extentnoiseSp = samples_struct.noisevect[iframe];
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


	return 0;
}






void do_PtNPS_nocorr(struct samples samples_struct, double *S, std::vector<std::string> noisevect, struct param_common dir,
		std::vector<std::string> det, long ndet, double f_lppix,double fsamp, bool flgdupl, long ns,
		long long *indpix, long NAXIS1, long NAXIS2, long long npix,
		long iframe,std::string fname, double *PtNPmatS, double *Mp, long *hits, int para_bolo_indice, int para_bolo_size)

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

	for (long idet=para_bolo_indice*ndet/para_bolo_size;idet<(para_bolo_indice+1)*ndet/para_bolo_size;idet++){

		field = det[idet];

		read_samptopix(samples_struct.dirfile_pointer, ns, samptopix, fname,field);


		// AS
		deproject(S,indpix,samptopix,ns,NAXIS1, NAXIS2,npix,Ps,flgdupl,factdupl);

		for (long ii=0;ii<(ns)/2+1;ii++){
			powered=pow(double(ii)/f_lppix, 16);
			bfilter[ii] = powered /(1.0+powered);
		}

		extentnoiseSp = noisevect[iframe];
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
