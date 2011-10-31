#include <ostream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>

#include "Corr_preprocess.h"
#include "dataIO.h"
#include "temporary_IO.h"
#include "inputFileIO.h"
#include "todprocess.h"
#include "map_making.h"
#include "struct_definition.h"
#include "covMatrix_IO.h"

#include <time.h>

#include <gsl/gsl_math.h>
#include <fftw3.h>

using namespace std;

#ifdef DEBUG
#include <ostream>
#include <sstream>
#endif

int write_tfAS(struct samples samples_struct, double *S, std::vector<std::string> det, long ndet,long long *indpix, long NAXIS1, long NAXIS2, long long npix,
		bool flgdupl, long ns, string filename, int para_bolo_indice, int para_bolo_size)
{



	double *Ps;
	long long *samptopix;

	fftw_plan fftplan;
	fftw_complex *fdata;

	samptopix = new long long[ns];

	Ps    = (double *) fftw_malloc(sizeof(double)*ns);
	fdata = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (ns/2+1));

	fill(samptopix,samptopix+ns,0);
	fill(Ps,Ps+ns,0.0);
	for(long ii =0; ii<ns/2+1;ii++){
		fdata[ii][0]=0.0;
		fdata[ii][1]=0.0;
	}

	int factdupl = 1;
	if(flgdupl==1)  factdupl = 2;

	fftplan = fftw_plan_dft_r2c_1d(ns, Ps, fdata, FFTW_ESTIMATE);

	for (long idet1=para_bolo_indice*ndet/para_bolo_size;idet1<(para_bolo_indice+1)*ndet/para_bolo_size;idet1++){

		//Read pointing data
		if(read_samptopix(samples_struct.dirfile_pointer, filename, det[idet1], samptopix,ns))
			return 1;

		deproject(S,indpix,samptopix,ns,NAXIS1, NAXIS2,npix,Ps,flgdupl,factdupl);

		//Fourier transform of the data
		fftw_execute_dft_r2c(fftplan, Ps, fdata);

		if(write_fdata(samples_struct.dirfile_pointer, ns, fdata, "fPs_", idet1, filename, det))
			return 1;

	}
	fftw_destroy_plan(fftplan);

	delete[] samptopix;
	delete[] Ps;
	delete[] fdata;

	return 0;
}

int write_ftrProcesdata(double *S, struct param_saneProc proc_param, struct samples samples_struct, struct param_sanePos pos_param,
		string tmp_dir,std::vector<std::string> det,long ndet, long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2,
		long long npix,	long long npixsrc, long long addnpix, double fhp_pix, long ns, long iframe, int para_bolo_indice, int para_bolo_size, std::string fname)
{



	double *data, *data_lp, *Ps=NULL;
	int *flag=NULL;
	long long *samptopix;

	fftw_plan fftplan;
	fftw_complex *fdata;

	string field1, dirfile_filename;

#ifdef DEBUG
	ofstream file;
	file.open(fname.c_str(), ios::out | ios::app);
	if(!file.is_open()){
		cerr << "File [" << file << "] Invalid." << endl;
		return 1;
	}
#endif

	data      = new double[ns];
	samptopix = new long long[ns];
	flag      = new int[ns];


	data_lp   = (double *) fftw_malloc(sizeof(double)*ns);
	fdata     = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*(ns/2+1));

	fill(data_lp,data_lp+ns,0.0);
	fill(samptopix,samptopix+ns,0);

	int factdupl = 1;
	if(pos_param.flgdupl==1)		factdupl = 2;

	dirfile_filename = samples_struct.basevect[iframe];

	fftplan = fftw_plan_dft_r2c_1d(ns, data_lp, fdata, FFTW_ESTIMATE);


	for (long idet1=para_bolo_indice*ndet/para_bolo_size;idet1<(para_bolo_indice+1)*ndet/para_bolo_size;idet1++){

#ifdef DEBUG
		cout << "[ " << para_bolo_indice << " ] progression write_ftr : " << 100.0*(1.0-((double)(para_bolo_indice+1)-(double)idet1*(double)para_bolo_size/(double)ndet)) << " %" << endl;
		ostringstream oss;
		oss << tmp_dir + "fdata_" << iframe << "_" << det[idet1] << ".bi";
		time_t rawtime;
		struct tm * timeinfo;
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		file << "Writing file : " << oss.str() << " at " << asctime (timeinfo) << endl;
#endif

		field1 = det[idet1];

		fill(data_lp,data_lp+ns,0.0);
			for (long ii=0;ii<ns/2+1;ii++){
					fdata[ii][0] = 0.0;
					fdata[ii][1] = 0.0;
				}

		if(read_data_from_dirfile(samples_struct.dirfile_pointer, dirfile_filename, field1, data, ns))
			return 1;
		if(read_flag_from_dirfile(samples_struct.dirfile_pointer, dirfile_filename, field1, flag, ns))
			return 1;

		if (S != NULL){
			//// Read pointing
			if(read_samptopix(samples_struct.dirfile_pointer,  dirfile_filename, field1, samptopix, ns))
				return 1;

			Ps        = new double[ns];
			fill(Ps,Ps+ns,0.0);


			deproject(S,indpix,samptopix,ns,NAXIS1, NAXIS2,npix,Ps,2,factdupl);

			//TODO : Ps should not be here...  remove the signal before or make the deproject inside MapMakePreProcess
			//TODO : write fdata inside MapMakePreProcess.. or create a function same is true in sanePS
			//********************  pre-processing of data ********************//
			MapMakePreProcessData(data,  flag, ns, proc_param, fhp_pix, data_lp, Ps);

			delete[] Ps;
		}else
			MapMakePreProcessData(data,  flag, ns, proc_param, fhp_pix, data_lp, NULL);


		//Fourier transform of the data
		fftw_execute_dft_r2c(fftplan, data_lp, fdata);

		//write fourier transform to disk
		if(write_fdata(samples_struct.dirfile_pointer, ns, fdata, "fdata_", idet1, dirfile_filename, det))
			return 1;


	} // idet1

	fftw_destroy_plan(fftplan);

	delete [] flag;
	delete [] data;
	delete [] data_lp;
	delete [] samptopix;
	delete [] fdata;


#ifdef DEBUG
	file.close();
#endif

	return 0;
}

int do_PtNd(struct samples samples_struct, double *PNd, string prefixe,
		std::vector<std::string> det, long ndet, double fhp_pix, double fsamp, long ns, int para_bolo_indice, int para_bolo_size,
		long long *indpix, long NAXIS1, long NAXIS2, long long npix, long iframe,
		double *Mp, long *hits,std::string fname)
{

	//	long  nbins;
	string field1, field2;
	string extentNoiseSp;

	string nameSpfile;

	long long *samptopix;
	double *ell, *SpN, *bfilter, *bfilter_, *Nk, *Nd;

	double powered;

	fftw_plan fftplan;
	fftw_complex *fdata, *Ndf;

#ifdef DEBUG
	ofstream file;
	file.open(fname.c_str(), ios::out | ios::app);
	if(!file.is_open()){
		cerr << "File [" << file << "] Invalid." << endl;
		return 1;
	}
#endif


	samptopix = new long long[ns];
	bfilter = new double[ns/2+1];
	bfilter_ = new double[ns/2+1];
	Nk = new double[ns/2+1];
	fdata = new fftw_complex[ns/2+1];

	Nd   = (double *) fftw_malloc(sizeof(double)*ns);
	Ndf  = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*(ns/2+1));

	double **SpN_all;

	// This is a butterworth filter.... why not use butterworth()
	// Cause we want 1/butterworth() + we don't want to deal with fourier transform here !
	for (long ii=0;ii<ns/2+1;ii++){
		powered=gsl_pow_int(double(ii)/fhp_pix,16);
		bfilter[ii] = powered /(1.0+powered);
	}
	for (long ii=0;ii<ns/2+1;ii++)
		bfilter_[ii] = 1.0/(bfilter[ii]+0.000001);

	fill(Nd,Nd+ns,0.0);
	fill(Nk,Nk+(ns/2+1),0.0);
	fill(samptopix,samptopix+ns,0);


	fftplan = fftw_plan_dft_c2r_1d(ns, Ndf, Nd, FFTW_ESTIMATE);


	for (long idet1=para_bolo_indice*ndet/para_bolo_size;idet1<(para_bolo_indice+1)*ndet/para_bolo_size;idet1++){



#ifdef DEBUG
		cout << "[ " << para_bolo_indice << " ] progression do_ptNd : " << 100.0*(1.0-((double)(para_bolo_indice+1)-(double)idet1*(double)para_bolo_size/(double)ndet)) << " %" << endl;
		ostringstream oss;
		oss << "frame : " << iframe << " bolo : " << det[idet1];
		time_t rawtime;
		struct tm * timeinfo;
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		file << "do_ptnd : " << oss.str() << " at " << asctime (timeinfo) << endl;
#endif

		field1 = det[idet1];

#ifdef DEBUG
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		file << "before read samptopix : at " << asctime (timeinfo) << endl;
#endif

		//Read pointing data
		if(read_samptopix(samples_struct.dirfile_pointer, samples_struct.basevect[iframe], field1, samptopix, ns))
			return 1;
#ifdef DEBUG
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		file << "After read samptopix : at " << asctime (timeinfo) << endl;
#endif

		//**************************************** Noise power spectrum
		string extname = "_InvNoisePS";
		string suffix = samples_struct.basevect[iframe] + extname;

		//read noise PS file for idet1
		long ndet2 = samples_struct.ndet[iframe];
		long nbins = samples_struct.nbins[iframe];

#ifdef DEBUG
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		file << "before read Invnoise : at " << asctime (timeinfo) << endl;
#endif

		if(read_InvNoisePowerSpectra(samples_struct.dirfile_pointer, field1,  suffix, nbins, ndet2, &ell, &SpN_all))
			return 1;

#ifdef DEBUG
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		file << "after read Invnoise : at " << asctime (timeinfo) << endl;
#endif
		//		if(ndet!=ndet2) cout << "Error. The number of detector in noisePower Spectra file must be egal to input bolofile number\n";


		SpN = new double[nbins];
		fill(SpN,SpN+nbins,0.0);

		//Init N-1d
		for (long ii=0;ii<ns/2+1;ii++){
			Ndf[ii][0] = 0.0;
			Ndf[ii][1] = 0.0;
		}


		for (long idet2=0;idet2<ndet;idet2++){
			field2 = det[idet2];

			fill(Nd,Nd+ns,0.0);
			fill(Nk,Nk+(ns/2+1),0.0);

			//read Fourier transform of the data
			if(read_fdata(samples_struct.dirfile_pointer,  samples_struct.basevect[iframe], det[idet2], prefixe, fdata, ns))
				return 1;

			//****************** Cross power spectrum of the noise  ***************//
			for (int ii=0;ii<nbins;ii++)
				SpN[ii] = SpN_all[idet2][ii];


#ifdef DEBUG
			time ( &rawtime );
			timeinfo = localtime ( &rawtime );
			file << "before Invbinned : at " << asctime (timeinfo) << endl;
#endif
			// TODO : Why do we need to reinterpolate the noise power spectrum here ?
			// interpolate logarithmically the noise power spectrum
			InvbinnedSpectrum2log_interpol(ell,SpN,bfilter_,nbins,ns,fsamp,Nk, NULL);

#ifdef DEBUG
			time ( &rawtime );
			timeinfo = localtime ( &rawtime );
			file << "after Invbinned : at " << asctime (timeinfo) << endl;
#endif


			for (long jj=0;jj<ns/2+1;jj++){
				if (isnan(Nk[jj])) {
					printf("A NaN has been found in Nk : iframe %ld, det1 %ld, det2 %ld\n",iframe, idet1, idet2);
					return 1;
				}
			}

			//********************************* compute N^-1 d  ***********************//
			for (long ii=0;ii<ns/2+1;ii++){
				Ndf[ii][0] += (fdata[ii][0]*Nk[ii]);
				Ndf[ii][1] += (fdata[ii][1]*Nk[ii]);

			}



#ifdef DEBUG
			time ( &rawtime );
			timeinfo = localtime ( &rawtime );
			file << "before compute diag : at " << asctime (timeinfo) << endl;
#endif
			//Compute weight map for preconditioner
			if ((Mp != NULL) && (idet2 == idet1))
				compute_diagPtNPCorr(Nk,samptopix,ns,NAXIS1, NAXIS2,indpix,npix,fhp_pix,Mp);
#ifdef DEBUG
			time ( &rawtime );
			timeinfo = localtime ( &rawtime );
			file << "after compute diag : at " << asctime (timeinfo) << endl;
#endif


		}// end of idet2 loop

		// dEBUG

#ifdef DEBUG
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		file << "before fft : at " << asctime (timeinfo) << endl;
#endif
		fftw_execute_dft_c2r(fftplan, Ndf, Nd);

#ifdef DEBUG
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		file << "after fft : at " << asctime (timeinfo) << endl;
#endif

		for (long ii=0;ii<ns;ii++){
			PNd[indpix[samptopix[ii]]] += Nd[ii]; // Nd real
		}

		//compute hit counts
		if (hits != NULL){
			for (long ii=0;ii<ns;ii++){
				hits[indpix[samptopix[ii]]] += 1;
			}
		}


		delete[] ell;
		delete[] SpN;
		free_dmatrix(SpN_all,0,ndet-1,0,nbins-1);


	}// end of idet1 loop

	fftw_destroy_plan(fftplan);

	delete[] samptopix;
	delete[] Nd;
	delete[] bfilter;
	delete[] bfilter_;
	delete[] Nk;
	delete[] fdata;
	delete[] Ndf;


#ifdef DEBUG
	file.close();
#endif

	return 0;
}
