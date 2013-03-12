#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

#include <ostream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>

#include "CorrPreprocess.h"
#include "DataIO.h"
#include "TemporaryIO.h"
#include "InputFileIO.h"
#include "TodProcess.h"
#include "MapMaking.h"
#include "StructDefinition.h"
#include "CovMatrixIO.h"

#include <time.h>

#include <gsl/gsl_math.h>



#include <fftw3.h>

using namespace std;

int write_tfAS(struct samples samples_struct, double *S, long long *indpix, long NAXIS1, long NAXIS2, long long npix,
		bool flgdupl, long iframe, int sub_rank, int sub_size) {

	double *Ps;
	long long *samptopix;

	fftw_plan fftplan;
	fftw_complex *fdata;

	std::vector<std::string> det = samples_struct.bolo_list[iframe];
	long ndet                    = (long) det.size();
	long ns				         = samples_struct.nsamples[iframe];

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

	for (long idet1 = floor(sub_rank*ndet*1.0/sub_size); idet1<floor((sub_rank+1)*ndet*1.0/sub_size); idet1++){

		//Read pointing data
		if(readSampToPix(samples_struct.dirfile_pointers[iframe], samples_struct.basevect[iframe] , det[idet1], samptopix,ns))
			return 1;

		deproject(S,indpix,samptopix,ns,NAXIS1, NAXIS2,npix,Ps,flgdupl,factdupl);

		//Fourier transform of the data
		fftw_execute_dft_r2c(fftplan, Ps, fdata);

		if(writeFdata(samples_struct.dirfile_pointers[iframe], ns, fdata, "fPs_", idet1, samples_struct.basevect[iframe], det))
			return 1;

	}
	fftw_destroy_plan(fftplan);

	delete[] samptopix;
	delete[] Ps;
	delete[] fdata;

	return 0;
}

int write_ftrProcesdata(double *S, struct param_saneProc proc_param, struct samples samples_struct, struct param_sanePos pos_param,
		string tmp_dir, long long *indpix, long long *indpsrc, long NAXIS1, long NAXIS2,
		long long npix,	long long npixsrc, long long addnpix, long iframe, int sub_rank, int sub_size) {


	double *data, *Ps=NULL;
	int *flag=NULL;
	long long *samptopix;

	double *bfilter;

	std::vector<std::string> det = samples_struct.bolo_list[iframe];
	long ndet                    = (long) det.size();
	long ns				         = samples_struct.nsamples[iframe];
	double fsamp                 = samples_struct.fsamp[iframe];
	double fhp_pix               = samples_struct.fhp[iframe]* double(ns) / fsamp;

	bfilter = new double[ns / 2 + 1];
	butterworth_filter(ns, fhp_pix, 8, bfilter);

	fftw_plan fftplan;
	fftw_complex *fdata;

	string field1, dirfile_filename;

	samptopix = new long long[ns];
	flag      = new int[ns];


	data      = (double *) fftw_malloc(sizeof(double)*ns);
	fdata     = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*(ns/2+1));

	fill(data,data+ns,0.0);
	fill(samptopix,samptopix+ns,0);

	int factdupl = 1;
	if(pos_param.flgdupl==1)		factdupl = 2;

	dirfile_filename = samples_struct.basevect[iframe];

	fftplan = fftw_plan_dft_r2c_1d(ns, data, fdata, FFTW_ESTIMATE);

	for (long idet1 = floor(sub_rank*ndet*1.0/sub_size); idet1<floor((sub_rank+1)*ndet*1.0/sub_size); idet1++){

		field1 = det[idet1];

		for (long ii=0;ii<ns/2+1;ii++){
			fdata[ii][0] = 0.0;
			fdata[ii][1] = 0.0;
		}

		if(readDataFromDirfile(samples_struct.dirfile_pointers[iframe], dirfile_filename, field1, data, ns))
			return 1;
		if(readFlagFromDirfile(samples_struct.dirfile_pointers[iframe], dirfile_filename, field1, flag, ns))
			return 1;

		if (S != NULL){
			//// Read pointing
			if(readSampToPix(samples_struct.dirfile_pointers[iframe],  dirfile_filename, field1, samptopix, ns))
				return 1;

			Ps        = new double[ns];
			fill(Ps,Ps+ns,0.0);


			deproject(S,indpix,samptopix,ns,NAXIS1, NAXIS2,npix,Ps,2,factdupl);

			//********************  pre-processing of data ********************//
			MapMakePreProcessData(data,  flag, ns, proc_param, bfilter, Ps);

			delete[] Ps;
		} else
			MapMakePreProcessData(data,  flag, ns, proc_param, bfilter, NULL);

		//
		//		if (field1=="PLWD8" && iframe == 0)
		//		{
		//			FILE *fp;
		//			if ((fp = fopen("PLWD8", "a+")) != NULL) {
		//				for (long iFile = 0; iFile < ns; iFile++)
		//					if (S != NULL){
		//						fprintf(fp,"%i %f %f\n",iFile, data_lp[iFile], Ps[iFile]);
		//					}else {
		//						fprintf(fp,"%i %f %f \n",iFile, data_lp[iFile], 0);
		//					}
		//						fclose(fp);
		//			} else {
		//				cerr << "ERROR : Could not open file " << endl;
		//				return 1;
		//			}
		//		}


		//Fourier transform of the data
		fftw_execute_dft_r2c(fftplan, data, fdata);

		//write fourier transform to disk
		if(writeFdata(samples_struct.dirfile_pointers[iframe], ns, fdata, "fData_", idet1, dirfile_filename, det))
			return 1;


	} // idet1

	fftw_destroy_plan(fftplan);

	delete [] flag;
	delete [] data;
	delete [] samptopix;
	delete [] fdata;
	delete [] bfilter;

	return 0;
}

int do_PtNd(struct samples samples_struct, double *PNd, string prefixe,
		int sub_rank, int sub_size,
		long long *indpix, long NAXIS1, long NAXIS2, long long npix, long iframe,
		double *Mp, long *hits)
{

	//	long  nbins;
	string field1, field2;
	string extentNoiseSp;

	string nameSpfile;

	long long *samptopix;
	double *ell, *SpN, *bfilter_InvSquared, *Nk, *Nd;
	double *km;

	fftw_plan fftplan;
	fftw_complex *fdata, *Ndf;

	std::vector<std::string> det = samples_struct.bolo_list[iframe];
	long ndet                    = (long) det.size();
	long ns				         = samples_struct.nsamples[iframe];
	double fsamp                 = samples_struct.fsamp[iframe];
	double fhp_pix               = samples_struct.fcut[iframe] * double(ns) / fsamp;

	samptopix = new long long[ns];
	bfilter_InvSquared = new double[ns/2+1];
	Nk = new double[ns/2+1];
	fdata = new fftw_complex[ns/2+1];

	Nd   = (double *) fftw_malloc(sizeof(double)*ns);
	Ndf  = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*(ns/2+1));

	double **SpN_all;

	butterworth_filter(ns, fhp_pix, 8, bfilter_InvSquared);
	for (long ii=0;ii<ns/2+1;ii++)
		bfilter_InvSquared[ii] = gsl_pow_2(1.0/(bfilter_InvSquared[ii]+0.000001));

	fftplan = fftw_plan_dft_c2r_1d(ns, Ndf, Nd, FFTW_ESTIMATE);

	for (long idet1 = floor(sub_rank*ndet*1.0/sub_size); idet1<floor((sub_rank+1)*ndet*1.0/sub_size); idet1++){

		field1 = det[idet1];

		//Read pointing data
		if(readSampToPix(samples_struct.dirfile_pointers[iframe], samples_struct.basevect[iframe], field1, samptopix, ns))
			return 1;

		//**************************************** Noise power spectrum

		//read noise PS file for idet1
		long ndet2 = samples_struct.ndet[iframe];
		long nbins = samples_struct.nbins[iframe];

		if(read_InvNoisePowerSpectra(samples_struct.dirfile_pointers[iframe], field1,  samples_struct.basevect[iframe], nbins, ndet2, &ell, &SpN_all))
			return 1;

		km = new double[nbins];
		//	logSpN = new double[nbins];

		for (long ii=0;ii<nbins;ii++)
			km[ii] = exp((log(ell[ii+1])+log(ell[ii]))/2.0)*ns/fsamp;

		delete[] ell;

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
			if(readFdata(samples_struct.dirfile_pointers[iframe],  samples_struct.basevect[iframe], det[idet2], prefixe, fdata, ns))
				return 1;

			//****************** Cross power spectrum of the noise  ***************//
			for (int ii=0;ii<nbins;ii++)
				SpN[ii] = SpN_all[idet2][ii] / (double) ns;


			// interpolate logarithmically the inversed noise power spectrum
			InvbinnedSpectrum2log_interpol(km,SpN,bfilter_InvSquared,nbins,ns,Nk, NULL);

			for (long jj=0;jj<ns/2+1;jj++){
				if (isnan(Nk[jj])) {
					printf("A NaN has been found in Nk : iframe %ld, det1 %s, det2 %s\n",iframe, det[idet1].c_str(), det[idet2].c_str());
					return 1;
				}
			}

			//********************************* compute N^-1 d  ***********************//
			for (long ii=0;ii<ns/2+1;ii++){
				Ndf[ii][0] += (fdata[ii][0]*Nk[ii]);
				Ndf[ii][1] += (fdata[ii][1]*Nk[ii]);

			}

			//Compute weight map for preconditioner
			if ((Mp != NULL) && (idet2 == idet1))
				compute_diagPtNPCorr(Nk,samptopix,ns,NAXIS1, NAXIS2,indpix,npix,fhp_pix,Mp);

		}// end of idet2 loop

		delete [] km;
		delete[] SpN;
		free_dmatrix(SpN_all,0,ndet-1,0,nbins-1);

		fftw_execute_dft_c2r(fftplan, Ndf, Nd);

		for (long ii=0;ii<ns;ii++){
			PNd[indpix[samptopix[ii]]] += Nd[ii]; // Nd real
		}

		//compute hit counts
		if (hits != NULL){
			for (long ii=0;ii<ns;ii++){
				hits[indpix[samptopix[ii]]] += 1;
			}
		}




	}// end of idet1 loop

	fftw_destroy_plan(fftplan);

	delete[] Nk;
	delete[] bfilter_InvSquared;
	delete[] samptopix;
	delete[] Nd;
	delete[] fdata;
	delete[] Ndf;

	return 0;
}
