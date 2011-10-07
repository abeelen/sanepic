#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "todprocess.h"
#include "map_making.h"

#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort_double.h>

#define NR_END 1
#define FREE_ARG char*

using namespace std;


long long *data_compare;



int compare_long_long (const void *a, const void *b)
{
	const long long *da = (const long long *) a;
	const long long *db = (const long long *) b;

	return (*da > *db) - (*da < *db);
}

int compare_global_array_long_long (const void *array_1, const void *array_2)
{

	const long long *long_array_1 = (const long long *) array_1;
	const long long *long_array_2 = (const long long *) array_2;

	return (data_compare[*long_array_1] > data_compare[*long_array_2]) - (data_compare[*long_array_1] < data_compare[*long_array_2]);
}

void compute_PtNmd(double *data, double *Nk, long ndata, long NAXIS1, long NAXIS2,
		long long *indpix, long long *samptopix, long long npix, double *PNd){


	double *Nd;
	fftw_complex  *fdata, *Ndf;
	fftw_plan fftplan;


	fdata = new fftw_complex[ndata/2+1];
	Ndf = new fftw_complex[ndata/2+1];
	Nd = new double[ndata];


	//Fourier transform of the data
	fftplan = fftw_plan_dft_r2c_1d(ndata, data, fdata, FFTW_ESTIMATE);
	fftw_execute(fftplan);
	fftw_destroy_plan(fftplan);


	for (long k=0;k<ndata/2+1;k++){
		Ndf[k][0] = fdata[k][0]/Nk[k]/(double)ndata/(double)ndata;
		Ndf[k][1] = fdata[k][1]/Nk[k]/(double)ndata/(double)ndata;
	}


	fftplan = fftw_plan_dft_c2r_1d(ndata, Ndf, Nd, FFTW_ESTIMATE);
	fftw_execute(fftplan);
	fftw_destroy_plan(fftplan);

	for (long ii=0;ii<ndata;ii++){
		if ((ii < 0) || (ii >= ndata)){
			PNd[npix-2] += Nd[ii];
		} else {
			PNd[indpix[samptopix[ii]]] += Nd[ii];
		}
	}



	delete [] Ndf;
	delete [] fdata;
	delete [] Nd;


}

void compute_diagPtNP(double *Nk, long long *samptopix, long ndata,
		long  NAXIS1, long NAXIS2, long long *indpix,
		long npix, double f_lppix, double *dPtNP){


	long long kk2, ii2, ndataf;
	long long *pixpos;
	long long count, count_;
	long long *pixtosamp;

	//fft stuff
	fftw_complex  *Nk_;
	double *N_;
	fftw_plan fftplan;

	Nk_ = new fftw_complex[ndata/2+1];
	N_ = new double[ndata];
	pixpos = new long long[ndata];
	data_compare = new long long[ndata];
	pixtosamp = new long long[ndata];


	// N^-1
	for (long k=0;k<ndata/2+1;k++){
		Nk_[k][0] = 1.0/abs(Nk[k])/(double)ndata/(double)ndata;
		Nk_[k][1] = 0.0;
	}
	fftplan = fftw_plan_dft_c2r_1d(ndata, Nk_, N_, FFTW_ESTIMATE);
	fftw_execute(fftplan);

	for (long ii=0;ii<ndata;ii++){
		if ((ii < 0) || (ii >= ndata)){
			pixpos[ii] = npix-2;
		} else {
			pixpos[ii] = indpix[samptopix[ii]];
		}
	}

	qsort(pixtosamp,ndata,sizeof(long long),compare_global_array_long_long);
	qsort(data_compare,ndata,sizeof(long long),compare_long_long);

	//	gsl_sort_index ((size_t*)pixtosamp, (double*)pixpos, 1, ndata); // replace qsort with compare_global_array_long_long

	//	for(long ii=0; ii< ndata; ii++)
	//		data_compare[ii]=pixpos[pixtosamp[ii]];

	ndataf = (ndata)/MAX(2,int(f_lppix+0.5));

	count = 0;

	for (long ipix=data_compare[0];ipix<npix;ipix++){

		count_ = count;

		while((count < ndata) && (data_compare[count] == ipix))
			count++;

		if (count-count_ > 0){
			for (long long ii=count_;ii<count;ii++){
				ii2 = pixtosamp[ii];
				if ((ipix == npix-1) || (ipix == npix-2)){ //This is just to avoid spending to much time computing this pixel which could contain a lot of data
					dPtNP[ipix] += N_[0];
					//printf("TEST");
				} else {
					for (long long kk=count_;kk<count;kk++){
						kk2 = pixtosamp[kk];
						if (abs(kk2-ii2) < ndataf)
							dPtNP[ipix] += N_[abs(ii2-kk2)];
					}
				}
			}
		}
	}

	// memory dealloc
	delete[] N_;
	delete[] Nk_;
	delete[] pixpos;
	delete[] pixtosamp;
	delete[] data_compare;


	//clean up
	fftw_destroy_plan(fftplan);

}


void compute_diagPtNPCorr(double *Nk, long long *samptopix, long ndata,
		long NAXIS1, long NAXIS2, long long *indpix,
		long long npix, double f_lppix, double *dPtNP){


	long long kk2, ii2, ndataf;
	//	long long *pixpos;
	long long count, count_;
	long long *pixtosamp;
	//	size_t *pixtosamp;

	//fft stuff
	fftw_complex  *Nk_;
	double *N_;
	fftw_plan fftplan;

	//	time_t rawtime;
	//	struct tm * timeinfo;
	//	time ( &rawtime );
	//	timeinfo = localtime ( &rawtime );
	//	file << "do_ptnd : " << oss.str() << " at " << asctime (timeinfo) << endl;

	Nk_ = new fftw_complex[ndata/2+1];
	N_ = new double[ndata];
	//	pixpos = new long long [ndata];
	data_compare = new long long[ndata];
	pixtosamp = new long long[ndata];
	//	pixtosamp = new size_t[ndata];

	//fill(N_,N_+ndata,0.0);
	//fill(pixpos,pixpos+ndata,0);

	// N^-1
	for (long k=0;k<ndata/2+1;k++){
		Nk_[k][0] = abs(Nk[k])/sqrt((double)ndata);
		Nk_[k][1] = 0.0;
	}
	fftplan = fftw_plan_dft_c2r_1d(ndata, Nk_, N_, FFTW_ESTIMATE);
	fftw_execute(fftplan);


	for (long ii=0;ii<ndata;ii++)
		data_compare[ii] = indpix[samptopix[ii]];

	for (long ii=0;ii<ndata;ii++)
		pixtosamp[ii] = ii;

	//	for (long ii=0;ii<ndata;ii++){
	//		data_compare[ii] = pixpos[ii];
	//	}

	qsort(pixtosamp,ndata,sizeof(long long),compare_global_array_long_long);
	qsort(data_compare,ndata,sizeof(long long),compare_long_long);

	//	gsl_sort_index ((size_t*)pixtosamp, (double*)pixpos, 1, ndata); // replace qsort with compare_global_array_long_long

	//	for(long ii=0; ii< ndata; ii++)
	//		data_compare[ii]=pixpos[pixtosamp[ii]];

	ndataf = (ndata)/MAX(2,int(f_lppix+0.5));

	count = 0;

	for (long long ipix=data_compare[0];ipix<npix;ipix++){

		count_ = count;

		while((count < ndata) && (data_compare[count] == ipix))
			count++;

		if (count-count_ > 0){
			for (long ii=count_;ii<count;ii++){
				ii2 = pixtosamp[ii];
				if ((ipix == npix-1) || (ipix == npix-2)){ //This is just to avoid spending to much time computing this pixel
					dPtNP[ipix] += N_[0];
				} else {
					for (long long kk=count_;kk<count;kk++){
						kk2 = pixtosamp[kk];
						if (abs(kk2-ii2) < ndataf)
							dPtNP[ipix] += N_[abs(ii2-kk2)];
					}
				}
			}
		}
	}


	delete[] N_;
	delete[] Nk_;
	//	delete[] pixpos;
	delete[] pixtosamp;
	delete[] data_compare;


	//clean up
	fftw_destroy_plan(fftplan);


}

void MapMakePreProcessData(double *data,  int *flag, long ns, struct param_saneProc proc_param, double f_lppix, double *data_lp,  double *Ps){


	double aa, bb;

	double *data_out;
	double *bfilter;

	data_out = new double[ns];
	bfilter = new double[ns/2+1];

	fill(data_out,data_out+ns,0.0);
	fill(bfilter,bfilter+ns/2+1,0.0);


	//TODO : TEST : Change the removal of the map here => means ??
	// Optimized the memory management here....
	// This routine CHANGES *data : this is really wanted ! => YES

	if (Ps != NULL)
		for (long ii=0;ii<ns;ii++)
			data[ii] = data[ii] - Ps[ii];


	//*********************************************************************
	if (proc_param.fill_gap)
		fillgaps2(data,ns,data_out,flag,40);


	if(proc_param.remove_polynomia)
		//remove polynomia to correct from time varying calibration
		remove_poly(data_out,ns,proc_param.poly_order,data_lp,flag);

	//linear prediction
	for (long ii=0;ii<ns;ii++)
		data[ii] = data_lp[ii];


	if (proc_param.remove_linear){
		/// remove a baseline
		aa = (data_lp[ns-1]-data[0])/double(ns);
		bb = data_lp[0];
		for (long ii=0;ii<ns;ii++)
			data_lp[ii] -= aa*(double)ii+bb;
	}

	//Butterworth filter (if necessary)
	if (proc_param.highpass_filter){
		butterworth(data_lp,ns,f_lppix,8,data_out,bfilter,1,proc_param.napod,0);
	} /*else{
		for (long ii=0;ii<(ns)/2+1;ii++)
			bfilter[ii] = 1.0;
	}*/



	//******************* process gaps
	if (proc_param.fill_gap)
		fillgaps2(data_out,ns,data_lp,flag,40);

	if (Ps != NULL)
		for (long ii=0;ii<ns;ii++)
			data_lp[ii] = data_lp[ii] + Ps[ii];



	delete [] data_out;
	delete [] bfilter;


}


///// measure power spectrum of the uncorrelated part of the noise
void noisepectrum_estim(double *data, long ns, double *ell, int nbins, double fsamp, double *bfilter, double *Nell, double *Nk){

	//TODO : should be almost the same as noisecrossspectrum_estim : i.e. take fdata as input
	int qq;
	double totapod;

	double *datatemp, *datatemp2, *apodwind, *bfiltertemp;
	int *count;

	fftw_complex  *fdata;
	fftw_plan fftplan;

	datatemp = new double[ns];
	datatemp2 = new double[ns];
	bfiltertemp = new double[ns/2+1];
	count = new int[nbins];
	fdata = new fftw_complex[ns/2+1];

	//TODO : Why removing a polynomial here ???? SHOULD NOT
	//TODO : apodization done twice ?? SHOULD NOT, but see below

	remove_poly(data,ns,4,datatemp2,0);
	apodwind = apodwindow(ns,ns/10);
	for (long ii=0;ii<ns;ii++)
		datatemp[ii] = datatemp2[ii]*apodwind[ii];


	//Fourier transform the data
	fftplan = fftw_plan_dft_r2c_1d(ns, datatemp, fdata, FFTW_ESTIMATE);
	fftw_execute(fftplan);
	fftw_destroy_plan(fftplan);


	// TODO: apodisation done twice ??
	// TODO: resamblance to factapod in other routines, can we merge that ? (factor ns)
	totapod = 0.0;
	for (long ii=0;ii<ns;ii++)
		totapod += apodwind[ii]*apodwind[ii];


	//power spectrum
	for (long k=0;k<ns/2+1;k++){
		Nk[k] = gsl_pow_2(fdata[k][0]) + gsl_pow_2(fdata[k][1]);
		Nk[k] = Nk[k]/(totapod/(double)ns)/(double)ns;
	}


	//bin power spectrum
	for (int q=0;q<nbins;q++){
		Nell[q] = 0.0;
		count[q] = 0;
	}

	qq=0;
	for (long k=0;k<ns/2+1;k++){
		if (k >= ell[qq+1]/fsamp*(double)ns)
			if (qq<nbins-1)
				qq++;
		Nell[qq] += Nk[k];
		count[qq] += 1;
	}

	for (int q=0;q<nbins;q++)
		Nell[q] /= double(count[q]);

	for(long ii=0;ii<ns/2+1;ii++){
		if (bfilter == NULL){
			bfiltertemp[ii] = 1.0;
		} else {
			bfiltertemp[ii] = bfilter[ii];
		}
	}


	// interpolate logarithmically the spectrum and filter
	binnedSpectrum2log_interpol(ell,Nell,bfiltertemp,nbins,ns,fsamp,Nk,NULL);


	//clean up
	delete [] count;
	delete [] datatemp;
	delete [] datatemp2;
	delete [] fdata;
	delete [] bfiltertemp;
	delete []  apodwind;


}



void noisecrosspectrum_estim(fftw_complex *fdata1, fftw_complex *fdata2, int ns, double *ell, int nbins, double fsamp, double *bfilter, double *Nell, double *Nk){

	// TODO : no apodization factor correction here ??
	// TODO : should be almost the same as noisespectrum_estim
	int qq;

	double *bfiltertemp;
	int *count;


	bfiltertemp = new double[ns/2+1];
	count = new int[nbins];


	//power spectrum
	for (long k=0;k<ns/2+1;k++){
		Nk[k] = fdata1[k][0]*fdata2[k][0] + fdata1[k][1]*fdata2[k][1];
		Nk[k] = Nk[k]/(double)ns;
	}


	//bin power spectrum
	for (int q=0;q<nbins;q++){
		Nell[q] = 0.0;
		count[q] = 0;
	}


	qq=0;
	for (long k=0;k<ns/2+1;k++){
		if (k >= ell[qq+1]/fsamp*(double)ns)
			if (qq<nbins-1)
				qq++;
		Nell[qq] += Nk[k];
		count[qq] += 1;
	}


	for (int q=0;q<nbins;q++)
		Nell[q] /= double(count[q]);



	for(long ii=0;ii<ns/2+1;ii++){
		if (bfilter == NULL){
			bfiltertemp[ii] = 1.0;
		} else {
			bfiltertemp[ii] = bfilter[ii];
		}
	}


	// interpol logarithmically the spectrum and filter
	binnedSpectrum2log_interpol(ell,Nell,bfiltertemp,nbins,ns,fsamp,Nk,NULL);



	//clean up
	delete [] count;
	delete [] bfiltertemp;

}

int readNSpectrum(string nameSpfile, double *bfilter, long ns, double fsamp, double *Nk){

	FILE *fp;

	int nbins;
	//long ii;
	double dummy1, dummy2;

	double *SpN;
	double *ell;
	int result;


	if ((fp = fopen(nameSpfile.c_str(),"r")) == NULL){
		cerr << "ERROR: Can't find noise power spectra file, check -k or -K in command line. Exiting. \n";
		return 1;
	}
	result = fscanf(fp,"%d",&nbins);


	SpN = new double[nbins];
	ell = new double[nbins+1];

	for (int ii=0;ii<nbins;ii++){
		result = fscanf(fp,"%lf %lf",&dummy1,&dummy2);
		ell[ii] = dummy1;
		SpN[ii] = dummy2;
	}
	result = fscanf(fp,"%lf",&dummy1);
	ell[nbins] = dummy1;
	fclose(fp);


	// interpolate logarithmically the noise power spectrum
	binnedSpectrum2log_interpol(ell,SpN,bfilter,nbins,ns,fsamp,Nk);



	delete[] SpN;
	delete[] ell;

	return 0;

}


void deproject(double *S, long long *indpix, long long *samptopix, long long ndata, long NAXIS1, long NAXIS2, long long npix, double *Ps, int flgdupl, int factdupl, long ntotscan, long long *indpsrc, long long npixsrc){


	//	double a, b;
	//long ll;

	for (long long ii=0;ii<ndata;ii++){

		Ps[ii] = S[indpix[samptopix[ii]]];

		// in case we replaced flagged data
		if ((flgdupl == 2) && (samptopix[ii] >= NAXIS1*NAXIS2) && (samptopix[ii] < 2*NAXIS1*NAXIS2) && (indpix[samptopix[ii] - NAXIS1*NAXIS2] >= 0)){
			if (indpix[samptopix[ii] - NAXIS1*NAXIS2] >= 0){
				Ps[ii] = S[indpix[samptopix[ii]-NAXIS1*NAXIS2]];
			}
			//			else {
			//					cout << "there " << endl;
			//				a = 0.0;
			//				b = 0.0;
			//				if (ntotscan){
			//					for (long iframe=0;iframe<ntotscan;iframe++){
			//						if (indpix[factdupl*NAXIS1*NAXIS2 + indpsrc[samptopix[ii] - NAXIS1*NAXIS2] + iframe*npixsrc] >= 0){
			//							a += S[indpix[factdupl*NAXIS1*NAXIS2 + indpsrc[samptopix[ii] - NAXIS1*NAXIS2] + iframe*npixsrc]];
			//							b++;
			//						}
			//					}
			//				}
			//				if (b > 0.5)
			//					Ps[ii] = a/b;
			//			}
		}

	}

}
