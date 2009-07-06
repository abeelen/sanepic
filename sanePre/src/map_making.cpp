#include <iostream>

#include <cmath>
//#include "/fir_data/patanch/numrec/inc/nrutil.h"
#include "todprocess.h"
#include "map_making.h"
#include <fftw3.h>
#include <time.h>

#define NR_END 1
#define FREE_ARG char*

using namespace std;


long *data_compare;


void compute_PtNmd(double *data, double *Nk, long ndata, int nn,
		long *indpix, long *samptopix, int npix, double *PNd){

	//long ii, k;
	//long ll;

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



	//for (ii=-marge;ii<ndata-marge;ii++){
	// if ((ii < 0) || (ii >= ndata-2*marge)){
	//   PNd[indpix[factdupl*nn*nn]] += Nd[ii+marge];
	// } else {
	//   if (rejectsamp[ii] == 2)
	//	PNd[indpix[factdupl*nn*nn]] += Nd[ii+marge];
	//   if (rejectsamp[ii] == 0){
	//	ll = indpix[yy[ii]*nn + xx[ii]];
	//	PNd[ll] += Nd[ii+marge];
	//   }
	//   if (rejectsamp[ii] == 1){
	//	if (flgdupl){
	//	  ll = indpix[(yy[ii]*nn + xx[ii])+nn*nn];
	//	  PNd[ll] += Nd[ii+marge];
	//	} else {
	//	  PNd[indpix[nn*nn+1]] += Nd[ii+marge];
	//	}
	//	if (rejectsamp[ii] == 3){
	//	  PNd[indpix[factdupl*nn*nn+1]] += Nd[ii+marge];
	//	}
	//   }
	// }
	//}




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





/*

void compute_PtNmd_corr(double *data, double *Nk, unsigned char *rejectsamp, unsigned char *binsamp,
		long ndata, int *xx, int *yy, int nn,
		long *indpix, int npix, double *PNd){

	long ii, k, ll;

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


	for (k=0;k<ndata/2+1;k++){
		Ndf[k][0] = fdata[k][0]*Nk[k];
		Ndf[k][1] = fdata[k][1]*Nk[k];
	}

	fftplan = fftw_plan_dft_c2r_1d(ndata, Ndf, Nd, FFTW_ESTIMATE);
	fftw_execute(fftplan);
	fftw_destroy_plan(fftplan);



	for (ii=0;ii<ndata;ii++){
		if ((ii < 0) || (ii >= ndata)){
			PNd[npix-2] += Nd[ii];
		} else {
			if (rejectsamp[ii] == 0){
				if (binsamp[ii] == 1){
					PNd[npix-2] += Nd[ii];
				} else {
					ll = indpix[yy[ii]*nn + xx[ii]];
					PNd[ll] += Nd[ii];
				}
			} else {
				PNd[npix-1] += Nd[ii];
			}
		}
	}

	delete [] Ndf;
	delete [] fdata;
	delete [] Nd;


}

 */



/*
void compute_PtNmfftd_corr(fftw_complex *fdata, double *Nk, unsigned char *rejectsamp, unsigned char *binsamp,
		long ndata, long marge, int *xx, int *yy, int nn,
		long *indpix, int npix, double *PNd){

	long ii, k, ll;

	double *Nd;
	fftw_complex *Ndf;
	fftw_plan fftplan;


	Ndf = new fftw_complex[ndata/2+1];
	Nd = new double[ndata];


	for (k=0;k<ndata/2+1;k++){
		Ndf[k][0] = fdata[k][0]*Nk[k];
		Ndf[k][1] = fdata[k][1]*Nk[k];
	}

	fftplan = fftw_plan_dft_c2r_1d(ndata, Ndf, Nd, FFTW_ESTIMATE);
	fftw_execute(fftplan);
	fftw_destroy_plan(fftplan);




	for (ii=0;ii<ndata;ii++){
		if ((ii < 0) || (ii >= ndata)){
			PNd[npix-2] += Nd[ii];
		} else {
			if (rejectsamp[ii] == 0){
				if (binsamp[ii] == 1){
					PNd[npix-2] += Nd[ii];
				} else {
					ll = indpix[yy[ii]*nn + xx[ii]];
					PNd[ll] += Nd[ii];
				}
			} else {
				PNd[npix-1] += Nd[ii];
			}
		}
	}

	delete [] Ndf;
	delete [] Nd;


}


 */





void compute_PtNP(double *Nk, unsigned char *rejectsamp, unsigned char *binsamp, long ndata,
		int *xx, int *yy, int nn, long *indpix,
		int npix, double f_lppix, double *PtNP){


	long jj, ll, ll2, indPtNP;
	int *pixpos;
	long *jj_sqr;


	//fft stuff
	fftw_complex  *Nk_;
	double *N_;
	fftw_plan fftplan;

	Nk_ = new fftw_complex[ndata/2+1];
	N_ = new double[ndata];
	pixpos = new int[ndata];
	jj_sqr = new long[npix];


	// N^-1
	for (long k=0;k<ndata/2+1;k++){
		Nk_[k][0] = 1.0/Nk[k]/(double)ndata/(double)ndata;
		Nk_[k][1] = 0.0;
	}
	fftplan = fftw_plan_dft_c2r_1d(ndata, Nk_, N_, FFTW_ESTIMATE);
	fftw_execute(fftplan);




	for (long ii=0;ii<ndata;ii++){
		if ((ii < 0) || (ii >= ndata)){
			pixpos[ii] = npix-2;
		} else {
			if (rejectsamp[ii] == 0){
				if (binsamp[ii] == 1){
					pixpos[ii] = npix-2;
				} else {
					pixpos[ii] = indpix[yy[ii]*nn + xx[ii]];
				}
			}
			else {
				pixpos[ii] = npix-1;
			}
		}
	}


	for (long ii=0;ii<npix;ii++)
		jj_sqr[ii] = ii*(ii+1)/2;



	for (long ii=0;ii<ndata;ii++){
		ll = pixpos[ii];
		ll2 = ll*(ll+1)/2;
		if (ii){
			for (long kk=MAX(ii-(ndata)/MAX(2,int(f_lppix+0.5)),0);kk<ii;kk++){
				//for (kk=0;kk<=ii;kk++){
				jj = pixpos[kk];
				if (ll < jj){
					indPtNP = jj*(jj+1)/2 + ll;
					//indPtNP = jj_sqr[jj] + ll;
				}else{
					indPtNP = ll2 + jj;
				}
				PtNP[indPtNP] += N_[ii-kk];
				//if (ii == kk) PtNP[indPtNP] -= N_[ii-kk]/2.0;
			}
		}
		PtNP[ll2+ll] += N_[0]/2.0;
		if ((ii % 20000) == 0)
			printf("%lf \n",pow((double)ii/double(ndata),2));
	}


	delete[] N_;
	delete[] Nk_;
	delete[] pixpos;
	delete[] jj_sqr;

	//clean up
	fftw_destroy_plan(fftplan);



}





/*

void compute_PtNP_frac(double *Nk, unsigned char *rejectsamp, unsigned char *binsamp, long ndata,
		int *xx, int *yy, int nn, long *indpix,
		int npix, double f_lppix, double *PtNP, int nfrac, int ifrac){


	long ii, k, jj, kk, ll, ii2, indPtNP;
	long *pixpos;
	long *jj_sqr;
	long ndataf;
	long count;
	long *pixtosamp_select;


	long indmin = ifrac*npix/nfrac;
	long indmax = (ifrac+1)*npix/nfrac;


	//fft stuff
	fftw_complex  *Nk_;
	double *N_;
	fftw_plan fftplan;

	Nk_ = new fftw_complex[ndata/2+1];
	N_ = new double[ndata];
	pixpos = new long[ndata];
	jj_sqr = new long[npix];
	pixtosamp_select = new long[ndata];


	// N^-1
	for (k=0;k<ndata/2+1;k++){
		Nk_[k][0] = 1.0/Nk[k]/(double)ndata/(double)ndata;
		Nk_[k][1] = 0.0;
	}
	fftplan = fftw_plan_dft_c2r_1d(ndata, Nk_, N_, FFTW_ESTIMATE);
	fftw_execute(fftplan);






	for (ii=0;ii<ndata;ii++){
		if ((ii < 0) || (ii >= ndata)){
			pixpos[ii] = npix-2;
		} else {
			if (rejectsamp[ii] == 0){
				if (binsamp[ii] == 1){
					pixpos[ii] = npix-2;
				} else {
					pixpos[ii] = indpix[yy[ii]*nn + xx[ii]];
				}
			}
			else {
				pixpos[ii] = npix-1;
			}
		}
	}


	for (ii=0;ii<npix;ii++)
		jj_sqr[ii] = ii*(ii+1)/2;





	count=0;
	for (ii=0;ii<ndata;ii++){
		if ((pixpos[ii] >= indmin) && (pixpos[ii] < indmax)){
			pixtosamp_select[count] = ii;
			count++;
		}
	}



	ndataf = (ndata)/MAX(2,int(f_lppix+0.5));

	for (ii=0;ii<count;ii++){
		ii2 = pixtosamp_select[ii];
		ll = pixpos[ii2]-indmin;
		if (ll > indmax-indmin)
			printf("ALERT ll = %ld\n",ll);
		for (kk=MAX(ii2-ndataf,0);kk<MIN(ii2+ndataf,ndata);kk++){ //MIN(ii2+(ndata)/MAX(2,int(f_lppix+0.5)),ndata);kk++){
			jj = pixpos[kk];
			indPtNP = ll*npix + jj;
			PtNP[indPtNP] += N_[abs(ii2-kk)];
		}

		if ((int((double)ii/(double)count*(double)ndata) % 20000) == 0)
			printf("%lf \n",(double)ii/double(count));
	}



	delete[] N_;
	delete[] Nk_;
	delete[] pixpos;
	delete[] jj_sqr;
	delete[] pixtosamp_select;


	//clean up
	fftw_destroy_plan(fftplan);



}



 */




void compute_diagPtNP(double *Nk, long *samptopix, long ndata,
		int nn, long *indpix,
		int npix, double f_lppix, double *dPtNP){


	long kk2, ii2, ndataf;
	long *pixpos;
	long count, count_;
	long *pixtosamp;


	//fft stuff
	fftw_complex  *Nk_;
	double *N_;
	fftw_plan fftplan;

	Nk_ = new fftw_complex[ndata/2+1];
	N_ = new double[ndata];
	pixpos = new long[ndata];
	//pixtosamp = new long[ndata];


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



	data_compare = new long[ndata];
	pixtosamp = new long[ndata];



	for (long ii=0;ii<ndata;ii++)
		pixtosamp[ii] = ii;

	for (long ii=0;ii<ndata;ii++)
		data_compare[ii] = pixpos[ii];



	qsort(pixtosamp,ndata,sizeof(long),compare_global_array_long);
	qsort(data_compare,ndata,sizeof(long),compare_long);



	ndataf = (ndata)/MAX(2,int(f_lppix+0.5));

	count = 0;

	for (long ipix=data_compare[0];ipix<npix;ipix++){

		count_ = count;

		while((count < ndata) && (data_compare[count] == ipix))
			count++;

		if (count-count_ > 0){
			for (long ii=count_;ii<count;ii++){
				ii2 = pixtosamp[ii];
				if ((ipix == npix-1) || (ipix == npix-2)){ //This is just to avoid spending to much time computing this pixel which could contain a lot of data
					dPtNP[ipix] += N_[0];
					//printf("TEST");
				} else {
					for (long kk=count_;kk<count;kk++){
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
	delete[] pixpos;
	delete[] pixtosamp;
	delete[] data_compare;


	//clean up
	fftw_destroy_plan(fftplan);


}






void compute_diagPtNPCorr(double *Nk, long *samptopix, long ndata,
		int nn, long *indpix,
		int npix, double f_lppix, double *dPtNP){


	long kk2, ii2, ndataf;
	long *pixpos;
	long count, count_;
	long *pixtosamp;


	//fft stuff
	fftw_complex  *Nk_;
	double *N_;
	fftw_plan fftplan;

	Nk_ = new fftw_complex[ndata/2+1];
	N_ = new double[ndata];
	pixpos = new long[ndata];



	// N^-1
	for (long k=0;k<ndata/2+1;k++){
		Nk_[k][0] = abs(Nk[k]);
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




	data_compare = new long[ndata];
	pixtosamp = new long[ndata];



	for (long ii=0;ii<ndata;ii++)
		pixtosamp[ii] = ii;

	for (long ii=0;ii<ndata;ii++)
		data_compare[ii] = pixpos[ii];



	qsort(pixtosamp,ndata,sizeof(long),compare_global_array_long);
	qsort(data_compare,ndata,sizeof(long),compare_long);



	ndataf = (ndata)/MAX(2,int(f_lppix+0.5));

	count = 0;

	for (long ipix=data_compare[0];ipix<npix;ipix++){

		count_ = count;

		while((count < ndata) && (data_compare[count] == ipix))
			count++;

		if (count-count_ > 0){
			for (long ii=count_;ii<count;ii++){
				ii2 = pixtosamp[ii];
				if ((ipix == npix-1) || (ipix == npix-2)){ //This is just to avoid spending to much time computing this pixel
					dPtNP[ipix] += N_[0];
					//printf("TEST");
				} else {
					for (long kk=count_;kk<count;kk++){
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
	delete[] pixpos;
	delete[] pixtosamp;
	delete[] data_compare;


	//clean up
	fftw_destroy_plan(fftplan);


}

/*


void compute_diagPtNPCorr_msk(double *Nk, unsigned char *mask, long iframe,
		unsigned char *rejectsamp, unsigned char *binsamp,
		long ndata, long marge, int *xx, int *yy, int nn,
		long *indpix, int npix, double f_lppix, double *dPtNP){


	long ii, k, kk, kk2, ipix, ii2, ndataf;
	long *pixpos;
	long count, count_;
	long *pixtosamp;


	//fft stuff
	fftw_complex  *Nk_;
	double *N_;
	fftw_plan fftplan;

	Nk_ = new fftw_complex[ndata/2+1];
	N_ = new double[ndata];
	pixpos = new long[ndata];
	//pixtosamp = new long[ndata];


	// N^-1
	for (k=0;k<ndata/2+1;k++){
		Nk_[k][0] = abs(Nk[k]);
		Nk_[k][1] = 0.0;
	}
	fftplan = fftw_plan_dft_c2r_1d(ndata, Nk_, N_, FFTW_ESTIMATE);
	fftw_execute(fftplan);






	for (ii=-marge;ii<ndata-marge;ii++){
		if ((ii < 0) || (ii >= ndata-2*marge)){
			pixpos[ii+marge] = npix-2;
		} else {
			if (rejectsamp[ii] == 0){
				if (binsamp[ii] == 1){
					pixpos[ii+marge] = npix-2;
				} else {
					if (mask[yy[ii]*nn + xx[ii]] == 1){
						pixpos[ii+marge] = indpix[yy[ii]*nn + xx[ii]];
					} else {
						pixpos[ii+marge] = indpix[(iframe + 1) * nn*nn + (yy[ii]*nn + xx[ii])];
						//printf("%d\n",indpix[(iframe + 1) * nn*nn + (yy[ii]*nn + xx[ii])]);
					}
				}
			}
			else {
				pixpos[ii+marge] = npix-1;
			}
		}
	}


	data_compare = new long[ndata];
	pixtosamp = new long[ndata];



	for (ii=0;ii<ndata;ii++)
		pixtosamp[ii] = ii;

	for (ii=0;ii<ndata;ii++)
		data_compare[ii] = pixpos[ii];



	qsort(pixtosamp,ndata,sizeof(long),compare_global_array_long);
	qsort(data_compare,ndata,sizeof(long),compare_long);



	ndataf = (ndata)/MAX(2,int(f_lppix+0.5));

	count = 0;


	//printf("NPIX = %d\n",npix);



	for (ipix=data_compare[0];ipix<npix;ipix++){

		count_ = count;

		while((count < ndata) && (data_compare[count] == ipix))
			count++;

		if (count-count_ > 0){
			for (ii=count_;ii<count;ii++){
				ii2 = pixtosamp[ii];
				if ((ipix == npix-2) || (ipix == npix-1)){ //This is just to avoid spending to much time computing this pixel
					dPtNP[ipix] += N_[0];
					//printf("TEST");
				} else {
					for (kk=count_;kk<count;kk++){
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
	delete[] pixpos;
	delete[] pixtosamp;
	delete[] data_compare;


	//clean up
	fftw_destroy_plan(fftplan);


}




 */



/*

void compute_diagPtNPCorr_new(double *Nk, unsigned char *rejectsamp,
		unsigned char *binsamp, long ndata,
		long marge, int *xx, int *yy, int nn, long *indpix,
		int npix, int npixmap, double f_lppix, double *dPtNP, long *countreject){


	long ii, k, kk, kk2, ipix, ii2, ndataf;
	long *pixpos;
	long count, count_;
	long *pixtosamp;


	//fft stuff
	fftw_complex  *Nk_;
	double *N_;
	fftw_plan fftplan;

	Nk_ = new fftw_complex[ndata/2+1];
	N_ = new double[ndata];
	pixpos = new long[ndata];
	//pixtosamp = new long[ndata];


	// N^-1
	for (k=0;k<ndata/2+1;k++){
		Nk_[k][0] = abs(Nk[k]);
		Nk_[k][1] = 0.0;
	}
	fftplan = fftw_plan_dft_c2r_1d(ndata, Nk_, N_, FFTW_ESTIMATE);
	fftw_execute(fftplan);






	for (ii=-marge;ii<ndata-marge;ii++){
		if ((ii < 0) || (ii >= ndata-2*marge)){
			pixpos[ii+marge] = npix-2;
		} else {
			if (rejectsamp[ii] == 0){
				if (binsamp[ii] == 1){
					pixpos[ii+marge] = npix-2;
				} else {
					pixpos[ii+marge] = indpix[yy[ii]*nn + xx[ii]];
				}
			}
			else {
				pixpos[ii+marge] = npixmap-1+ *countreject;
 *countreject = *countreject+1;
			}
		}
	}


	data_compare = new long[ndata];
	pixtosamp = new long[ndata];



	for (ii=0;ii<ndata;ii++)
		pixtosamp[ii] = ii;

	for (ii=0;ii<ndata;ii++)
		data_compare[ii] = pixpos[ii];



	qsort(pixtosamp,ndata,sizeof(long),compare_global_array_long);
	qsort(data_compare,ndata,sizeof(long),compare_long);



	ndataf = (ndata)/MAX(2,int(f_lppix+0.5));

	count = 0;

	for (ipix=data_compare[0];ipix<npix;ipix++){

		count_ = count;

		while((count < ndata) && (data_compare[count] == ipix))
			count++;

		if (count-count_ > 0){
			for (ii=count_;ii<count;ii++){
				ii2 = pixtosamp[ii];
				if ((ipix == npix-2) || (ipix == npix-1)){ //This is just to avoid spending to much time computing this pixel
					dPtNP[ipix] += N_[0];
					//printf("TEST");
				} else {
					for (kk=count_;kk<count;kk++){
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
	delete[] pixpos;
	delete[] pixtosamp;
	delete[] data_compare;


	//clean up
	fftw_destroy_plan(fftplan);


}






 */

/*

void compute_PtNP_corr(double *Nk, unsigned char *rejectsamp1, unsigned char *rejectsamp2,
		unsigned char *binsamp1, unsigned char *binsamp2,
		long ndata, long marge, int *xx1, int *yy1, int *xx2, int *yy2,
		int nn, long *indpix, int npix, double f_lppix, double *PtNP){





	long ii, k, jj, kk, ll, ll2, indPtNP;
	int *pixpos1, *pixpos2;

	//fft stuff
	fftw_complex  *Nk_;
	double *N_;
	fftw_plan fftplan;

	Nk_ = new fftw_complex[ndata/2+1];
	N_ = new double[ndata];
	pixpos1 = new int[ndata];
	pixpos2 = new int[ndata];



	// N^-1
	for (k=0;k<ndata/2+1;k++){
		Nk_[k][0] = Nk[k];
		Nk_[k][1] = 0.0;
		//printf("Nk[%d] = %lf\n",k,Nk[k]*(double)ndata*(double)ndata);
	}
	fftplan = fftw_plan_dft_c2r_1d(ndata, Nk_, N_, FFTW_ESTIMATE);
	fftw_execute(fftplan);




	for (ii=-marge;ii<ndata-marge;ii++){
		if ((ii < 0) || (ii >= ndata-2*marge)){
			pixpos1[ii+marge] = npix-2;
			pixpos2[ii+marge] = npix-2;
		} else {
			if (rejectsamp1[ii] == 0){
				if (binsamp1[ii] == 1){
					pixpos1[ii+marge] = npix-2;
				} else {
					pixpos1[ii+marge] = indpix[yy1[ii]*nn + xx1[ii]];
				}
			}
			else {
				pixpos1[ii+marge] = npix-1;
			}
			if (rejectsamp2[ii] == 0){
				if (binsamp2[ii] == 1){
					pixpos2[ii+marge] = npix-2;
				} else {
					pixpos2[ii+marge] = indpix[yy2[ii]*nn + xx2[ii]];
				}
			}
			else {
				pixpos2[ii+marge] = npix-1;
			}
		}
	}





	for (ii=0;ii<ndata;ii++){
		ll = pixpos1[ii];
		ll2 = ll*(ll+1)/2;
		for (kk=MAX(ii-(ndata)/MAX(2,int(f_lppix+0.5)),0);kk<=ii;kk++){
			//for (kk=0;kk<=ii;kk++){
			jj = pixpos2[kk];
			if (ll < jj){
				indPtNP = jj*(jj+1)/2 + ll;
			}else{
				indPtNP = ll2 + jj;
			}
			PtNP[indPtNP] += N_[ii-kk];
			if (ii == kk) PtNP[indPtNP] -= N_[ii-kk]/2.0;
		}
		if ((ii % 20000) == 0)
			printf("%lf \n",pow((double)ii/double(ndata),2));
	}


	delete[] N_;
	delete[] Nk_;
	delete[] pixpos1;
	delete[] pixpos2;

	//clean up
	fftw_destroy_plan(fftplan);





}






 */




void MapMakPreProcessData(double *data, unsigned char *flag, /*double *calp, */long ns, int napod,
		int orderpoly, double f_lppix, double *data_lp, double *bfilter, bool NORMLIN, bool NOFILLGAP,bool remove_polynomia, double *Ps){


	//long ii;
	double aa, bb;

	double *data_out, *data_out_lp;

	data_out = new double[ns];
	data_out_lp = new double[ns];


	//*********************************************************************

	if (NOFILLGAP == 0){
		//fill gaps with straight line
		fillgaps(data,ns,data_out,flag,0);
		for (long ii=0;ii<ns;ii++)
			data[ii] = data_out[ii];
	}

	if(remove_polynomia){
		//remove polynomia
		remove_poly(data,ns,orderpoly,data_out,0);



		//Used to correct from time varying calibration
		// now this is done to copy data_out in data
		for (long ii=0;ii<ns;ii++)
			data[ii] = data_out[ii]/**calp[ii/20]*/;
	}


	//linear prediction
	for (long ii=0;ii<ns;ii++)
		data_lp[ii] = data[ii];




	if (NORMLIN == 0){
		/// remove a baseline
		aa = (data_lp[ns-1]-data[0])/double(ns);
		bb = data_lp[0];
		for (long ii=0;ii<ns;ii++)
			data_lp[ii] -= aa*(double)ii+bb;
	}


	//Butterworth filter (if necessary)
	if (f_lppix > 0.0){
		butterworth(data_lp,ns,f_lppix,8,data_out_lp,bfilter,1,napod,0);
		for (long ii=0;ii<(ns);ii++)
			data_lp[ii] = data_out_lp[ii];
	} else{
		for (long ii=0;ii<(ns)/2+1;ii++)
			bfilter[ii] = 1.0;
	}



	if (Ps != NULL)
		for (long ii=0;ii<ns;ii++)
			data_lp[ii] = data_lp[ii] - Ps[ii];

	//******************* process gaps
	if (NOFILLGAP == 0){
		for (long ii=0;ii<ns;ii++)
			data_out[ii] = data_lp[ii];
		fillgaps(data_out,ns,data,flag,0);
		for (long ii=0;ii<ns;ii++)
			data_lp[ii] = data[ii];
	}


	if (Ps != NULL){
		//linear prediction

		for (long ii=0;ii<ns;ii++)
			//if (flag[ii] == 0)
			data_lp[ii] = data_lp[ii] + Ps[ii];
	}



	delete [] data_out;
	delete [] data_out_lp;



}







void noisepectrum_estim(double *data, int ns, double *ell, int nbins, double fsamp, double *bfilter, double *Nell, double *Nk){

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


	remove_poly(data,ns,4,datatemp2,0);
	apodwind = apodwindow(ns,ns/10);
	for (long ii=0;ii<ns;ii++)
		datatemp[ii] = datatemp2[ii]*apodwind[ii];


	//Fourier transform the data
	fftplan = fftw_plan_dft_r2c_1d(ns, datatemp, fdata, FFTW_ESTIMATE);
	fftw_execute(fftplan);


	totapod = 0.0;
	for (long ii=0;ii<ns;ii++)
		totapod += apodwind[ii]*apodwind[ii];


	//power spectrum
	for (long k=0;k<ns/2+1;k++){
		Nk[k] = pow(fdata[k][0],2) + pow(fdata[k][1],2);
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


	// interpol logarithmically the spectrum and filter
	binnedSpectrum2log_interpol(ell,Nell,bfiltertemp,nbins,ns,fsamp,Nk,NULL);



	//clean up
	delete [] count;
	delete [] datatemp;
	delete [] datatemp2;
	delete [] fdata;
	delete [] bfiltertemp;
	delete []  apodwind;

	fftw_destroy_plan(fftplan);

}







void noisecrosspectrum_estim(fftw_complex *fdata1, fftw_complex *fdata2, int ns, double *ell, int nbins, double fsamp, double *bfilter, double *Nell, double *Nk){

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








void readNSpectrum(char *nameSpfile, double *bfilter, long ns, double fsamp, double *Nk){

	FILE *fp;

	int nbins;
	//long ii;
	double dummy1, dummy2;

	double *SpN;
	double *ell;



	if ((fp = fopen(nameSpfile,"r")) == NULL){
		cerr << "ERROR: Can't find noise power spectra file, check -k or -K in command line. Exiting. \n";
		exit(1);
	}
	fscanf(fp,"%d",&nbins);


	SpN = new double[nbins];
	ell = new double[nbins+1];

	for (int ii=0;ii<nbins;ii++){
		fscanf(fp,"%lf %lf",&dummy1,&dummy2);
		ell[ii] = dummy1;
		SpN[ii] = dummy2;
	}
	fscanf(fp,"%lf",&dummy1);
	ell[nbins] = dummy1;
	fclose(fp);


	// interpolate logarithmically the noise power spectrum
	binnedSpectrum2log_interpol(ell,SpN,bfilter,nbins,ns,fsamp,Nk);



	delete[] SpN;
	delete[] ell;



}





/*
void readalldata(long ff, long ns, string field, string ra_field, string dec_field, string phi_field,
		string scerr_field, string flpoint_field, string dirfile,
		string bextension, string fextension, string cextension, double *data,
		double *calp, double *ra, double *dec, double *phi,
		double *scerr, unsigned char *flpoint, unsigned char *flag, int shift_data_to_point)
{


	long ii;

	char type = 'd';

	string bolofield;
	string flagfield;
	string calfield;


	bolofield = field+bextension;



	if (cextension != "NOCALP")
		calfield  = field+cextension;
	if (fextension != "NOFLAG")
		flagfield = field+fextension;



	read_data_std(dirfile, ff, shift_data_to_point, ns, data, bolofield,     type);
	if (cextension != "NOCALP"){
		read_data_std(dirfile, ff, 0, ns/20, calp, calfield, type);
	} else {
		//printf("NOCALP\n");
		for (ii=0;ii<ns/20;ii++)
			calp[ii] = 1.0;
	}
	read_data_std(dirfile, ff, 0, ns, ra,   ra_field,  type);
	read_data_std(dirfile, ff, 0, ns, dec,  dec_field, type);
	read_data_std(dirfile, ff, 0, ns, phi,  phi_field, type);
	read_data_std(dirfile, ff, 0, ns, scerr, scerr_field, type);
	read_data_std(dirfile, ff, 0, ns, flpoint, flpoint_field, 'c');
	for (ii=0;ii<ns;ii++)
		if (isnan(ra[ii]) || isnan(dec[ii]))
			flpoint[ii] = 1;
	if (fextension != "NOFLAG"){
		read_data_std(dirfile, ff, shift_data_to_point, ns, flag, flagfield,  'c');
	} else {
		//printf("NOFLAG\n");
		for (ii=0;ii<ns;ii++)
			flag[ii] = 0;
	}




}
 */

void deproject(double *S, long *indpix, long *samptopix, long ndata, long nn, long npix, double *Ps, int flgdupl, int factdupl, long ntotscan, long *indpsrc, long npixsrc){

	long ii,iframe;
	double a, b;
	//long ll;




	for (ii=0;ii<ndata;ii++){
		if ((ii < 0) || (ii >= ndata)){
			Ps[ii] = S[npix-2];
		} else {
			Ps[ii] = S[indpix[samptopix[ii]]];

			// in case we replaced flagged data
			if ((flgdupl == 2) && (samptopix[ii] >= nn*nn) && (samptopix[ii] < 2*nn*nn) && (indpix[samptopix[ii] - nn*nn] >= 0)){
				if (indpix[samptopix[ii] - nn*nn] >= 0){
					Ps[ii] = S[indpix[samptopix[ii]-nn*nn]];
				} else {
					a = 0.0;
					b = 0.0;
					if (ntotscan){
						for (iframe=0;iframe<ntotscan;iframe++){
							if (indpix[factdupl*nn*nn + indpsrc[samptopix[ii] - nn*nn] + iframe*npixsrc] >= 0){
								a += S[indpix[factdupl*nn*nn + indpsrc[samptopix[ii] - nn*nn] + iframe*npixsrc]];
								b++;
							}
						}
					}
					if (b > 0.5)
						Ps[ii] = a/b;
				}
			}
		}
	}

}




/*
void deproject_msk(double *S, unsigned char *mask, long *indpix, int *xx, int *yy, unsigned char *rejectsamp, unsigned char *binsamp, long ndata, long marge, long nn, long npix, long iframe, double *Ps){

  long ii, ll;

  for (ii=-marge;ii<ndata-marge;ii++){
    if ((ii < 0) || (ii >= ndata-2*marge)){
      Ps[ii+marge] = S[npix-2];
    } else {
      if (rejectsamp[ii] == 0){
	if (binsamp[ii] == 1){
	  Ps[ii+marge] = S[npix-2];
	} else {
	  if (mask[yy[ii]*nn + xx[ii]] == 1){
	    ll = indpix[yy[ii]*nn + xx[ii]];
	    Ps[ii+marge] = S[ll];
	  } else {
	    ll = indpix[(iframe + 1) * nn*nn + (yy[ii]*nn + xx[ii])];
	    Ps[ii+marge] = S[ll];
	  }
	}
      } else {
	Ps[ii+marge] = 0.0;
      }
    }
  }

}


 */


/*

void deproject_new(double *S, long *indpix, int *xx, int *yy, unsigned char *rejectsamp, unsigned char *binsamp, long ndata, long  , long nn, long npix, long npixmap, double *Ps, long *countreject){

  long ii, ll;

  for (ii=-marge;ii<ndata-marge;ii++){
    if ((ii < 0) || (ii >= ndata-2*marge)){
      Ps[ii+marge] = S[npix-1];
    } else {
      if (rejectsamp[ii] == 0){
	if (binsamp[ii] == 1){
	  Ps[ii+marge] = S[npix-1];
	} else {
	  ll = indpix[yy[ii]*nn + xx[ii]];
	  Ps[ii+marge] = S[ll];
	}
      } else {
	Ps[ii+marge] = S[npixmap-1+*countreject];
 *countreject = *countreject + 1;
      }
    }
  }

}



 */



int compare_global_array_long (const void *a, const void *b)
{

	const long *da = (const long *) a;
	const long *db = (const long *) b;

	return (data_compare[*da] > data_compare[*db]) - (data_compare[*da] < data_compare[*db]);
}

