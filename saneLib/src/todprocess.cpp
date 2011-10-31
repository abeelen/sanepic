


#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <fstream>
#include <cmath>
#include <fftw3.h>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>

#include "todprocess.h"
#include "cholesky.h"
#include "fitpoly.h"

using namespace std;





void init2D_double(double **A, long im, long jm, long nx, long ny, double val){


	for (long ii=im;ii<nx+im;ii++)
		for (long jj=jm;jj<ny+jm;jj++)
			A[ii][jj] = val;

}

void remove_poly(double y[], long ndata, int norder, double* yout, int* flag)
{
	long j;
	long ndint;
	double *sx, *sy;
	double* a;

	sx = new double[ndata];
	sy = new double[ndata];
	a = new double[norder+1];

	j=0;
	for (long i = 0; i < ndata; i++) {
		if(flag != NULL && (flag[i]!=0)) continue;
		sx[j] = (double)i;
		sy[j] = y[i];
		j++;
	}
	ndint = j;

	fitpoly(norder, ndint, sx, sy, a);


	double value =0.0;
	//remove best fit poly
	//	for (long i=0;i<ndata;i++) yout[i] = y[i];
	for (long i=0;i<ndata;i++){
		value = gsl_poly_eval (a, norder+1, (double)i);
		yout[i] = y[i] - value ;
		//	for (int pp=0;pp<=norder;pp++)
		//		yout[i] -= a[pp]*gsl_pow_int((double)i,pp);
	}



	delete [] sx;
	delete [] sy;
	delete [] a;


}



void butterworth(double y[], int ndata, double f_hp, int orderB, double *yout,
		double *bfilter, bool apodize, int napod, bool overwrite)
{

	double *apodwind;
	double powered;

	fftw_complex *fdata;
	fftw_plan fftplan;

	fdata = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*(ndata/2+1));

	//apodize if asked, and define plan for fft
	if (apodize){
		apodwind = apodwindow(ndata,napod);
		if (overwrite){
			for (int ii=0;ii<ndata;ii++) y[ii] = apodwind[ii] * y[ii];
			fftplan = fftw_plan_dft_r2c_1d(ndata, y, fdata, FFTW_ESTIMATE);
		} else{
			for (int ii=0;ii<ndata;ii++) yout[ii] = apodwind[ii] * y[ii];
			fftplan = fftw_plan_dft_r2c_1d(ndata, yout, fdata, FFTW_ESTIMATE);
		}
		delete[] apodwind;
	} else{
		fftplan = fftw_plan_dft_r2c_1d(ndata, y, fdata, FFTW_ESTIMATE);
	}


	fftw_execute(fftplan);
	fftw_destroy_plan(fftplan);

	//filter
	for (int ii=0;ii<ndata/2+1;ii++){
		powered = gsl_pow_int(double(ii)/f_hp, 2*orderB) ;
		bfilter[ii] = powered/(1.0+powered);//pow(double(ii)/f_hp, 2*orderB) /(1.0+pow(double(ii)/f_hp, 2*orderB));
		fdata[ii][0] = fdata[ii][0]*bfilter[ii]/ndata;
		fdata[ii][1] = fdata[ii][1]*bfilter[ii]/ndata;
	}

	if (overwrite){
		fftplan = fftw_plan_dft_c2r_1d(ndata, fdata, y, FFTW_ESTIMATE);
	}else{
		fftplan = fftw_plan_dft_c2r_1d(ndata, fdata, yout, FFTW_ESTIMATE);
	}
	fftw_execute(fftplan);
	fftw_destroy_plan(fftplan);

	delete [] fdata;


}




double* apodwindow(int ns, int nn)
{

	double *apodis;

	apodis = new double[ns];

	for (int ii=0;ii<ns;ii++){
		apodis[ii] = 1.0;
	}

	if (nn!=1){
		for (int ii=0;ii<nn;ii++){
			apodis[ii] = (sin(double(ii)/(nn-1.0)*M_PI - M_PI/2.0) + 1.0)/2.0;
		}
		for (int ii=ns-nn;ii<ns;ii++){
			apodis[ii] = (sin(double(ns-ii-1)/(nn-1.0)*M_PI - M_PI/2.0) + 1.0)/2.0;
		}
	}

	return apodis;

}



void binnedSpectrum2log_interpol(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode)
{

	// ell is an array of double, units are Hz

	int counttemp, f_hp;
	double ellmin, ellmax, kmin, kmax, a, b;
	double *ellm;
	double N_flp;

	// interpolate logarithmically the noise power spectrum

	ellm = new double[nbins];
	for (int ii=0;ii<nbins;ii++)
		ellm[ii] = exp((log(ell[ii+1])+log(ell[ii]))/2.0);

	counttemp = 0;
	ellmin = ellm[0];
	ellmax = ellm[1];
	kmin = ellmin*ns/fsamp;
	kmax = ellmax*ns/fsamp;

	a = (log(SpN[1]) - log(SpN[0]))/(log(kmax)-log(kmin));
	b = log(SpN[0]);


	if (mode == NULL){
		for (int k=0;k<ns/2+1;k++){
			while (double(k) > ellm[counttemp]*ns/fsamp && counttemp < nbins-1){
				counttemp++;
			}
			if (counttemp > 0){
				ellmin = ellm[counttemp-1];
				ellmax = ellm[counttemp];
				kmin = ellmin*ns/fsamp;
				kmax = ellmax*ns/fsamp;
				if ((abs(SpN[counttemp]) > 0) || (SpN[counttemp] > 0)){
					a = (log(SpN[counttemp]) -
							log(SpN[counttemp-1]))/(log(kmax)-log(kmin));
					b = log(SpN[counttemp-1]);
				} else {
					a = 0;
					b = 0;
				}
			}
			if (SpN[counttemp] > 0)
				Nk[k] = exp(a*(log((double)k)-log(kmin))+b)/double(ns);
			else {
				Nk[k] = 0.0;
			}
			//apply filter
			Nk[k] *= gsl_pow_2(bfilter[k]);

		}
		Nk[0] = Nk[1];



		f_hp = 0;
		while (bfilter[f_hp] <= 0.5) f_hp++;
		f_hp++;


		//give a lower limit to the spectrum
		for (int k=0;k<f_hp;k++) if (Nk[k] < Nk[f_hp]) Nk[k] = Nk[f_hp];
		//for (k=0;k<f_hp;k++) Nk[k] = Nk[f_hp];

		// suppress effect of aafilter on the Noise Sp
		for (int k=0;k<ns/2+1;k++) if (Nk[k] < Nk[ns/5]) Nk[k] = Nk[ns/5];




	} else {////// compute noise power spectrum for a given mode
		while (*mode > ellm[counttemp]*ns/fsamp && counttemp < nbins-1){
			counttemp++;
		}
		if (counttemp > 0){
			ellmin = ellm[counttemp-1];
			ellmax = ellm[counttemp];
			kmin = ellmin*ns/fsamp;
			kmax = ellmax*ns/fsamp;
			a = (log(SpN[counttemp]) -
					log(SpN[counttemp-1]))/(log(kmax)-log(kmin));
			b = log(SpN[counttemp-1]);
		}
		*Nk = exp(a*(log(*mode)-log(kmin))+b)/double(ns);
		*Nk *= gsl_pow_2(bfilter[int(*mode)]);


		f_hp = 0;
		while (bfilter[f_hp] <= 0.5) f_hp++;
		f_hp++;


		if (*mode < f_hp){

			counttemp = 0;
			while (f_hp > ellm[counttemp]*ns/fsamp && counttemp < nbins-1){
				counttemp++;
			}
			if (counttemp > 0){
				ellmin = ellm[counttemp-1];
				ellmax = ellm[counttemp];
				kmin = ellmin*ns/fsamp;
				kmax = ellmax*ns/fsamp;
				a = (log(SpN[counttemp]) -
						log(SpN[counttemp-1]))/(log(kmax)-log(kmin));
				b = log(SpN[counttemp-1]);
			}
			N_flp = exp(a*(log((double)f_hp)-log(kmin))+b)/double(ns);
			N_flp *= gsl_pow_2(bfilter[f_hp]);



			//give a lower limit to the spectrum
			if (*Nk < N_flp) *Nk = N_flp;

			// suppress effect of aafilter on the Noise Sp
			for (int k=0;k<ns/2+1;k++) if (Nk[k] < Nk[ns/5]) Nk[k] = Nk[ns/5];
		}
	}

	delete [] ellm;

}






void InvbinnedSpectrum2log_interpol(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode)
{
	// ell is an array of double, units are Hz

	int counttemp, f_hp;
	double ellmin, ellmax, kmin, kmax, a, b;
	//	double lkmin, lkmax;
	double *ellm;
	//	double *logSpN;
	//	double N_flp;

	a = 0.0;
	b = 0.0;


	ellm = new double[nbins];
	//	logSpN = new double[nbins];

	fill(ellm,ellm+nbins,0.0);
	//	fill(logSpN,logSpN+nbins,0.0);


	for (int ii=0;ii<nbins;ii++)
		ellm[ii] = exp((log(ell[ii+1])+log(ell[ii]))/2.0);


	counttemp = 0;
	ellmin = ellm[0];
	ellmax = ellm[1];
	kmin = ellmin*ns/fsamp;
	kmax = ellmax*ns/fsamp;
	//	lkmin = log(kmin);
	//	lkmax = log(kmax);


	//	if (mode == NULL){



	for (long k=1;k<=(long)kmin;k++){

		Nk[k] = SpN[0]/double(ns);

	}



	/////////////  linear interpolation
	// because, by convention, some spectrum can be negatives !!!
	for (int ibin=0;ibin<nbins-1;ibin++){
		ellmin = ellm[ibin];
		ellmax = ellm[ibin+1];
		kmin = ellmin*ns/fsamp;
		kmax = ellmax*ns/fsamp;
		if (abs(SpN[ibin]) > 0){
			a = (SpN[ibin+1] - SpN[ibin])/(kmax-kmin)/double(ns);
			b = SpN[ibin]/double(ns);
		} else {
			a = 0.0;
			b = 0.0;
		}
		for (long k=long(kmin+1);k<=long(kmax);k++){
			Nk[k] = (a*((double)k-kmin)+b);

			//apply filter
			Nk[k] *= gsl_pow_2(bfilter[k]);

		}
	}



	for (long k=long(kmax);k<ns/2+1;k++){

		Nk[k] = SpN[nbins-1]*gsl_pow_2(bfilter[k])/double(ns);

	}


	Nk[0] = Nk[1];



	f_hp = 0;
	while (bfilter[f_hp] > 2.0) f_hp++;
	f_hp++;



	//give a lower limit to the spectrum
	for (int k=0;k<f_hp;k++) Nk[k] = Nk[f_hp];

	// suppress effect of aafilter on the Noise Sp
	for (long k=ns/20;k<ns/2+1;k++) Nk[k] = Nk[ns/20];




	//	} else {////// compute noise power spectrum for a given mode
	//		while (*mode > ellm[counttemp]*ns/fsamp && counttemp < nbins-1){
	//			counttemp++;
	//		}
	//		if (counttemp > 0){
	//			ellmin = ellm[counttemp-1];
	//			ellmax = ellm[counttemp];
	//			kmin = ellmin*ns/fsamp;
	//			kmax = ellmax*ns/fsamp;
	//			a = (log(SpN[counttemp]) -
	//					log(SpN[counttemp-1]))/(log(kmax)-log(kmin));
	//			b = log(SpN[counttemp-1]);
	//		}
	//		*Nk = exp(a*(log(*mode)-log(kmin))+b)/double(ns);
	//		*Nk *= gsl_pow_2(bfilter[int(*mode)]);
	//
	//
	//		f_hp = 0;
	//		while (bfilter[f_hp] > 2.0) f_hp++;
	//		f_hp++;
	//
	//
	//		if (*mode < f_hp){
	//
	//			counttemp = 0;
	//			while (f_hp > ellm[counttemp]*ns/fsamp && counttemp < nbins-1){
	//				counttemp++;
	//			}
	//			if (counttemp > 0){
	//				ellmin = ellm[counttemp-1];
	//				ellmax = ellm[counttemp];
	//				kmin = ellmin*ns/fsamp;
	//				kmax = ellmax*ns/fsamp;
	//				a = (log(SpN[counttemp]) -
	//						log(SpN[counttemp-1]))/(log(kmax)-log(kmin));
	//				b = log(SpN[counttemp-1]);
	//			}
	//			N_flp = exp(a*(log((double)f_hp)-log(kmin))+b)/double(ns);
	//			N_flp *= gsl_pow_2(bfilter[f_hp]);
	//
	//
	//
	//			//give a lower limit to the spectrum
	//			if (*Nk < N_flp) *Nk = N_flp;
	//
	//		}
	//	}





	delete [] ellm;
	//	delete [] logSpN;

}

void fillgaps2(double data[], long ns, double* yout,  int* flag, int taille){

	int indic=0;
	std::vector<long> buf_0, buf_1;
	long p=0;
	double *sx, *sy, *a;
	long p_beg_init = 0, p_end_init = 0;

	a = new double[2];
	a[0]=0.0;
	a[1]=0.0;

	sx= new double[taille];
	sy= new double[taille];

	for(long ff=0;ff<ns;ff++)
		if(flag[ff]==0)
			yout[ff]=data[ff];

	buf_0.push_back(0);
	buf_1.push_back(0);

	if(flag[0]!=0){
		buf_0.push_back(0);
		indic=1;
		p=0;
	}

	for(long ii=1; ii<ns;ii++) {
		if( indic == 1 && flag[ii] == 0 ){
			buf_1.push_back(ii-p);
			indic=0;
			p=ii;
		}
		if ( indic == 0 && flag[ii] != 0 ) {
			buf_0.push_back(ii-p);
			indic=1;
			p=ii;
		}
	}


	if(indic==0)
		buf_0.push_back(ns-p);
	else
		buf_1.push_back(ns-p);

	buf_0.erase(buf_0.begin());
	buf_1.erase(buf_1.begin());


	long n=(long)buf_0.size();

	long flag_precedent=0;
	long init=0;

	if(buf_0[0]==0){
		p_beg_init=0;
		p_end_init=buf_1[0];
	}

	long i=0;
	long reste;
	long p_beg=0, p_end=0;
	long miss=0;

	while(i<n-1){

		fill(sx,sx+taille,0.0);
		fill(sy,sy+taille,0.0);

		miss=0;

		if(buf_0[i]>=taille/2){ // enough samples before gap
			p_beg=buf_0[i]-taille/2+flag_precedent; // get first indice in data array (we don't take flagged samples)
			for(int gg=0;gg<taille/2;gg++){ // fill x and y arrays for fitpoly with the good samples found
				sx[gg]=gg;
				sy[gg]=data[p_beg+gg];
			}
			reste=0; // we have taille/2 samples before gap so no need to get more than taille/2 after gap : reste=0

		}else{ // not enough samples before gap

			p_beg=flag_precedent+init; // get first indice in data array (we don't take flagged samples)
			reste=taille/2-buf_0[i]; // reste = number of samples to take after gap because there were not enough before
			for(int gg=0;gg<buf_0[i];gg++){ // fill x and y arrays for fitpoly with the good samples found
				sx[gg]=gg;
				sy[gg]=data[p_beg+gg];
			}
		}
		long old_i=i;
		long sum_0=0;
		if(buf_0[i+1]>=(taille/2+reste)){ // enough samples just after the gap
			p_end=taille+buf_1[i]+p_beg-init-1;
			for(int gg=taille/2-reste;gg<taille;gg++){
				sx[gg]=buf_1[i]+gg;
				sy[gg]=data[p_beg+buf_1[i]-init+gg];
			}
		}else{ // not enough samples just after the gap, consider a new gap AND the good samples after it to have enough good samples
			while(sum_0<taille/2+reste){ // check how many gap we need to "over-fill"
				i++;
				if(i>n-1){
					i--;
					break;
				}
				sum_0 += buf_0[i];
			}
			if(sum_0>taille/2+reste){ // we have found enough good samples
				p_end=taille+p_beg-1;
				for(long ff=old_i;ff<i;ff++)
					p_end+=buf_1[ff];

			}else{ // end of flag array and we have not enough (<taille) good samples for fitpoly

				p_end=taille/2-reste+sum_0+p_beg-1; // case last buf_0 < taille/2+reste
				for(long ff=old_i;ff<i;ff++)
					p_end+=buf_1[ff];
			}
			i--;

			long indice=0;
			for(long gg=p_beg;gg<=p_end;gg++){
				if(flag[gg]>0){					continue;
				}
				sx[indice]=gg-p_beg;
				sy[indice]=data[gg];
				indice++;
			}
			miss = taille-indice;
		}


		fitpoly(1, taille-miss, sx, sy, a);

		for (long j=p_beg;j<=p_end;j++){
			if(flag[j]!=0)
				yout[j] = a[0]+a[1]*(j-p_beg);
		}

		if(init!=0)
			for(long j=p_beg_init;j<p_end_init;j++)
				yout[j] = a[0]+a[1]*(j-p_beg_init);


		i++;
		if(i<n-1){
			flag_precedent=0;
			for(long gg=0;gg<i;gg++)
				flag_precedent+=buf_0[gg] + buf_1[gg];
		}

		init=0;

	}

	if(n==(long)buf_1.size()){
		if(buf_0[n]>=taille){
			p_beg=ns-buf_1[n]-taille;
			p_end=p_beg+taille-1;
			for(int gg=0;gg<taille;gg++){ // fill x and y arrays for fitpoly with the good samples found
				sx[gg]=gg;
				sy[gg]=data[p_beg+gg];
			}

			fitpoly(1, taille, sx, sy, a);
			for (long j=p_beg;j<=ns-1;j++){
				if(flag[j]>0)
					yout[j] = a[0]+a[1]*(j-p_beg);
			}

		}else{
			p_end=ns-1;
			for (int j=p_beg;j<=p_end;j++){
				if(flag[j]>0)
					yout[j] = a[0]+a[1]*(j-p_beg);
			}
			// use last poly and fit with it
		}
	}

	delete [] sx;
	delete [] sy;
	delete [] a;


}
