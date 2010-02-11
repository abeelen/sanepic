


#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <fstream>
#include <cmath>
#include <fftw3.h>
#include <vector>
#include <gsl/gsl_math.h>

#include "todprocess.h"

using namespace std;





void init2D_double(double **A, long im, long jm, long nx, long ny, double val){

	//long ii, jj;

	for (long ii=im;ii<nx+im;ii++)
		for (long jj=jm;jj<ny+jm;jj++)
			A[ii][jj] = val;

}



/*
void  minmax(double* data, int ndata, double *min, double *max, int *posmin, int *posmax, unsigned char *flag)
{
  int j;

 *min = data[0];
 *max = data[0];
 *posmin = 0;
 *posmax = 0;
  j=0;
  if (flag != NULL){
    while (flag[j] != 0){
 *min = data[j+1];
 *max = data[j+1];
 *posmin = j+1;
 *posmax = j+1;
      j++;
    }
  }

  for(int k = 1; k < ndata; k++) {
    if ((flag == NULL || flag[k] == 0) && (isnan(data[k]) == 0)){
      if(data[k] < *min){
 *min = data[k];
 *posmin = k;
      } else if(data[k] > *max){
 *max = data[k];
 *posmax = k;
      }
    }
  }
}
 */




void polynomia(double x, double y[], int dma)
{
	//int i;

	for(int i=1;i<=dma;i++){
		y[i] = pow(x,i-1);
	}

}



void dpolyfit(double x[], double y[], int ndata, int norder, double *a)
{
	//int i;
	int ma;
	double chisq;
	double *sig, *b;
	int *ia;
	double** covar;

	ma = norder+1;

	sig = new double[ndata];
	ia  = new int[ma];
	covar = dmatrix(1,ma,1,ma);
	fill(sig,sig+ndata,0.0);
	fill(ia,ia+ma,0);


	//initialize sig to 1
	for (int i=0;i<ndata;i++){
		sig[i] = 1.0;
	}
	//set to estimate all parameters
	for (int i=0;i<ma;i++){
		ia[i]  = 1;
	}

	b=a-1;
	dlfit(x-1,y-1,sig-1,ndata,b,ia-1,ma,covar,&chisq,polynomia);
	a = b+1;

	delete [] sig;
	delete [] ia;
	//delete(sig);
	//delete(ia);
	free_dmatrix(covar,1,ma,1,ma);

}

void remove_poly(double y[], int ndata, int norder, double* yout, short* flag)
{
	int j;
	int ndint;
	double *sx, *sy;
	double* a;

	sx = new double[ndata];
	sy = new double[ndata];
	a = new double[norder+1];

	j=0;
	for (int i = 0; i < ndata; i++) {
		//if(flag != NULL && flag[i]) continue;
		if(flag != NULL && (flag[i]==1)) continue;
		sx[j] = i;
		sy[j] = y[i];
		j++;
	}
	ndint = j;

	dpolyfit(sx,sy,ndint,norder,a);

	//remove best fit poly
	for (int i=0;i<ndata;i++) yout[i] = y[i];
	for (int i=0;i<ndata;i++)
		for (int pp=0;pp<=norder;pp++)
			yout[i] -= a[pp]*gsl_pow_int((double)i,pp);//pow((double)i,pp);


	delete [] sx;
	delete [] sy;
	delete [] a;
	/* delete(sx);
  delete(sy);
  delete(a);*/

}



void butterworth(double y[], int ndata, double f_lp, int orderB, double *yout,
		double *bfilter, bool apodize, int napod, bool overwrite)
{

	// int ii;
	double *apodwind;
	double powered;

	fftw_complex *fdata;
	fftw_plan fftplan;

	//fftw_complex fdata[ns/2+1];
	fdata = new fftw_complex[ndata/2+1];

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
		powered = gsl_pow_int(double(ii)/f_lp, 2*orderB) ;
		bfilter[ii] = powered/(1.0+powered);//pow(double(ii)/f_lp, 2*orderB) /(1.0+pow(double(ii)/f_lp, 2*orderB));
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

	//int ii;
	double *apodis;

	apodis = new double[ns];

	for (int ii=0;ii<ns;ii++){
		apodis[ii] = 1.0;
	}

	if (nn){
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

	/////////////////
	//
	// ell is an array of double, units are Hz
	//
	////////////////


	int counttemp, f_lp;
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
			if ((SpN[counttemp] > 0) || (SpN[counttemp] > 0))
				Nk[k] = exp(a*(log((double)k)-log(kmin))+b)/double(ns);
			else {
				Nk[k] = 0.0;
			}
			//apply filter
			Nk[k] *= gsl_pow_2(bfilter[k]);

		}
		Nk[0] = Nk[1];



		f_lp = 0;
		while (bfilter[f_lp] <= 0.5) f_lp++;
		f_lp++;


		//give a lower limit to the spectrum
		for (int k=0;k<f_lp;k++) if (Nk[k] < Nk[f_lp]) Nk[k] = Nk[f_lp];
		//for (k=0;k<f_lp;k++) Nk[k] = Nk[f_lp];

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


		f_lp = 0;
		while (bfilter[f_lp] <= 0.5) f_lp++;
		f_lp++;


		if (*mode < f_lp){

			counttemp = 0;
			while (f_lp > ellm[counttemp]*ns/fsamp && counttemp < nbins-1){
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
			N_flp = exp(a*(log((double)f_lp)-log(kmin))+b)/double(ns);
			N_flp *= gsl_pow_2(bfilter[f_lp]);



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

	//TODO : Comment all variable, understand what this is really doing, test & optimize
	/////////////////
	//
	// ell is an array of double, units are Hz
	//
	////////////////


	int counttemp, f_lp;
	double ellmin, ellmax, kmin, kmax, lkmin, lkmax, a, b;
	double *ellm, *logSpN;
	double N_flp;

	a = 0.0;
	b = 0.0;


	/*
  counttemp = 0;

  if (mode == NULL){
    for (k=0;k<ns/2+1;k++){
      while (double(k) > ell[counttemp]*ns/fsamp && counttemp < nbins){
	counttemp++;
      }
      if (counttemp > 0){
	Nk[k] = SpN[counttemp-1]/double(ns);
      } else {
	Nk[k] = SpN[0]/double(ns);
      }
      //apply filter
      Nk[k] *= pow(bfilter[k],2);

    }
    Nk[0] = Nk[1];



    f_lp = 0;
    while (bfilter[f_lp] <= 0.5) f_lp++;
    f_lp++;

    //give a lower limit to the spectrum
    //for (k=0;k<f_lp;k++) if (Nk[k] < Nk[f_lp]) Nk[k] = Nk[f_lp];
    for (k=0;k<f_lp;k++) Nk[k] = Nk[f_lp];


  } else {////// compute noise power spectrum for a given mode
    while (*mode > ell[counttemp]*ns/fsamp && counttemp < nbins){
      counttemp++;
    }
    if (counttemp > 0){
	 *Nk = SpN[counttemp-1]/double(ns);
    } else {
	 *Nk = SpN[0]/double(ns);
    }

    f_lp = 0;
    while (bfilter[f_lp] <= 0.5) f_lp++;
    f_lp++;

    if (*mode < f_lp){

      counttemp = 0;
      while (f_lp > ell[counttemp]*ns/fsamp && counttemp < nbins){
	counttemp++;
      }
      if (counttemp > 0){
	N_flp = SpN[counttemp-1]/double(ns);
      } else {
	N_flp = SpN[0]/double(ns);
      }
      N_flp *= pow(bfilter[f_lp],2);



      //give a lower limit to the spectrum
      if (*Nk < N_flp) *Nk = N_flp;

    }
  }

	 */





	ellm = new double[nbins];
	logSpN = new double[nbins];

	fill(ellm,ellm+nbins,0.0);
	fill(logSpN,logSpN+nbins,0.0);


	for (int ii=0;ii<nbins;ii++)
		ellm[ii] = exp((log(ell[ii+1])+log(ell[ii]))/2.0);


	counttemp = 0;
	ellmin = ellm[0];
	ellmax = ellm[1];
	kmin = ellmin*ns/fsamp;
	kmax = ellmax*ns/fsamp;
	lkmin = log(kmin);
	lkmax = log(kmax);


	if (mode == NULL){



		for (long k=1;k<=(long)kmin;k++){

			Nk[k] = SpN[0]/double(ns);

		}




		//		if (0){
		//
		//			for (int ii=0;ii<nbins;ii++)
		//				logSpN[ii] = log(abs(SpN[ii]));
		//			//////////////////////////////////////////  log-log interpolation
		//			for (int ibin=0;ibin<nbins-1;ibin++){
		//				ellmin = ellm[ibin];
		//				ellmax = ellm[ibin+1];
		//				kmin = ellmin*ns/fsamp;
		//				kmax = ellmax*ns/fsamp;
		//				lkmin = log(kmin);
		//				lkmax = log(kmax);
		//				if (abs(SpN[ibin]) > 0){
		//					a = (logSpN[ibin+1] - logSpN[ibin])/(lkmax-lkmin);
		//					b = logSpN[ibin];
		//				} else {
		//					a = 0;
		//					b = 0;
		//				}
		//				for (long k=(long)kmin;k<long(kmax);k++){
		//					if (SpN[ibin] > 0)
		//						Nk[k] = exp(a*(log((double)k)-lkmin)+b)/double(ns);
		//					else {
		//						if (SpN[ibin] < 0){
		//							Nk[k] = -exp(a*(log((double)k)-lkmin)+b)/double(ns);
		//						} else {
		//							Nk[k] = 0.0;
		//						}
		//					}
		//					//apply filter
		//					Nk[k] *= gsl_pow_2(bfilter[k]);
		//
		//				}
		//			}
		//		}else{





		/////////////  linear interpolation
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
		//		}



		for (long k=long(kmax);k<ns/2+1;k++){

			Nk[k] = SpN[nbins-1]*gsl_pow_2(bfilter[k])/double(ns);

		}


		Nk[0] = Nk[1];








		f_lp = 0;
		while (bfilter[f_lp] > 2.0) f_lp++;
		f_lp++;



		//give a lower limit to the spectrum
		//for (k=0;k<f_lp;k++) if (Nk[k] < Nk[f_lp]) Nk[k] = Nk[f_lp];
		for (int k=0;k<f_lp;k++) Nk[k] = Nk[f_lp];

		// suppress effect of aafilter on the Noise Sp
		for (long k=ns/20;k<ns/2+1;k++) Nk[k] = Nk[ns/20];




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


		f_lp = 0;
		while (bfilter[f_lp] > 2.0) f_lp++;
		f_lp++;


		if (*mode < f_lp){

			counttemp = 0;
			while (f_lp > ellm[counttemp]*ns/fsamp && counttemp < nbins-1){
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
			N_flp = exp(a*(log((double)f_lp)-log(kmin))+b)/double(ns);
			N_flp *= gsl_pow_2(bfilter[f_lp]);



			//give a lower limit to the spectrum
			if (*Nk < N_flp) *Nk = N_flp;

			// suppress effect of aafilter on the Noise Sp
			//for (k=0;k<ns/2+1;k++) if (Nk[k] < Nk[ns/100]) Nk[k] = Nk[ns/100];
		}
	}





	delete [] ellm;
	delete [] logSpN;

}



void InvbinnedSpectrum2bis(double* ell, double* SpN, double* bfilter, int nbins, int ns, double fsamp, double* Nk, double* mode)
{

	/////////////////
	//
	// ell is an array of double, units are Hz
	//
	////////////////


	int counttemp, f_lp;
	double ellmin, ellmax, kmin, kmax = 0, a, b;
	double *ellm, *logSpN;
	double SpN_max;

	a = 0.0;
	b = 0.0;


	ellm = new double[nbins];
	//ellm2= new double[nbins];
	logSpN = new double[nbins];



	for (int ii=0;ii<nbins;ii++){
		//ellm[ii] = exp((log(ell[ii+1])+log(ell[ii]))/2.0);
		ellm[ii] = sqrt(ell[ii+1])*sqrt(ell[ii]); // ajout mat 19/08 : faster ?
		//cout << "ellm : " << ellm[ii] << " " << "ellm2 : " << ellm2[ii] << endl;
	}


	counttemp = 0;
	ellmin = ellm[0];
	ellmax = ellm[1];
	kmin = ellmin*ns/fsamp;
	//kmax = ellmax*ns/fsamp;
	//lkmin = log(kmin);
	//lkmax = log(kmax);




	for (long k=1;k<=(long)kmin;k++){

		Nk[k] = SpN[0]/double(ns);

	}
	//cout << endl;

	f_lp = 0;
	while (bfilter[f_lp] > 2.0) f_lp++;
	f_lp++;

	/////////////  linear interpolation
	for (int ibin=0;ibin<nbins-1;ibin++){
		ellmin = ellm[ibin];
		ellmax = ellm[ibin+1];
		kmin = ellmin*ns/fsamp;
		kmax = ellmax*ns/fsamp;
		if((long)kmax<f_lp) // ajout mat 21/08
			continue;
		if((long)kmin>(ns/20)) // ajout mat 21/08
			break;

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
			//cout << bfilter[k] << " ";

		}
	}

	//cout << endl;
	//getchar();

	SpN_max=SpN[nbins-1];
	//cout << "bfilter avec kmax = " <<  long(kmax) << endl;
	//for (long k=long(kmax);k<ns/2+1;k++){
	for (long k=long(kmax);k<ns/20+1;k++){

		Nk[k] = SpN_max*gsl_pow_2(bfilter[k])/double(ns);
		//cout << bfilter[k] << " ";

	}

	//cout << endl;
	//getchar();

	Nk[0] = Nk[1];



	//give a lower limit to the spectrum
	for (int k=0;k<f_lp;k++) if (Nk[k] < Nk[f_lp]) Nk[k] = Nk[f_lp];
	//	for (int k=0;k<f_lp;k++) Nk[k] = Nk[f_lp];

	// suppress effect of aafilter on the Noise Sp
	for (long k=ns/20+1;k<ns/2+1;k++) Nk[k] = Nk[ns/20];






	delete [] ellm;
	delete [] logSpN;

}







/*
void deconv_antialias(double y[], int ndata, double f_lp, double* yout, bool apodize)
{

  int k, j;
    double omega, filter, ofilter;

  double omega_nyquist = 2.0*M_PI*0.5*HUNDRED_HZ_CLK;
  double omega_lp = 2.0*M_PI*f_lp;
  double omega_max = omega_lp + 0.2*omega_nyquist;

  const double stage[] = {208.0/ADC_CLK, 175.0/ADC_CLK,
			 147.0/ADC_CLK, 123.0/ADC_CLK};

  Vector v(0, ndata-1, 0.0);

  for(j = 0; j <= ndata-1; j++) {
    v[j] = y[j];

    // apodize with a Hanning window function
    if(apodize) v[j] *= 1.0 - cos(2.0*M_PI * (double)j/(double)ndata);
  }

  FFT(v);

  // Filter here
  for (ofilter = 1.0,  k = 0; k <= ndata/2 ; k++) {
    omega = 2.0*(double)k/double(ndata) * omega_nyquist;

    if (omega < omega_lp) {
      for(filter = 1.0, j = 0; j < 4; j++) {
	if (omega != 0) filter /= sin(0.5*omega*stage[j])/(0.5*omega*stage[j]);
      }
      ofilter = filter;


    } else if (omega < omega_max){
      // Low pass filter
      filter = ofilter *
	(0.5 + 0.5*cos(M_PI*(omega - omega_lp)/(omega_max - omega_lp)));
    } else {
      filter = 0.0;
    }


    // negative frequency treatment
    if(k == 0 || k == ndata/2) {
      //cout << omega/(2*M_PI) << " " << fabs(v[k]) << " ";
      v[k] *= filter;
      //cout <<  fabs(v[k]);
    } else {
      //cout << omega/(2*M_PI) << " " << sqrt(SQR(v[k]) + SQR(v[out_size - k])) << " ";
      v[k] *= filter;
      v[ndata - k] *= filter;
      //cout << sqrt(SQR(v[k]) + SQR(v[out_size - k]));
    }
    //cout <<  " " << filter << endl;
  }

  InverseFFT(v);
  v /= (double)ndata;


  for(j=0 ; j <= ndata-1; j++) {
    if(apodize) v[j] /= 1.0 - cos(2.0*M_PI * (double)j/(double)ndata) + 1E-20;
    yout[j] = v[j];
  }
}
 */




void cutdata(double y[], int indm, int indp, double *yout)
{
	//int i;

	for (int i=indm;i<=indp;i++){
		yout[i-indm] = y[i];
	}
}


void cutdata(unsigned char y[], int indm, int indp, unsigned char *yout)
{
	// int i;

	for (int i=indm;i<=indp;i++){
		yout[i-indm] = y[i];
	}
}


void mergedata(double y1[], int ndata1, double y2[], int ndata2, double *yout)
{
	//int i;

	for (int i=0;i<ndata1;i++)
		yout[i] = y1[i];
	for (int i=0;i<ndata2;i++)
		yout[i+ndata1] = y2[i];

}


void dindgen(int nn, double *y)
{
	//int i;

	for (int i=0;i<nn;i++)
		y[i]=(double)i;

}




void fillgaps(double y[], int ndata, double* yout, int* flag, double sign)
{
	// data are assumed to vary linearly in every window

	int i, j, count, countp, countm, seedpass;
	bool sp;

	seedpass = 0;

	//// interval to cut when an event is detected
	//int margcut_m = 15;
	//int margcut_p = 30;
	int margfit   = 20;

	double *xx, *yy, *seriep=NULL, *seriem=NULL, *tempdata1=NULL, *tempdata2=NULL, *xx2;
	double *a;
	double *valtemp;

	a = new double[2];
	a[0]=0.0;
	a[1]=0.0;


	//init random generator
	valtemp = randg(1,0);
	delete [] valtemp;


	////copy data
	for (i=0;i<ndata;i++)
		yout[i] = y[i];

	count = 0;
	sp = 0;
	countm = 0;
	while ((countm<margfit) && (flag[countm] != 0))
		countm++;


	for (i=0;i<ndata;i++){
		if (flag[i] == 1){
			count++;
			sp = 0;
		}
		else{
			sp = 1;
		}

		if (sp && count){
			countp = 0;
			while ((countp < margfit) && (countp+i<ndata-1) && (flag[i+countp] != 1)){
				countp++;
				if (i+countp >= ndata) printf("SDHFIDF\n");
			}

			xx = new double[countp+countm];
			yy = new double[countp+countm];
			xx2 = new double[count];
			fill(xx,xx+(countp+countm),0.0);
			fill(yy,yy+(countp+countm),0.0);
			fill(xx2,xx2+count,0.0);




			if (countm > 0){
				seriem = new double[countm];
				tempdata1 = new double[countm];
				fill(seriem,seriem+countm,0.0);
				fill(tempdata1,tempdata1+countm,0.0);
				dindgen(countm,seriem);
				cutdata(y,i-count-countm,i-count-1,tempdata1);
			}
			if (countp > 0){
				seriep = new double[countp];
				tempdata2 = new double[countp];
				fill(seriep,seriep+countp,0.0);
				fill(tempdata2,tempdata2+countp,0.0);
				dindgen(countp,seriep);
				cutdata(y,i,i+countp-1,tempdata2);
			}



			if (countm && countp){
				mergedata(seriem,countm,seriep,countp,xx);
				mergedata(tempdata1,countm,tempdata2,countp,yy);
			} else {
				if (countm){
					for (j=0;j<countm;j++){
						xx[j] = seriem[j];
						yy[j] = tempdata1[j];
					}
				}
				if (countp){
					for (j=0;j<countp;j++){
						xx[j] = seriep[j];
						yy[j] = tempdata2[j];
					}
				}
			}


			if (countp){
				for (j=0;j<countp;j++){
					xx[countm+j] += double(countm+count);
				}
			}

			dpolyfit(xx,yy,countp+countm,1,a);


			dindgen(count,xx2);
			for (j=0;j<count;j++){
				xx2[j] += double(countm);
			}
			for (j=0;j<count;j++){
				valtemp = randg(1,-1);
				yout[i+j-count] = a[0]+a[1]*xx2[j] + sign*valtemp[0];
				delete [] valtemp;
			}


			if (countp){
				delete[] seriep;
				delete[] tempdata2;
			}
			if (countm){
				delete[] seriem;
				delete[] tempdata1;
			}
			delete[] xx;
			delete[] yy;
			delete[] xx2;

			countm = countp;

			countp = 0;
			count = 0;
			sp = 0;

		}
	}

	delete[] a;

}



double randg_value(long nombre, int seedpass) {

	double nombre_hasard=0.0;
	time_t temps;
	temps = time(NULL);

	unsigned int seed = 0;

	if (seedpass == 0) seed = (unsigned int) temps;
	if (seedpass != 0 && seedpass != -1) seed = (unsigned int) seedpass;
	if (seedpass != -1) srandom(seed);

	//	nombre_hasard= new double[nombre];

	for (long i=0;i<nombre/2;i++) {
		cout << "hmm problem" << endl;
		exit(0);
		//double t1 = (double(rand())/RAND_MAX);
		//double t2 = (double(rand())/RAND_MAX);
		//	nombre_hasard[2*i]=sqrt(-2*log(t1))*cos(2*M_PI*t2);
		///	nombre_hasard[2*i+1]=sqrt(-2*log(t1))*sin(2*M_PI*t2);
	}

	if (nombre/2!=nombre/2.) {
		double t1 = (double(rand())/RAND_MAX);
		double t2 = (double(rand())/RAND_MAX);
		nombre_hasard=sqrt(-2*log(t1))*cos(2*M_PI*t2);
	}//

	//cout << "valtemp" << nombre_hasard << endl;
	return nombre_hasard;
}



double* randg(long nombre, int seedpass) {

	double* nombre_hasard;
	time_t temps;
	temps = time(NULL);

	unsigned int seed = 0;

	if (seedpass == 0) seed = (unsigned int) temps;
	if (seedpass != 0 && seedpass != -1) seed = (unsigned int) seedpass;
	if (seedpass != -1) srandom(seed);

	nombre_hasard= new double[nombre];

	for (long i=0;i<nombre/2;i++) {
		double t1 = (double(rand())/RAND_MAX);
		double t2 = (double(rand())/RAND_MAX);
		nombre_hasard[2*i]=sqrt(-2*log(t1))*cos(2*M_PI*t2);
		nombre_hasard[2*i+1]=sqrt(-2*log(t1))*sin(2*M_PI*t2);
	}

	if (nombre/2!=nombre/2.) {
		double t1 = (double(rand())/RAND_MAX);
		double t2 = (double(rand())/RAND_MAX);
		nombre_hasard[nombre-1]=sqrt(-2*log(t1))*cos(2*M_PI*t2);
	}//


	return nombre_hasard;
}


double* rand(long nombre, int seed) {

	double* nombre_hasard;
	time_t temps;
	temps = time(NULL);
	if (seed == 0) seed = (unsigned int) temps;
	if (seed != -1) srandom(seed);

	nombre_hasard= new double[nombre];

	for (long i=0;i<nombre;i++) nombre_hasard[i]=(double(rand())/RAND_MAX);

	return nombre_hasard;
}


