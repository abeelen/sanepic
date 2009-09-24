#ifndef MAP_MAKING
#define MAP_MAKING

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <fftw3.h>

//#include "/fir_data/patanch/numrec/inc/nrutil.h"
//#include <matpack/matpack.h>

#define D2PI 6.2831853071795864769252867665590057683943387987502

/* dmod(A,B) - A modulo B (double) */
#define dmod(A,B) ((B)!=0.0?((A)*(B)>0.0?(A)-(B)*floor((A)/(B))\
                                        :(A)+(B)*floor(-(A)/(B))):(A))

void compute_PtNmd(double *data, double *Nk, long ndata, int nn,
		   long *indpix, long *samptopix, int npix, double *PNd);


/*void compute_PtNmd_corr(double *data, double *Nk, unsigned char *rejectsamp, unsigned char *binsamp,
			long ndata, int *xx, int *yy, int nn,
			long *indpix, int npix, double *PNd);*/

/*void compute_PtNmfftd_corr(fftw_complex *fdata, double *Nk, unsigned char *rejectsamp, unsigned char *binsamp,
		   long ndata, long marge, int *xx, int *yy, int nn,
			   long *indpix, int npix, double *PNd);*/
/*
void compute_PtNP(double *Nk, unsigned char *rejectsamp, unsigned char *binsamp, long ndata,
		  int *xx, int *yy, int nn, long *indpix,
		  int npix, double f_lppix, double *PtNP);*/
/*
void compute_PtNP_frac(double *Nk, unsigned char *rejectsamp, unsigned char *binsamp, long ndata,
		   int *xx, int *yy, int nn, long *indpix,
		       int npix, double f_lppix, double *PtNP, int nfrac, int ifrac);*/

/*
void compute_PtNP_corr(double *Nk, unsigned char *rejectsamp1, unsigned char *rejectsamp2,
		       unsigned char *binsamp1, unsigned char *binsamp2,
		       long ndata, int *xx1, int *yy1, int *xx2, int *yy2,
		       int nn, long *indpix, int npix, double f_lppix, double *PtNP);*/



void compute_diagPtNP(double *Nk, long *samptopix, long ndata,
		      int nn, long *indpix,
		      int npix, double f_lppix, double *dPtNP);



void compute_diagPtNPCorr(double *Nk, long *samptopix, long ndata,
			  int nn, long *indpix,
			  int npix, double f_lppix, double *dPtNP);


/*
void compute_diagPtNPCorr_msk(double *Nk, unsigned char *mask, long iframe,
			      unsigned char *rejectsamp, unsigned char *binsamp,
			      long ndata, long marge, int *xx, int *yy, int nn,
			      long *indpix, int npix, double f_lppix, double *dPtNP);


void compute_diagPtNPCorr_new(double *Nk, unsigned char *rejectsamp,
			      unsigned char *binsamp, long ndata, long marge,
			      int *xx, int *yy, int nn, long *indpix, int npix, int npixmap,
			      double f_lppix, double *dPtNP, long *countreject);*/




void MapMakPreProcessData(double *data,  short *flag, /*double *calp,*/ long ns, int napod,
			  int orderpoly, double f_lppix, double *data_lp, double *bfilter, bool NORMLIN,
			  bool NOFILLGAP,bool remove_polynomia, double *Ps = NULL);


void noisepectrum_estim(double *data, int ns, double *ell, int nbins, double fsamp,
			double *bfilter, double *Nell, double *Nk);


//void noisecrosspectrum_estim(double *data1, double *data2, int ns, double *ell, int nbins, double fsamp, double *bfilter, double *Nell, double *Nk, bool APOD, bool APPOLY);

void noisecrosspectrum_estim(fftw_complex *fdata1, fftw_complex *fdata2, int ns, double *ell, int nbins, double fsamp, double *bfilter, double *Nell, double *Nk);

//void noisepectrum_estim(double *data, int ns, double *ell, int nbins, double fsamp, double *bfilter,
//		double *Nell, double *Nk);



//void readNSpectrum(char *nameSpfile, double *bfilter, long ns, double fsamp, double *Nk);
void readNSpectrum(string nameSpfile, double *bfilter, long ns, double fsamp, double *Nk);


/*
void readalldata(long ff, long ns, string field, string ra_field, string dec_field,
		 string phi_field, string scerr_field, string flpoint_field,
		 string dirfile, string bextension, string fextension, string cextension,
		 double *data, double *calp, double *ra, double *dec,
		 double *phi, double *scerr, unsigned char *flpoint,
		 unsigned char *flag, int shift_data_to_point);*/





void deproject(double *S, long *indpix, long *samptopix, long ndata, long nn,
	       long npix, double *Ps, int flgdupl, int factdupl, long ntotscan = 0,
	       long *indpsrc = NULL, long npixsrc = 0);


//void deproject_msk(double *S, unsigned char *mask, long *indpix, int *xx, int *yy, unsigned char *rejectsamp, unsigned char *binsamp, long ndata, long marge, long nn, long npix, long iframe, double *Ps);


//void deproject_new(double *S, long *indpix, int *xx, int *yy, unsigned char *rejectsamp, unsigned char *binsamp, long ndata, long marge, long nn, long npix, long npixmap, double *Ps, long *countreject);



int compare_global_array_long (const void *array_1, const void *array_2);


#endif
