#include <iostream>

#include <math.h>
//#include "/fir_data/patanch/numrec/inc/nrutil.h"
#include "todprocess.h"
#include "map_making.h"
#include <fftw3.h>
#include <time.h>

#define NR_END 1
#define FREE_ARG char*

using namespace std;


long *data_compare;


double slaDranrm ( double angle )
/*
**  - - - - - - - - - -
**   s l a D r a n r m
**  - - - - - - - - - -
**
**  Normalize angle into range 0-2 pi.
**
**  (double precision)
**
**  Given:
**     angle     double      the angle in radians
**
**  The result is angle expressed in the range 0-2 pi (double).
**
**  Defined in slamac.h:  D2PI, dmod
**
**  Last revision:   19 March 1996
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
   double w;

   w = dmod ( angle, D2PI );
   return ( w >= 0.0 ) ? w : w + D2PI;
}



void slaDs2tp ( double ra, double dec, double raz, double decz,
                double *xi, double *eta, int *j )
/*
**  - - - - - - - - -
**   s l a D s 2 t p
**  - - - - - - - - -
**
**  Projection of spherical coordinates onto tangent plane
**  ('gnomonic' projection - 'standard coordinates').
**
**  (double precision)
**
**  Given:
**     ra,dec      double   spherical coordinates of point to be projected
**     raz,decz    double   spherical coordinates of tangent point
**
**  Returned:
**     *xi,*eta    double   rectangular coordinates on tangent plane
**     *j          int      status:   0 = OK, star on tangent plane
**                                    1 = error, star too far from axis
**                                    2 = error, antistar on tangent plane
**                                    3 = error, antistar too far from axis
**
**  Last revision:   18 July 1996
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
#define TINY 1e-6
{
   double sdecz, sdec, cdecz, cdec, radif, sradif, cradif, denom;


/* Trig functions */
   sdecz = sin ( decz );
   sdec = sin ( dec );
   cdecz = cos ( decz );
   cdec = cos ( dec );
   radif = ra - raz;
   sradif = sin ( radif );
   cradif = cos ( radif );

/* Reciprocal of star vector length to tangent plane */
   denom = sdec * sdecz + cdec * cdecz * cradif;

/* Handle vectors too far from axis */
   if ( denom > TINY ) {
      *j = 0;
   } else if ( denom >= 0.0 ) {
      *j = 1;
      denom = TINY;
   } else if ( denom > -TINY ) {
      *j = 2;
      denom = -TINY;
   } else {
      *j = 3;
   }

   /* Compute tangent plane coordinates (even in dubious cases) */
   *xi = cdec * sradif / denom;
   *eta = ( sdec * cdecz - cdec * sdecz * cradif ) / denom;
}


void slaDtp2s ( double xi, double eta, double raz, double decz,
                double *ra, double *dec )
/*
**  - - - - - - - - -
**   s l a D t p 2 s
**  - - - - - - - - -
**
**  Transform tangent plane coordinates into spherical.
**
**  (double precision)
**
**  Given:
**     xi,eta      double   tangent plane rectangular coordinates
**     raz,decz    double   spherical coordinates of tangent point
**
**  Returned:
**     *ra,*dec    double   spherical coordinates (0-2pi,+/-pi/2)
**
**  Called:  slaDranrm
**
**  Last revision:   3 June 1995
**
**  Copyright P.T.Wallace.  All rights reserved.
*/
{
  double sdecz, cdecz, denom;

  sdecz = sin ( decz );
  cdecz = cos ( decz );
  denom = cdecz - eta * sdecz;
  *ra = slaDranrm ( atan2 ( xi, denom ) + raz );
  *dec = atan2 ( sdecz + eta * cdecz, sqrt ( xi * xi + denom * denom ) );
}






void sph_coord_to_sqrmap(double pixdeg, double *ra, double *dec, double *phi,
			 double *offsets, int ns, int *xx, int *yy, int *nn,
			 double *coordscorner, double *tancoord, double *tanpix,
			 bool fixcoord, double radius, double *offmap, double *radecsrc)
{
  // ra is in hours, dec in degrees

  int ii, jj, j;
  double ra_rad, dec_rad, dec_min, dec_max, ra_min, ra_max, dec_mean, ra_mean;
  double x1, y1, xm, ym, ix, iy, nnf;
  double sinphi, cosphi, cosphi_phas, sinphi_phas;
  double *ra_bolo, *dec_bolo, *pixsrc;
  double offxx, offyy, ra0, dec0, ra00, dec00;
  int temp1, temp2;
  double dummy1, dummy2;


  double degtorad =  M_PI/180.0;


  ra_bolo = new double[ns];
  dec_bolo = new double[ns];
  pixsrc = new double[2];

  double ang_sc2array = 0.0;

  double dl_rad = pixdeg/180.0*M_PI;




  // find the pointing solution of the detector
  for (ii=0;ii<ns;ii++){
    offxx = cos((phi[ii]-ang_sc2array)/180.0*M_PI)*offsets[0]
      - sin((phi[ii]-ang_sc2array)/180.0*M_PI)*offsets[1] + offmap[0];
    offyy = sin((phi[ii]-ang_sc2array)/180.0*M_PI)*offsets[0]
      + cos((phi[ii]-ang_sc2array)/180.0*M_PI)*offsets[1] + offmap[1];

    xm = offxx/pixdeg;
    ym = offyy/pixdeg;

    slaDtp2s ( -xm*dl_rad, ym*dl_rad, ra[ii]/12.0*M_PI, dec[ii]/180.0*M_PI, &ra_rad, &dec_rad );
    ra_bolo[ii] = ra_rad*12.0/M_PI;
    dec_bolo[ii] = dec_rad*180.0/M_PI;
  }



  // find coordinates min and max
  minmax(ra_bolo, ns, &ra_min, &ra_max, &temp1, &temp2, NULL);
  minmax(dec_bolo, ns, &dec_min, &dec_max, &temp1, &temp2, NULL);



  /// add a small interval of 2 arcmin
  ra_min = ra_min - 2.0/60.0/180.0*12.0/cos((dec_max+dec_min)/2.0/180.0*M_PI);
  ra_max = ra_max + 2.0/60.0/180.0*12.0/cos((dec_max+dec_min)/2.0/180.0*M_PI);
  dec_min = dec_min - 2.0/60.0;
  dec_max = dec_max + 2.0/60.0;



  ///////// save or read coordinates of the box
  if ((coordscorner != NULL) && (fixcoord == 0)){
    coordscorner[0] = ra_min;
    coordscorner[1] = ra_max;
    coordscorner[2] = dec_min;
    coordscorner[3] = dec_max;
  }

  if (fixcoord){
    ra_min = coordscorner[0];
    ra_max = coordscorner[1];
    dec_min = coordscorner[2];
    dec_max = coordscorner[3];
  }


  // define tangent point
  dec_mean = (dec_max+dec_min)/2.0/180.0 * M_PI;
  ra_mean = (ra_max+ra_min)/2.0/12.0 * M_PI;



  // calculate the size of the map if necessary
  double stepb = 1.0/60;
  if (radius <= 0){
    nnf=0.0;
    for (ii=int(ra_min/12.0*180.0/stepb);ii<int(ra_max/12.0*180.0/stepb);ii++){
      for (jj=int(dec_min/stepb);jj<int(dec_max/stepb);jj++){
	slaDs2tp ((double)ii/180.0*M_PI*stepb,(double)jj/180.0*M_PI*stepb,ra_mean,dec_mean,&ix   ,&iy  ,&j);

	if (fabs(ix)/dl_rad > nnf/2.0)
	  nnf = (2.0*fabs(ix))/dl_rad;
	if (fabs(iy)/dl_rad > nnf/2.0)
	  nnf = (2.0*fabs(iy))/dl_rad;
      }
    }
    *nn = int(nnf)+1+int(2.0*stepb/pixdeg);

  } else {
    *nn = int(radius*2.0/pixdeg) + 1;
  }



  ///// in case telescope coordinates
  if ((radecsrc != NULL) && (radecsrc[0] >= -100) && (radecsrc[1] >= -100)){
    ra_rad = radecsrc[0]/12.0 * M_PI;
    dec_rad = radecsrc[1]/180.0 * M_PI;

    slaDs2tp (ra_rad,dec_rad,ra_mean,dec_mean,&x1,&y1,&j);
    pixsrc[0] = -x1/dl_rad + double(*nn/2) + 0.5;
    pixsrc[1] = y1/dl_rad + double(*nn/2) + 0.5;

    // header info
    tancoord[0] = radecsrc[0] * 15.; // convert to deg
    tancoord[1] = radecsrc[1];
    tanpix[0] = (int)pixsrc[0] + 1;// + 0.5;
    tanpix[1] = (int)pixsrc[1] + 1;// + 0.5;

  }



  ////////////////////////////////////////////////////////////////////
  // compute x and y for each sample
  ////////////////////////////////////////////////////////////////////


 for (ii=0;ii<ns;ii++){

    ra_rad = ra[ii]/12.0 * M_PI;
    dec_rad = dec[ii]/180.0 * M_PI;


    cosphi = cos(phi[ii]*degtorad);
    sinphi = sin(phi[ii]*degtorad);
    cosphi_phas = cos((phi[ii]-ang_sc2array)*degtorad);
    sinphi_phas = sin((phi[ii]-ang_sc2array)*degtorad);

    offxx = cosphi_phas*offsets[0] - sinphi_phas*offsets[1] + offmap[0];
    offyy = sinphi_phas*offsets[0] + cosphi_phas*offsets[1] + offmap[1];



    slaDs2tp (ra_rad,dec_rad,ra_mean,dec_mean,&x1,&y1,&j);
    x1 = -x1/dl_rad + double(*nn/2) + 0.5;
    y1 = y1/dl_rad + double(*nn/2) + 0.5;



    if ((radecsrc != NULL) && (radecsrc[0] >= -100) && (radecsrc[1] >= -100)){//telescope coordinates
      xm = pixsrc[0]  + (x1-pixsrc[0]) * cosphi + (y1-pixsrc[1]) * sinphi + offsets[0]/pixdeg;
      ym = pixsrc[1]  + (y1-pixsrc[1]) * cosphi - (x1-pixsrc[0]) * sinphi + offsets[1]/pixdeg;
    }else{
      xm = x1 + offxx/pixdeg;
      ym = y1 + offyy/pixdeg;
    }


    xx[ii] = int(xm);
    yy[ii] = int(ym);

  }
  //////////////////////////////////////////////////////////////




  ////////////////////////////////////////////////
  //coordinates of the edge pixels
  //dummy1 = (double(*nn/2)+0.5)*dl_rad;
  dummy1 = (double(*nn/2))*dl_rad;
  slaDtp2s (0.0, -dummy1, ra_mean, dec_mean, &ra0, &dec0 );
  slaDtp2s (0.0, double(*nn)*dl_rad - dummy1, ra_mean, dec_mean, &ra00, &dec00 );
  slaDtp2s (-dummy1,0.0, ra_mean, dec_mean, &ra0, &dummy2 );
  slaDtp2s (double(*nn)*dl_rad - dummy1,0.0, ra_mean, dec_mean, &ra00, &dummy2 );


  ra0 = ra0*12.0/M_PI;
  dec0 = dec0*180.0/M_PI;
  ra00 = ra00*12.0/M_PI;
  dec00 = dec00*180.0/M_PI;

  ////////////////////////////////////////////////



  if ((radecsrc == NULL) || (radecsrc[0] < -100) || (radecsrc[1] < -100)){
    tancoord[0] = ra_mean * 180.0 / M_PI;  // convert to deg
    tancoord[1] = dec_mean * 180.0 / M_PI; // convert to deg
    tanpix[0] = (int)(double(*nn)/2.0) + 1;// + 0.5;
    tanpix[1] = tanpix[0];
  }





  // some cleaning
  delete[] ra_bolo;
  delete[] dec_bolo;
  delete[] pixsrc;

}









void reproj_to_map(double *data, int *xx, int *yy, int ns, double **map, double **count, int nn, unsigned char *flag, double **map_f, double **count_f)
{

  int ii, jj;

  for (ii=0;ii<nn;ii++){
    for (jj=0;jj<nn;jj++){
      map[ii][jj] = 0.0;
      map_f[ii][jj] = 0.0;
      count[ii][jj] = 0.0;
      count_f[ii][jj] = 0.0;
    }
  }

  // cerr << "reproj here1?\n";

  for (ii=0;ii<ns;ii++){
    if (ii == 137909) cerr << ii << ", " << ns << endl;
    if ((flag == NULL) || ((flag[ii] & 1) == 0)){
      /*if (ii == 137909) {
	cerr << xx[ii] << endl;
	cerr << yy[ii] << endl;
	cerr << data[ii] << endl;
      }*/
      map[xx[ii]][yy[ii]] += data[ii];
      count[xx[ii]][yy[ii]] += 1.0;
    } else{
      map_f[xx[ii]][yy[ii]] += data[ii];
      count_f[xx[ii]][yy[ii]] += 1.0;
    }
  }
  //cerr << "reproj here2?\n";

  for (ii=0;ii<nn;ii++){
    for (jj=0;jj<nn;jj++){
      if (count[ii][jj]-0.5 > 0)
	map[ii][jj] = -map[ii][jj]/count[ii][jj];
      if (count_f[ii][jj]-0.5 > 0)
	map_f[ii][jj] = -map_f[ii][jj]/count_f[ii][jj];
    }
  }


}







void compute_PtNmd(double *data, double *Nk, long ndata, long marge, int nn,
		   long *indpix, long *samptopix, int npix, double *PNd){

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




  for (ii=-marge;ii<ndata-marge;ii++){
    if ((ii < 0) || (ii >= ndata-2*marge)){
      PNd[npix-2] += Nd[ii+marge];
    } else {
      PNd[indpix[samptopix[ii]]] += Nd[ii+marge];
    }
  }




  delete [] Ndf;
  delete [] fdata;
  delete [] Nd;


}







void compute_PtNmd_corr(double *data, double *Nk, unsigned char *rejectsamp, unsigned char *binsamp,
		   long ndata, long marge, int *xx, int *yy, int nn,
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



  for (ii=-marge;ii<ndata-marge;ii++){
    if ((ii < 0) || (ii >= ndata-2*marge)){
      PNd[npix-2] += Nd[ii+marge];
    } else {
      if (rejectsamp[ii] == 0){
	if (binsamp[ii] == 1){
	  PNd[npix-2] += Nd[ii+marge];
	} else {
	  ll = indpix[yy[ii]*nn + xx[ii]];
	  PNd[ll] += Nd[ii+marge];
	}
      } else {
	PNd[npix-1] += Nd[ii+marge];
      }
    }
  }

  delete [] Ndf;
  delete [] fdata;
  delete [] Nd;


}






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




  for (ii=-marge;ii<ndata-marge;ii++){
    if ((ii < 0) || (ii >= ndata-2*marge)){
      PNd[npix-2] += Nd[ii+marge];
    } else {
      if (rejectsamp[ii] == 0){
	if (binsamp[ii] == 1){
	  PNd[npix-2] += Nd[ii+marge];
	} else {
	  ll = indpix[yy[ii]*nn + xx[ii]];
	  PNd[ll] += Nd[ii+marge];
	}
      } else {
	PNd[npix-1] += Nd[ii+marge];
      }
    }
  }

  delete [] Ndf;
  delete [] Nd;


}








void compute_PtNP(double *Nk, unsigned char *rejectsamp, unsigned char *binsamp, long ndata,
		  long marge, int *xx, int *yy, int nn, long *indpix,
		  int npix, double f_lppix, double *PtNP){


  long ii, k, jj, kk, ll, ll2, indPtNP;
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
  for (k=0;k<ndata/2+1;k++){
	Nk_[k][0] = 1.0/Nk[k]/(double)ndata/(double)ndata;
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
	pixpos[ii+marge] = npix-1;
      }
    }
  }


  for (ii=0;ii<npix;ii++)
    jj_sqr[ii] = ii*(ii+1)/2;



  for (ii=0;ii<ndata;ii++){
    ll = pixpos[ii];
    ll2 = ll*(ll+1)/2;
    if (ii){
      for (kk=MAX(ii-(ndata)/MAX(2,int(f_lppix+0.5)),0);kk<ii;kk++){
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







void compute_PtNP_frac(double *Nk, unsigned char *rejectsamp, unsigned char *binsamp, long ndata,
		  long marge, int *xx, int *yy, int nn, long *indpix,
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
	pixpos[ii+marge] = npix-1;
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








void compute_diagPtNP(double *Nk, long *samptopix, long ndata,
		      long marge, int nn, long *indpix,
		      int npix, double f_lppix, double *dPtNP){


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
	Nk_[k][0] = 1.0/abs(Nk[k])/(double)ndata/(double)ndata;
	Nk_[k][1] = 0.0;
  }
  fftplan = fftw_plan_dft_c2r_1d(ndata, Nk_, N_, FFTW_ESTIMATE);
  fftw_execute(fftplan);





  for (ii=-marge;ii<ndata-marge;ii++){
    if ((ii < 0) || (ii >= ndata-2*marge)){
      pixpos[ii+marge] = npix-2;
    } else {
      pixpos[ii+marge] = indpix[samptopix[ii]];
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
	if ((ipix == npix-1) || (ipix == npix-2)){ //This is just to avoid spending to much time computing this pixel which could contain a lot of data
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






void compute_diagPtNPCorr(double *Nk, long *samptopix, long ndata,
			  long marge, int nn, long *indpix,
			  int npix, double f_lppix, double *dPtNP){


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
      pixpos[ii+marge] = indpix[samptopix[ii]];
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
	if ((ipix == npix-1) || (ipix == npix-2)){ //This is just to avoid spending to much time computing this pixel
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









void flag_conditions(unsigned char *flag, double *scerr, unsigned char *flpoint,
		     long ns, long napod, long marge, int *xx, int *yy, int nn, double errarcsec,
		     bool NOFILLGAP, unsigned char *rejectsamp){


  // define the rules for bad samples



  long ii;
  unsigned char *flagtmp;
  double *scerrtmp;
  unsigned char *flpointtmp;

  flagtmp = new unsigned char[ns];
  scerrtmp = new double[ns];
  flpointtmp = new unsigned char[ns];


  if (NOFILLGAP){
    for (ii=0;ii<ns;ii++){
      flagtmp[ii] = 0;
      scerrtmp[ii] = 0.0;
      flpointtmp[ii] = 0;
    }
  } else {
    for (ii=0;ii<ns;ii++){
      flagtmp[ii] = flag[ii];
      scerrtmp[ii] = scerr[ii];
      flpointtmp[ii] = flpoint[ii];
    }
  }





  for (ii=0;ii<ns;ii++){
    rejectsamp[ii] = 0;
    if ((flagtmp[ii] & 1) != 0 || (scerrtmp[ii] > errarcsec) || (flpointtmp[ii] & 1) != 0)
      rejectsamp[ii] = 1;
    if ((ii < napod-marge) || (ii >= ns-napod+marge))
      rejectsamp[ii] = 3;
  }



  for (ii=0;ii<ns;ii++){
    if ((xx[ii] < 0) || (yy[ii] < 0) || (xx[ii] >= nn) || (yy[ii] >= nn) || (NOFILLGAP && (flag[ii] & 1))){
      rejectsamp[ii] = 2;
      if ((xx[ii] < -100000) || (yy[ii] < -100000) || (xx[ii] > 100000) || (yy[ii] > 100000))
	rejectsamp[ii] = 3;
    }
  }

  delete [] flagtmp;
  delete [] scerrtmp;
  delete [] flpointtmp;

}









void MapMakPreProcessData(double *data, unsigned char *flag, double *calp, long ns, int marge, int napod,
			  int orderpoly, double f_lppix, double *data_lp, double *bfilter, bool NORMLIN, bool NOFILLGAP, double *Ps){


  long ii;
  double aa, bb;

  double *data_out, *data_out_lp;

  data_out = new double[ns];
  data_out_lp = new double[ns+2*marge];


  //*********************************************************************


  //  cout << "fill gap" << endl;
  if (NOFILLGAP == 0){
    //fill 5Bgaps with straight line
    fillgaps(data,ns,data_out,flag,0);
    for (ii=0;ii<ns;ii++)
      data[ii] = data_out[ii];
  }


  //  cout << "remove polynom" << endl;
  //remove polynomia
  remove_poly(data,ns,orderpoly,data_out,0);


  //  cout << "varying calibration " << endl;
  //correct from time varying calibration
  for (ii=0;ii<ns;ii++)
    data[ii] = data_out[ii]*calp[ii/20];



  //  cout << "linear prediction" << endl;
  //linear prediction
  for (ii=0;ii<ns;ii++)
    data_lp[ii+marge] = data[ii];
  if (marge) Pad(data_lp,marge,ns,ns+2*marge);



  //  cout << "baselining" << endl;
  if (NORMLIN == 0){
    /// remove a baseline
    aa = (data_lp[ns+2*marge-1]-data[0])/double(ns+2*marge);
    bb = data_lp[0];
    for (ii=0;ii<ns+2*marge;ii++)
      data_lp[ii] -= aa*(double)ii+bb;
  }



  //  cout << "butterworth filter" << endl;
  //Butterworth filter (if necessary)
  if (f_lppix > 0.0){
    butterworth(data_lp,ns+2*marge,f_lppix,8,data_out_lp,bfilter,1,napod,0);
   for (ii=0;ii<(ns+2*marge);ii++)
      data_lp[ii] = data_out_lp[ii];
  } else{
    for (ii=0;ii<(ns+2*marge)/2+1;ii++)
      bfilter[ii] = 1.0;
  }



  if (Ps != NULL)
    for (ii=marge;ii<ns+marge;ii++)
      data_lp[ii] = data_lp[ii] - Ps[ii];

  //  cout << "fill gap" << endl;
  //******************* process gaps
  if (NOFILLGAP == 0){
    for (ii=0;ii<ns;ii++)
      data_out[ii] = data_lp[ii+marge];
    fillgaps(data_out,ns,data,flag,0);
    for (ii=0;ii<ns;ii++)
      data_lp[ii+marge] = data[ii];
  }


  if (Ps != NULL){
    //    cout << "linear predidction " << endl;
    //linear prediction
    if (marge) Pad(data_lp,marge,ns,ns+2*marge);

    for (ii=marge;ii<ns+marge;ii++)
      //if (flag[ii-marge] == 0)
	data_lp[ii] = data_lp[ii] + Ps[ii];
  }



  delete [] data_out;
  delete [] data_out_lp;



}







void noisepectrum_estim(double *data, int ns, double *ell, int nbins, double fsamp, double *bfilter, double *Nell, double *Nk){

  int ii,k,q;
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
  for (ii=0;ii<ns;ii++)
    datatemp[ii] = datatemp2[ii]*apodwind[ii];

  //Fourier transform the data
  fftplan = fftw_plan_dft_r2c_1d(ns, datatemp, fdata, FFTW_ESTIMATE);
  fftw_execute(fftplan);




  totapod = 0.0;
  for (ii=0;ii<ns;ii++)
    totapod += apodwind[ii]*apodwind[ii];


  //power spectrum
  for (k=0;k<ns/2+1;k++){
    Nk[k] = pow(fdata[k][0],2) + pow(fdata[k][1],2);
    Nk[k] = Nk[k]/(totapod/(double)ns)/(double)ns;
  }



  //bin power spectrum
  for (q=0;q<nbins;q++){
    Nell[q] = 0.0;
    count[q] = 0;
  }


  q=0;
  for (k=0;k<ns/2+1;k++){
    if (k >= ell[q+1]/fsamp*(double)ns)
      if (q<nbins-1)
	q++;
    Nell[q] += Nk[k];
    count[q] += 1;
  }


  for (q=0;q<nbins;q++)
    Nell[q] /= double(count[q]);



  for(ii=0;ii<ns/2+1;ii++){
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

  int ii,k,q;

  double *bfiltertemp;
  int *count;


  bfiltertemp = new double[ns/2+1];
  count = new int[nbins];


  //power spectrum
  for (k=0;k<ns/2+1;k++){
    Nk[k] = fdata1[k][0]*fdata2[k][0] + fdata1[k][1]*fdata2[k][1];
    Nk[k] = Nk[k]/(double)ns;
  }


  //bin power spectrum
  for (q=0;q<nbins;q++){
    Nell[q] = 0.0;
    count[q] = 0;
  }


  q=0;
  for (k=0;k<ns/2+1;k++){
    if (k >= ell[q+1]/fsamp*(double)ns)
      if (q<nbins-1)
	q++;
    Nell[q] += Nk[k];
    count[q] += 1;
  }


  for (q=0;q<nbins;q++)
    Nell[q] /= double(count[q]);



  for(ii=0;ii<ns/2+1;ii++){
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








void readNSpectrum(char *nameSpfile, double *bfilter, long ns, long marge, double fsamp, double *Nk){

  FILE *fp;

  int nbins;
  long ii;
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

  for (ii=0;ii<nbins;ii++){
    fscanf(fp,"%lf %lf",&dummy1,&dummy2);
    ell[ii] = dummy1;
    SpN[ii] = dummy2;
  }
  fscanf(fp,"%lf",&dummy1);
  ell[nbins] = dummy1;
  fclose(fp);


  // interpolate logarithmically the noise power spectrum
  binnedSpectrum2log_interpol(ell,SpN,bfilter,nbins,ns+2*marge,fsamp,Nk);



  delete[] SpN;
  delete[] ell;



}






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


void correctFrameOffsets(int nfoff, long ff, double *offsets, foffset *foffsets, double *froffsets){

  int find ;

  for (find=0; find<nfoff-1; find++) {
    if (ff < (foffsets[find+1]).frame) break;
  }

  froffsets[0] = offsets[0] - (foffsets[find]).yaw;
  froffsets[1] = offsets[1] + (foffsets[find]).pitch;

}




void deproject(double *S, long *indpix, long *samptopix, long ndata, long marge, long nn, long npix, double *Ps, int flgdupl, int factdupl, long ntotscan, long *indpsrc, long npixsrc){

  long ii, ll, iframe;
  double a, b;





  for (ii=-marge;ii<ndata-marge;ii++){
    if ((ii < 0) || (ii >= ndata-2*marge)){
      Ps[ii+marge] = S[npix-2];
    } else {
      Ps[ii+marge] = S[indpix[samptopix[ii]]];

      // in case we replaced flagged data
      if ((flgdupl == 2) && (samptopix[ii] >= nn*nn) && (samptopix[ii] < 2*nn*nn) && (indpix[samptopix[ii] - nn*nn] >= 0)){
	if (indpix[samptopix[ii] - nn*nn] >= 0){
	  Ps[ii+marge] = S[indpix[samptopix[ii]-nn*nn]];
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
	    Ps[ii+marge] = a/b;
	}
      }
    }
  }






      //      if (rejectsamp[ii] == 0){
      //	if (binsamp[ii] == 1){
      //	  Ps[ii+marge] = S[npix-2];
      //	} else {
      //	  ll = indpix[yy[ii]*nn + xx[ii]];
      //	  Ps[ii+marge] = S[ll];
      //	}
      //      } else {
      //	if (flgdupl == 1){
      //	  ll = indpix[(yy[ii]*nn + xx[ii]) + nn*nn];
      //	  Ps[ii+marge] = S[ll];
      //	}
      //	if (flgdupl == 0) {
      //	  Ps[ii+marge] = S[npix-1];
      //	}
      //	if (flgdupl == 2) {
      //	  ll = indpix[(yy[ii]*nn + xx[ii])];
      //	  if (ll >= 0) Ps[ii+marge] = S[ll];
      //	  if (ll < 0) Ps[ii+marge] = S[indpix[(yy[ii]*nn + xx[ii]) + nn*nn]];
      //	}
      //      }
      //    }

      //  }

 }





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







void deproject_new(double *S, long *indpix, int *xx, int *yy, unsigned char *rejectsamp, unsigned char *binsamp, long ndata, long marge, long nn, long npix, long npixmap, double *Ps, long *countreject){

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







int compare_global_array_long (const void *a, const void *b)
{

  const long *da = (const long *) a;
  const long *db = (const long *) b;

  return (data_compare[*da] > data_compare[*db]) - (data_compare[*da] < data_compare[*db]);
}

