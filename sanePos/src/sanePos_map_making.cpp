
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>

#include <cstdlib>
#include <cstring>
#include <cstdio> // for printf

#include "sanePos_map_making.h"
#include "dataIO.h"

extern "C" {
#include "nrutil.h"
#include <wcslib/cel.h>
#include <wcslib/wcs.h>
}

using namespace std;


//void computeMapMinima(std::vector<string> bolonames, string *fits_table,
//		long iframe_min, long iframe_max, long *nsamples, double pixdeg,
//		double &ra_min,double &ra_max,double &dec_min,double &dec_max)

void computeMapMinima(std::vector<string> bolonames, struct samples samples_struct,
		long iframe_min, long iframe_max, double pixdeg,
		double &ra_min,double &ra_max,double &dec_min,double &dec_max)
{

	// Compute map extrema by projecting the bolometers offsets back into the sky plane
	// output or update (ra|dec)_(min|max)

	string fits_file;

	long ndet = bolonames.size();


	// Define default values
	ra_min  =  360;
	ra_max  = -360;
	dec_min =  360;
	dec_max = -360;

	for (long iframe=iframe_min;iframe<iframe_max;iframe++){
		// for each scan
		fits_file=samples_struct.fits_table[iframe];

		double *ra, *dec, *phi, **offsets;
		short *flpoint;

		long ns = samples_struct.nsamples[iframe];

		// read bolo offsets
		// TODO : This function should also return the PRJCODE to be used below...
		read_all_bolo_offsets_from_fits(fits_file, bolonames, offsets);

//		for (unsigned long idet = 0; idet < ndet; idet++){
//			cout << offsets[idet][0]*3600 << " " << offsets[idet][1]*3600 << endl;
//		}

		// read reference position
		long test_ns;
		read_ReferencePosition_from_fits(fits_file, ra, dec, phi, flpoint, test_ns);
		if (test_ns != ns) {
			cerr << "Read position does not correspond to frame position" << endl;
			cerr << "Check !!" << endl;
			exit(-1);
		}


		// find the pointing solution at each time stamp for each detector
		struct celprm celestial;
		celini(&celestial);

		// TODO: use the PRJCODE read from the file...
		tanset(&celestial.prj);

		for (long ii=0; ii <ns; ii++){

			celestial.ref[0] =  ra[ii]*15.;
			celestial.ref[1] =  dec[ii];
			celset(&celestial);

			double * offxx, *offyy, *lon, *lat, *ra_deg, *dec_deg;
			int * status;
			offxx   = new double[ndet];
			offyy   = new double[ndet];
			lon     = new double[ndet];
			lat     = new double[ndet];
			ra_deg  = new double[ndet];
			dec_deg = new double[ndet];

			status  = new int[ndet];

			for (long idet=0;idet<ndet;idet++){

				double sinphi = sin(phi[ii]/180.0*M_PI);
				double cosphi = cos(phi[ii]/180.0*M_PI);

				//TODO : check this -1 factor... just a stupid convention...
				offxx[idet] = (cosphi * offsets[idet][0]
				                                      - sinphi * offsets[idet][1])*-1;
				offyy[idet] =  sinphi * offsets[idet][0]
				                                      + cosphi * offsets[idet][1];

			}


			if (celx2s(&celestial, ndet, 1, 0, 1, offxx, offyy, lon, lat, ra_deg, dec_deg, status) == 1) {
				printf("   TAN(X2S) ERROR 1: %s\n", prj_errmsg[1]);
				continue;
			}

			delete [] offxx;;
			delete [] offyy;
			delete [] lon;
			delete [] lat;

			// find coordinates min and max
			double lra_max  = *max_element(ra_deg, ra_deg+ndet);
			double lra_min  = *min_element(ra_deg, ra_deg+ndet);
			double ldec_max = *max_element(dec_deg, dec_deg+ndet);
			double ldec_min = *min_element(dec_deg, dec_deg+ndet);


			if (ra_max < lra_max)    ra_max = lra_max;
			if (ra_min > lra_min)    ra_min = lra_min;
			if (dec_max < ldec_max) dec_max = ldec_max;
			if (dec_min > ldec_min) dec_min = ldec_min;


			delete [] ra_deg;
			delete [] dec_deg;
			delete [] status;

		}



		delete [] ra;
		delete [] dec;
		delete [] phi;
		delete [] flpoint;

		free_dmatrix(offsets,(long)0,ndet-1,(long)0,2-1);
	}

	//TODO : The interval has to be increased or some pixels will be outside the map... NOT UNDERSTOOD WHY...
	// add a small interval of 1 arcmin
	ra_min =  ra_min  - 6.0/60.0/cos((dec_max+dec_min)/2.0/180.0*M_PI);
	ra_max =  ra_max  + 6.0/60.0/cos((dec_max+dec_min)/2.0/180.0*M_PI);
	dec_min = dec_min - 6.0/60.0;
	dec_max = dec_max + 6.0/60.0;

	ra_min  = ra_min/15; // in hour
	ra_max  = ra_max/15;
	dec_min = dec_min;
	dec_max = dec_max;

}

void minmax_flag(double  *& array, short *& flag, long size, double & min_array, double &  max_array){

  // First unflagged data
  long ii=0;
  while(flag[ii] != 0 && ii < size)
    ii++;

  // Start values
  min_array = array[ii];
  max_array = array[ii];


  // Scan the array
  while(ii++ < size-1) 	{
    if (flag[ii] == 0){
      if (array[ii] > max_array) max_array = array[ii];
      if (array[ii] < min_array) min_array = array[ii];
    }
  }
}

void computeMapMinima_HIPE(std::vector<string> bolonames, string *fits_table,
		long iframe_min, long iframe_max, long *nsamples, double pixdeg,
		double &ra_min,double &ra_max,double &dec_min,double &dec_max){

	// Compute map extrema by projecting the bolometers offsets back into the sky plane
	// output (ra|dec)_(min|max)

	string fits_file;
	string field;

	long ndet = bolonames.size();

	double lra_min, lra_max;
	double ldec_min, ldec_max;

	// Define default values
	ra_min  =  360;
	ra_max  = -360;
	dec_min =  360;
	dec_max = -360;

	for (long iframe=iframe_min;iframe<iframe_max;iframe++){
		// for each scan
		fits_file=fits_table[iframe];


		long ns = nsamples[iframe];

		for (long idet=0; idet < ndet; idet++){

			field = bolonames[idet];

			double *ra, *dec;
			short *flag;
			long test_ns;

			read_ra_from_fits(fits_file, field, ra, test_ns);
			if (test_ns != ns) {
				cerr << "Read ra does not correspond to frame size : Check !!" << endl;
				exit(-1);
			}
			read_dec_from_fits(fits_file, field, dec, test_ns);
			if (test_ns != ns) {
				cerr << "Read dec does not correspond to frame size : Check !!" << endl;
				exit(-1);
			}

			read_flag_from_fits(fits_file, field, flag, test_ns);
			if (test_ns != ns) {
				cerr << "Read flag does not correspond to frame size : Check !!" << endl;
				exit(-1);
			}

			//				  for (int ii=0; ii<50; ii++)
			//		    cout << ii << " " << ra[ii] << " " << flag[ii] <<endl;
			//		 exit(0);
			minmax_flag(ra,flag,ns,lra_min,lra_max);
			minmax_flag(dec,flag,ns,ldec_min,ldec_max);

			// cout << field << " " << ns << " " <<ns2 <<endl;
			//			cout << field << " " << lra_max << " " << lra_min << " " << ldec_max << " " << ldec_min << endl;

			if (ra_max < lra_max)    ra_max = lra_max;
			if (ra_min > lra_min)    ra_min = lra_min;
			if (dec_max < ldec_max) dec_max = ldec_max;
			if (dec_min > ldec_min) dec_min = ldec_min;

			//			cout << " field : " << field << " ra_min  : " << lra_min << " ra_max : " << lra_max << " dec_min : " << ldec_min << " dec_max : " << ldec_max << endl;

			delete [] ra;
			delete [] dec;
			delete [] flag;

		}




	}


	/// add a small interval of 10 arcmin
	ra_min =  ra_min  - 1.0/60.0/cos((dec_max+dec_min)/2.0/180.0*M_PI);
	ra_max =  ra_max  + 1.0/60.0/cos((dec_max+dec_min)/2.0/180.0*M_PI);
	dec_min = dec_min - 1.0/60.0;
	dec_max = dec_max + 1.0/60.0;

	ra_min  = ra_min/15; // in hour
	ra_max  = ra_max/15;
	dec_min = dec_min;
	dec_max = dec_max;

}


void computeMapHeader(double pixdeg, char *ctype, char *prjcode, double * coordscorner,
		struct wcsprm &wcs, long &NAXIS1, long &NAXIS2){

	int NAXIS = 2; // image
	int wcsstatus;

	double ra_min = coordscorner[0]*15;
	double ra_max = coordscorner[1]*15;
	double dec_min = coordscorner[2];
	double dec_max = coordscorner[3];

	// Define tangent point
	double ra_mean  = (ra_max+ra_min)/2.0;      // RA in deg
	double dec_mean = (dec_max+dec_min)/2.0;

	// Construct the wcsprm structure
	wcs.flag = -1;
	wcsini(1, NAXIS, &wcs);

	// Pixel size in deg
	for (int ii = 0; ii < NAXIS; ii++) wcs.cdelt[ii] = (ii) ? pixdeg : -1*pixdeg ;
	for (int ii = 0; ii < NAXIS; ii++) strcpy(wcs.cunit[ii], "deg");

	// This will be the reference center of the map
	wcs.crval[0] = ra_mean;
	wcs.crval[1] = dec_mean;

	// Axis label
	if (strcmp(ctype, "EQ") == 0){
		char TYPE[2][5] = { "RA--", "DEC-"};
		char NAME[2][16] = {"Right Ascension","Declination"};

		for (int ii = 0; ii < NAXIS; ii++) {
			strcpy(wcs.ctype[ii], &TYPE[ii][0]);
			strncat(wcs.ctype[ii],"-",1);
			strncat(wcs.ctype[ii],prjcode, 3);
			strcpy(wcs.cname[ii], &NAME[ii][0]);
		}
	}

	if (strcmp(ctype,"GAL") == 0){
		char TYPE[2][5] = { "GLON", "GLAT"};
		char NAME[2][19] = {"Galactic Longitude", "Galactic Latitude"};

		for (int ii = 0; ii < NAXIS; ii++) {
			strcpy(wcs.ctype[ii], &TYPE[ii][0]);
			strncat(wcs.ctype[ii],"-",1);
			strncat(wcs.ctype[ii],prjcode, 3);
			strcpy(wcs.cname[ii], &NAME[ii][0]);
		}
	}

	// Set the structure to have the celestial projection routine in order to ....
	if ((wcsstatus = wcsset(&wcs))) {
		printf("wcsset ERROR %d: %s.\n", wcsstatus, wcs_errmsg[wcsstatus]);
	}


	// .... calculate the size of the map if necessary :
	// As the map could be distorded, deproject the grid into a plane to get the size of the map

	// TODO : instead of testing a gird of nStep x nStep points
	//        	just test the edges
	int nStep = 30;

	double *lon, *lat, *phi, *theta, *x, *y;
	int *status;

	lon = new double[nStep];
	lat = new double[nStep];

	for (int ii=0; ii<nStep; ii++) {
		lon[ii] = (ra_max-ra_min)*ii/(nStep-1)+ra_min;
		lat[ii] = (dec_max-dec_min)*ii/(nStep-1)+dec_min;
	}

	phi    = new double [nStep*nStep];
	theta  = new double [nStep*nStep];
	x      = new double [nStep*nStep];
	y      = new double [nStep*nStep];
	status = new int    [nStep*nStep];

	if (cels2x(&wcs.cel, nStep, nStep, 1, 1, lon, lat, phi, theta, x, y, status) == 1) {
		printf("ERROR 1: %s\n", prj_errmsg[1]);
	}

	// find coordinates min and max
	double x_max  = *max_element(x, x+(nStep*nStep));
	double x_min  = *min_element(x, x+(nStep*nStep));
	double y_max  = *max_element(y, y+(nStep*nStep));
	double y_min  = *min_element(y, y+(nStep*nStep));

	NAXIS1 = ceil((x_max-x_min)/pixdeg);
	NAXIS2 = ceil((y_max-y_min)/pixdeg);

	// Save it as the center of the image

	wcs.crpix[0] = NAXIS1*1./2;
	wcs.crpix[1] = NAXIS2*1./2;

	if ((wcsstatus = wcsset(&wcs))) {
		printf("wcsset ERROR %d: %s.\n", wcsstatus, wcs_errmsg[wcsstatus]);
	}

	//	wcsprt(&wcs);
	delete [] phi;
	delete [] theta;
	delete [] x;
	delete [] y;
	delete [] status;
	delete [] lon;
	delete [] lat;

}

////TODO : remove : not used anymore
//double slaDranrm ( double angle )
///*
// **  - - - - - - - - - -
// **   s l a D r a n r m
// **  - - - - - - - - - -
// **
// **  Normalize angle into range 0-2 pi.
// **
// **  (double precision)
// **
// **  Given:
// **     angle     double      the angle in radians
// **
// **  The result is angle expressed in the range 0-2 pi (double).
// **
// **  Defined in slamac.h:  D2PI, dmod
// **
// **  Last revision:   19 March 1996
// **
// **  Copyright P.T.Wallace.  All rights reserved.
// */
//{
//	double w;
//
//	w = dmod ( angle, D2PI );
//	return ( w >= 0.0 ) ? w : w + D2PI;
//}
//
//
//// TODO: removed : not used anymore
//void slaDs2tp ( double ra, double dec, double raz, double decz,
//		double *xi, double *eta, int *j )
///*
// **  - - - - - - - - -
// **   s l a D s 2 t p
// **  - - - - - - - - -
// **
// **  Projection of spherical coordinates onto tangent plane
// **  ('gnomonic' projection - 'standard coordinates').
// **
// **  (double precision)
// **
// **  Given:
// **     ra,dec      double   spherical coordinates of point to be projected
// **     raz,decz    double   spherical coordinates of tangent point
// **
// **  Returned:
// **     *xi,*eta    double   rectangular coordinates on tangent plane
// **     *j          int      status:   0 = OK, star on tangent plane
// **                                    1 = error, star too far from axis
// **                                    2 = error, antistar on tangent plane
// **                                    3 = error, antistar too far from axis
// **
// **  Last revision:   18 July 1996
// **
// **  Copyright P.T.Wallace.  All rights reserved.
// */
//#define TINY 1e-6
//{
//	double sdecz, sdec, cdecz, cdec, radif, sradif, cradif, denom;
//
//
//	/* Trig functions */
//	sdecz = sin ( decz );
//	sdec = sin ( dec );
//	cdecz = cos ( decz );
//	cdec = cos ( dec );
//	radif = ra - raz;
//	sradif = sin ( radif );
//	cradif = cos ( radif );
//
//	/* Reciprocal of star vector length to tangent plane */
//	denom = sdec * sdecz + cdec * cdecz * cradif;
//
//	/* Handle vectors too far from axis */
//	if ( denom > TINY ) {
//		*j = 0;
//	} else if ( denom >= 0.0 ) {
//		*j = 1;
//		denom = TINY;
//	} else if ( denom > -TINY ) {
//		*j = 2;
//		denom = -TINY;
//	} else {
//		*j = 3;
//	}
//
//	/* Compute tangent plane coordinates (even in dubious cases) */
//	*xi = cdec * sradif / denom;
//	*eta = ( sdec * cdecz - cdec * sdecz * cradif ) / denom;
//}
//
//// TODO: remove : not used anymore
//void slaDtp2s ( double xi, double eta, double raz, double decz,
//		double *ra, double *dec )
///*
// **  - - - - - - - - -
// **   s l a D t p 2 s
// **  - - - - - - - - -
// **
// **  Transform tangent plane coordinates into spherical.
// **
// **  (double precision)
// **
// **  Given:
// **     xi,eta      double   tangent plane rectangular coordinates
// **     raz,decz    double   spherical coordinates of tangent point
// **
// **  Returned:
// **     *ra,*dec    double   spherical coordinates (0-2pi,+/-pi/2)
// **
// **  Called:  slaDranrm
// **
// **  Last revision:   3 June 1995
// **
// **  Copyright P.T.Wallace.  All rights reserved.
// */
//{
//	double sdecz, cdecz, denom;
//
//	sdecz = sin ( decz );
//	cdecz = cos ( decz );
//	denom = cdecz - eta * sdecz;
//	*ra = slaDranrm ( atan2 ( xi, denom ) + raz );
//	*dec = atan2 ( sdecz + eta * cdecz, sqrt ( xi * xi + denom * denom ) );
//}
//
//
//
//
//
//// TODO: remove not used anymore
//void sph_coord_to_sqrmap(double pixdeg, double *ra, double *dec, double *phi,
//		double *offsets, int ns, int *xx, int *yy, int *nn,
//		double *coordscorner, double *tancoord, double *tanpix,
//		bool fixcoord, double radius, /* double *offmap ,*/ double *radecsrc, bool compute_xx_yy)
//{
//	// ra is in hours, dec in degrees
//
//	int j;
//	//int ii, jj, j;
//	double ra_rad, dec_rad, dec_min, dec_max, ra_min, ra_max, dec_mean, ra_mean;
//	double x1, y1, xm, ym, ix, iy, nnf;
//	double sinphi, cosphi, cosphi_phas, sinphi_phas;
//	double *ra_bolo, *dec_bolo, *pixsrc;
//	double offxx, offyy, ra0, dec0, ra00, dec00;
//	double dummy1, dummy2;
//
//
//	double degtorad =  M_PI/180.0;
//
//
//	ra_bolo = new double[ns];
//	dec_bolo = new double[ns];
//	pixsrc = new double[2];
//
//	double dl_rad = pixdeg/180.0*M_PI;
//
//	// find the pointing solution of the detector
//	for (long ii=0;ii<ns;ii++){
//		// Rotation of angle phi[ii] to the offsets
//
//		offxx = cos((phi[ii])/180.0*M_PI)*offsets[0]
//		                                          - sin((phi[ii])/180.0*M_PI)*offsets[1];
//		offyy = sin((phi[ii])/180.0*M_PI)*offsets[0]
//		                                          + cos((phi[ii])/180.0*M_PI)*offsets[1];
//
//		// conversion in pixel scale
//		xm = offxx/pixdeg;
//		ym = offyy/pixdeg;
//
//		// deprojection into spherical coordinates using a -TAN deprojection
//		slaDtp2s ( -xm*dl_rad, ym*dl_rad, ra[ii]/12.0*M_PI, dec[ii]/180.0*M_PI, &ra_rad, &dec_rad );
//
//		// Converting to hour/degree
//		ra_bolo[ii] = ra_rad*12.0/M_PI;
//		dec_bolo[ii] = dec_rad*180.0/M_PI;
//
//	}
//
//	// find coordinates min and max
//	ra_max  = *max_element(ra_bolo, ra_bolo+ns);
//	ra_min  = *min_element(ra_bolo, ra_bolo+ns);
//	dec_max = *max_element(dec_bolo, dec_bolo+ns);
//	dec_min = *min_element(dec_bolo, dec_bolo+ns);
//	//	minmax(ra_bolo, ns, &ra_min, &ra_max, &temp1, &temp2, NULL);
//	//	minmax(dec_bolo, ns, &dec_min, &dec_max, &temp1, &temp2, NULL);
//
//
//
//	/// add a small interval of 2 arcmin
//	ra_min = ra_min - 2.0/60.0/180.0*12.0/cos((dec_max+dec_min)/2.0/180.0*M_PI);
//	ra_max = ra_max + 2.0/60.0/180.0*12.0/cos((dec_max+dec_min)/2.0/180.0*M_PI);
//	dec_min = dec_min - 2.0/60.0;
//	dec_max = dec_max + 2.0/60.0;
//
//
//
//	///////// save or read coordinates of the box
//	if ((coordscorner != NULL) && (fixcoord == 0)){
//		coordscorner[0] = ra_min;
//		coordscorner[1] = ra_max;
//		coordscorner[2] = dec_min;
//		coordscorner[3] = dec_max;
//	}
//
//	if (fixcoord){
//		ra_min = coordscorner[0];
//		ra_max = coordscorner[1];
//		dec_min = coordscorner[2];
//		dec_max = coordscorner[3];
//	}
//
//
//	// define tangent point
//	dec_mean = (dec_max+dec_min)/2.0/180.0 * M_PI;
//	ra_mean = (ra_max+ra_min)/2.0/12.0 * M_PI;
//
//
//
//	// calculate the size of the map if necessary
//	double stepb = 1.0/60;
//	if (radius <= 0){
//		nnf=0.0;
//		for (int ii=int(ra_min/12.0*180.0/stepb);ii<int(ra_max/12.0*180.0/stepb);ii++){
//			for (int jj=int(dec_min/stepb);jj<int(dec_max/stepb);jj++){
//				slaDs2tp ((double)ii/180.0*M_PI*stepb,(double)jj/180.0*M_PI*stepb,ra_mean,dec_mean,&ix   ,&iy  ,&j);
//
//				if (fabs(ix)/dl_rad > nnf/2.0)
//					nnf = (2.0*fabs(ix))/dl_rad;
//				if (fabs(iy)/dl_rad > nnf/2.0)
//					nnf = (2.0*fabs(iy))/dl_rad;
//			}
//		}
//
//		*nn = int(nnf)+1+int(2.0*stepb/pixdeg);
//
//	} else {
//		*nn = int(radius*2.0/pixdeg) + 1;
//	}
//
//
//
//	///// in case telescope coordinates
//	if ((radecsrc != NULL) && (radecsrc[0] >= -100) && (radecsrc[1] >= -100)){
//		ra_rad = radecsrc[0]/12.0 * M_PI;
//		dec_rad = radecsrc[1]/180.0 * M_PI;
//
//		slaDs2tp (ra_rad,dec_rad,ra_mean,dec_mean,&x1,&y1,&j);
//		pixsrc[0] = -x1/dl_rad + double(*nn/2) + 0.5;
//		pixsrc[1] = y1/dl_rad + double(*nn/2) + 0.5;
//
//		// header info
//		tancoord[0] = radecsrc[0] * 15.; // convert to deg
//		tancoord[1] = radecsrc[1];
//		tanpix[0] = (int)pixsrc[0] + 1;// + 0.5;
//		tanpix[1] = (int)pixsrc[1] + 1;// + 0.5;
//
//	}
//
//
//
//	////////////////////////////////////////////////////////////////////
//	// compute x and y for each sample
//	////////////////////////////////////////////////////////////////////
//
//	if(compute_xx_yy){
//		for (long ii=0;ii<ns;ii++){
//
//			ra_rad = ra[ii]/12.0 * M_PI;
//			dec_rad = dec[ii]/180.0 * M_PI;
//
//
//			cosphi = cos(phi[ii]*degtorad);
//			sinphi = sin(phi[ii]*degtorad);
//			cosphi_phas = cos((phi[ii])*degtorad);
//			sinphi_phas = sin((phi[ii])*degtorad);
//
//			offxx = cosphi_phas*offsets[0] - sinphi_phas*offsets[1]; // + offmap[0];
//			offyy = sinphi_phas*offsets[0] + cosphi_phas*offsets[1]; // + offmap[1];
//
//
//
//			slaDs2tp (ra_rad,dec_rad,ra_mean,dec_mean,&x1,&y1,&j);
//			x1 = -x1/dl_rad + double(*nn/2) + 0.5;
//			y1 = y1/dl_rad + double(*nn/2) + 0.5;
//
//
//
//			if ((radecsrc != NULL) && (radecsrc[0] >= -100) && (radecsrc[1] >= -100)){//telescope coordinates
//				xm = pixsrc[0]  + (x1-pixsrc[0]) * cosphi + (y1-pixsrc[1]) * sinphi + offsets[0]/pixdeg;
//				ym = pixsrc[1]  + (y1-pixsrc[1]) * cosphi - (x1-pixsrc[0]) * sinphi + offsets[1]/pixdeg;
//			}else{
//				xm = x1 + offxx/pixdeg;
//				ym = y1 + offyy/pixdeg;
//			}
//
//
//			xx[ii] = int(xm);
//			yy[ii] = int(ym);
//
//		}
//	}
//	//////////////////////////////////////////////////////////////
//
//
//
//
//	////////////////////////////////////////////////
//	//coordinates of the edge pixels
//	//dummy1 = (double(*nn/2)+0.5)*dl_rad;
//	dummy1 = (double(*nn/2))*dl_rad;
//	slaDtp2s (0.0, -dummy1, ra_mean, dec_mean, &ra0, &dec0 );
//	slaDtp2s (0.0, double(*nn)*dl_rad - dummy1, ra_mean, dec_mean, &ra00, &dec00 );
//	slaDtp2s (-dummy1,0.0, ra_mean, dec_mean, &ra0, &dummy2 );
//	slaDtp2s (double(*nn)*dl_rad - dummy1,0.0, ra_mean, dec_mean, &ra00, &dummy2 );
//
//
//	/*ra0 = ra0*12.0/M_PI;
//	dec0 = dec0*180.0/M_PI;
//	ra00 = ra00*12.0/M_PI;
//	dec00 = dec00*180.0/M_PI;*/
//
//	////////////////////////////////////////////////
//
//
//
//	if ((radecsrc == NULL) || (radecsrc[0] < -100) || (radecsrc[1] < -100)){
//		tancoord[0] = ra_mean * 180.0 / M_PI;  // convert to deg
//		tancoord[1] = dec_mean * 180.0 / M_PI; // convert to deg
//		tanpix[0] = (int)(double(*nn)/2.0) + 1;// + 0.5;
//		tanpix[1] = tanpix[0];
//	}
//
//
//
//
//
//	// some cleaning
//	delete[] ra_bolo;
//	delete[] dec_bolo;
//	delete[] pixsrc;
//
//}


/*
void reproj_to_map(double *data, int *xx, int *yy, int ns, double **map, double **count, int nn, short *flag, double **map_f, double **count_f)
{

	//int ii, jj;

	for (int ii=0;ii<nn;ii++){
		for (int jj=0;jj<nn;jj++){
			map[ii][jj] = 0.0;
			map_f[ii][jj] = 0.0;
			count[ii][jj] = 0.0;
			count_f[ii][jj] = 0.0;
		}
	}

	// cerr << "reproj here1?\n";

	for (long ii=0;ii<ns;ii++){
		if (ii == 137909) cerr << ii << ", " << ns << endl;
		if ((flag == NULL) || ((flag[ii] & 1) == 0)){
			//if (ii == 137909) {
	//cerr << xx[ii] << endl;
	//cerr << yy[ii] << endl;
	//cerr << data[ii] << endl;
     // }
			map[xx[ii]][yy[ii]] += data[ii];
			count[xx[ii]][yy[ii]] += 1.0;
		} else{
			map_f[xx[ii]][yy[ii]] += data[ii];
			count_f[xx[ii]][yy[ii]] += 1.0;
		}
	}
	//cerr << "reproj here2?\n";

	for (int ii=0;ii<nn;ii++){
		for (int jj=0;jj<nn;jj++){
			if (count[ii][jj]-0.5 > 0)
				map[ii][jj] = -map[ii][jj]/count[ii][jj];
			if (count_f[ii][jj]-0.5 > 0)
				map_f[ii][jj] = -map_f[ii][jj]/count_f[ii][jj];
		}
	}


}

 */



////TODO: remove : not used anymore
//void flag_conditions(short *flag, /*double *scerr,*/ short *flpoint,
//		long ns, long napod, int *xx, int *yy, int nn, double errarcsec,
//		bool NOFILLGAP, unsigned char *rejectsamp){
//
//
//	// define the rules for bad samples
//
//
//
//	//long ii;
//	short *flagtmp;
//	//double *scerrtmp;
//	short *flpointtmp;
//
//
//	//flagtmp = new unsigned char[ns];
//	//scerrtmp = new double[ns];
//	//flpointtmp = new unsigned char[ns];
//	flagtmp = new short[ns];
//	flpointtmp = new short[ns];
//
//
//	if (NOFILLGAP){
//		for (long ii=0;ii<ns;ii++){
//			flagtmp[ii] = 0;
//			//scerrtmp[ii] = 0.0;
//			flpointtmp[ii] = 0;
//		}
//	} else {
//		for (long ii=0;ii<ns;ii++){
//			flagtmp[ii] = flag[ii];
//			//scerrtmp[ii] = scerr[ii];
//			flpointtmp[ii] = flpoint[ii];
//		}
//	}
//
//
//
//
//
//	for (long ii=0;ii<ns;ii++){
//		rejectsamp[ii] = 0;
//		//if ((flagtmp[ii] == 1) != 0 || (flpointtmp[ii] == 1) != 0) //|| (scerrtmp[ii] > errarcsec))
//		if ((flagtmp[ii] == 1) || (flpointtmp[ii] == 1))
//			rejectsamp[ii] = 1;
//		if ((ii < napod) || (ii >= ns-napod))
//			rejectsamp[ii] = 3;
//	}
//
//
//
//	for (long ii=0;ii<ns;ii++){
//		if ((xx[ii] < 0) || (yy[ii] < 0) || (xx[ii] >= nn) || (yy[ii] >= nn) || (NOFILLGAP && (flag[ii] & 1))){
//			rejectsamp[ii] = 2;
//			if ((xx[ii] < -100000) || (yy[ii] < -100000) || (xx[ii] > 100000) || (yy[ii] > 100000))
//				rejectsamp[ii] = 3;
//		}
//	}
//
//	delete [] flagtmp;
//	//delete [] scerrtmp;
//	delete [] flpointtmp;
//
//}







