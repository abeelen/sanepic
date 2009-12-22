
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
//		short *flpoint;

		long ns = samples_struct.nsamples[iframe];

		// read bolo offsets
		// TODO : This function should also return the PRJCODE to be used below...
		read_all_bolo_offsets_from_fits(fits_file, bolonames, offsets);

		//		for (unsigned long idet = 0; idet < ndet; idet++){
		//			cout << offsets[idet][0]*3600 << " " << offsets[idet][1]*3600 << endl;
		//		}

		// read reference position
		long test_ns;
		read_ReferencePosition_from_fits(fits_file, ra, dec, phi, test_ns);
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
//		delete [] flpoint;

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

int minmax_flag(double  *& array, short *& flag, long size, double & min_array, double &  max_array){

	// First unflagged data

	long ii=0;
	while(flag[ii] != 0 && ii < size)
		ii++;

	// Everything is flagged
	if (ii == size)
		return EXIT_FAILURE;

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

	return EXIT_SUCCESS;
}

void computeMapMinima_HIPE(std::vector<string> bolonames, struct samples samples_struct,
		long iframe_min, long iframe_max, double pixdeg,
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
		fits_file=samples_struct.fits_table[iframe];


		long ns = samples_struct.nsamples[iframe];

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
			//
			//			for (int ii=0; ii<50; ii++)
			//				cout << ii << " " << ra[ii] << " " << flag[ii] <<endl;
			//			exit(0);


			if( minmax_flag(ra,flag,ns,lra_min,lra_max) ||
					minmax_flag(dec,flag,ns,ldec_min,ldec_max) ){

				cerr << "WW - " << field << " has no usable data : Check !!" << endl;

			} else {

				if (ra_max < lra_max)    ra_max = lra_max;
				if (ra_min > lra_min)    ra_min = lra_min;
				if (dec_max < ldec_max) dec_max = ldec_max;
				if (dec_min > ldec_min) dec_min = ldec_min;
			}

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
		struct wcsprm *& wcs, long &NAXIS1, long &NAXIS2){

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
	wcs->flag = -1;
	wcsini(1, NAXIS, wcs);

	// Pixel size in deg
	for (int ii = 0; ii < NAXIS; ii++) wcs->cdelt[ii] = (ii) ? pixdeg : -1*pixdeg ;
	for (int ii = 0; ii < NAXIS; ii++) strcpy(wcs->cunit[ii], "deg");

	// This will be the reference center of the map
	wcs->crval[0] = ra_mean;
	wcs->crval[1] = dec_mean;

	// Axis label
	if (strcmp(ctype, "EQ") == 0){
		char TYPE[2][5] = { "RA--", "DEC-"};
		char NAME[2][16] = {"Right Ascension","Declination"};

		for (int ii = 0; ii < NAXIS; ii++) {
			strcpy(wcs->ctype[ii], &TYPE[ii][0]);
			strncat(wcs->ctype[ii],"-",1);
			strncat(wcs->ctype[ii],prjcode, 3);
			strcpy(wcs->cname[ii], &NAME[ii][0]);
		}
	}

	if (strcmp(ctype,"GAL") == 0){
		char TYPE[2][5] = { "GLON", "GLAT"};
		char NAME[2][19] = {"Galactic Longitude", "Galactic Latitude"};

		for (int ii = 0; ii < NAXIS; ii++) {
			strcpy(wcs->ctype[ii], &TYPE[ii][0]);
			strncat(wcs->ctype[ii],"-",1);
			strncat(wcs->ctype[ii],prjcode, 3);
			strcpy(wcs->cname[ii], &NAME[ii][0]);
		}
	}

	// Set the structure to have the celestial projection routine in order to ....
	if ((wcsstatus = wcsset(wcs))) {
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

	if (cels2x(&(wcs->cel), nStep, nStep, 1, 1, lon, lat, phi, theta, x, y, status) == 1) {
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

	wcs->crpix[0] = NAXIS1*1./2;
	wcs->crpix[1] = NAXIS2*1./2;

	if ((wcsstatus = wcsset(wcs))) {
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

