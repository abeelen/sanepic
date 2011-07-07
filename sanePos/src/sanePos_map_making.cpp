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
#include "inputFileIO.h"
#include "temporary_IO.h"
#include "todprocess.h"


extern "C" {
#include "nrutil.h"
#include <wcslib/cel.h>
#include <wcslib/wcs.h>
#include <wcslib/sph.h>
#include <wcslib/wcsmath.h>
#include <wcslib/wcstrig.h>
#include <wcslib/prj.h>
}

using namespace std;

int computeMapMinima(struct samples samples_struct, string dirfile,
		long iframe_min, long iframe_max,
		double &ra_min,double &ra_max,double &dec_min,double &dec_max, std::vector<std::vector<std::string> > bolo_vect)
{

	// Compute map extrema by projecting the bolometers offsets back into the sky plane
	// output or update (ra|dec)_(min|max)

	string fits_file;
	string field; // test


	// Define default values
	ra_min  =  360;
	ra_max  = -360;
	dec_min =  360;
	dec_max = -360;

	for (long iframe=iframe_min;iframe<iframe_max;iframe++){
		// for each scan
		fits_file=dirfile + samples_struct.fitsvect[iframe];

		std::vector<string> det_vect = bolo_vect[iframe];

		long ndet = (long)det_vect.size();

		double *ra, *dec, *phi, **offsets;

		long ns = samples_struct.nsamples[iframe];

		// read bolo offsets
		// TODO : This function should also return the PRJCODE to be used below...
		if(read_all_bolo_offsets_from_fits(fits_file, det_vect, offsets))
			return 1;

		// read reference position
		long test_ns;
		if(read_ReferencePosition_from_fits(fits_file, ra, dec, phi, test_ns))
			return 1;

		if (test_ns != ns) {
			cerr << "Read position does not correspond to frame position" << endl;
			cerr << "Check !!" << endl;
			return 1;
		}

		// find the pointing solution at each time stamp for each detector
		struct celprm celestial;
		celini(&celestial);

		// TODO: use the PRJCODE read from the file...
		tanset(&celestial.prj);

		for (long ii=0; ii <ns; ii++){

			celestial.ref[0] =  ra[ii]*15.0;
			celestial.ref[1] =  dec[ii];
			if(celset(&celestial))
				cout << "problem celset\n";

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


			if (celx2s(&celestial, ndet, 0, 1, 1, offxx, offyy, lon, lat, ra_deg, dec_deg, status) == 1) {
				printf("   TAN(X2S) ERROR 1: %s\n", prj_errmsg[1]);
				continue;
			}


			delete [] offxx;
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

		free_dmatrix(offsets,(long)0,ndet-1,(long)0,2-1);
	}

	// The interval has to be increased or some pixels will be outside the map...
	// add a small interval of 1 arcmin
	ra_min =  ra_min  - 6.0/60.0/cos((dec_max+dec_min)/2.0/180.0*M_PI);
	ra_max =  ra_max  + 6.0/60.0/cos((dec_max+dec_min)/2.0/180.0*M_PI);
	dec_min = dec_min - 6.0/60.0;
	dec_max = dec_max + 6.0/60.0;

	/// add a small interval of 2 arcmin
	//	ra_min = ra_min - 2.0/60.0/180.0*12.0/cos((dec_max+dec_min)/2.0/180.0*M_PI);
	//	ra_max = ra_max + 2.0/60.0/180.0*12.0/cos((dec_max+dec_min)/2.0/180.0*M_PI);
	//	dec_min = dec_min - 2.0/60.0;
	//	dec_max = dec_max + 2.0/60.0;

	ra_min  = ra_min/15; // in hour
	ra_max  = ra_max/15;
	dec_min = dec_min;
	dec_max = dec_max;

	return 0;

}

int minmax_flag(double  *& array, int *& flag, long size, double & min_array, double &  max_array){

	// get First unflagged data indice

	long ii=0;
	while(flag[ii] != 0 && ii < size)
		ii++;

	// Everything is flagged : Run saneCheck/Fix !
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

int computeMapMinima_HIPE(std::string tmp_dir, struct samples samples_struct,
		long iframe_min, long iframe_max,
		double &ra_min,double &ra_max,double &dec_min,double &dec_max, std::vector<std::vector<std::string> > bolo_vect){

	// Compute map extrema by projecting the bolometers position
	// output (ra|dec)_(min|max)

	//TODO: Wrong scheme.... :
	// - Shall project with the correct projection + fake center
	// - Find proper center
	// - Reproject all data
	// - find map extrema in projected space shall return the WCS keywords

	string base_file;
	string field;
	int drop_sanepos=0;



	double lra_min, lra_max;
	double ldec_min, ldec_max;

	// Define default values
	ra_min  =  360;
	ra_max  = -360;
	dec_min =  360;
	dec_max = -360;

	for (long iframe=iframe_min;iframe<iframe_max;iframe++){
		// for each scan
		base_file=samples_struct.basevect[iframe];

		std::vector<string> det_vect = bolo_vect[iframe];

		long ndet = (long)det_vect.size();

		long ns = samples_struct.nsamples[iframe];

		for (long idet=0; idet < ndet; idet++){

			field = det_vect[idet];

			double *ra, *dec;
			int *flag=NULL;

			//			if(read_ra_dec_from_fits(fits_file, field, ra, dec, test_ns))
			if(read_RA_from_dirfile(samples_struct.dirfile_pointer, base_file, field, ra, ns))
				return 1;
			if(read_DEC_from_dirfile(samples_struct.dirfile_pointer, base_file, field, dec, ns))
				return 1;

			if(read_flag_from_dirfile(samples_struct.dirfile_pointer, base_file, field, flag, ns))
				return 1;

			if( minmax_flag(ra,flag,ns,lra_min,lra_max) ||
					minmax_flag(dec,flag,ns,ldec_min,ldec_max) ){

				cerr << "WARNING - frame : " << base_file << " : " << field << " has no usable data : Check !!" << endl;
				drop_sanepos++;

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


	if(drop_sanepos>0)
		return 1;


	/// add a small interval of 10 arcmin
	ra_min =  ra_min  - 1.0/60.0/cos((dec_max+dec_min)/2.0/180.0*M_PI);
	ra_max =  ra_max  + 1.0/60.0/cos((dec_max+dec_min)/2.0/180.0*M_PI);
	dec_min = dec_min - 1.0/60.0;
	dec_max = dec_max + 1.0/60.0;

	ra_min  = ra_min/15; // in hour
	ra_max  = ra_max/15;
	dec_min = dec_min;
	dec_max = dec_max;

	return 0;

}


void computeMapHeader(double pixdeg, char *ctype, char *prjcode, double * coordscorner,
		struct wcsprm * &wcs, long &NAXIS1, long &NAXIS2){

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
	wcs = (struct wcsprm *) malloc(sizeof(struct wcsprm));
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

	// TODO : instead of testing a grid of nStep x nStep points
	//        	just test the edges
	int nStep = 30;

	double *lon, *lat, *phi, *theta, *x, *y;
	int *status;

	lon = new double[nStep];
	lat = new double[nStep];

	//TODO: This scheme will fail over the Pole,
	//one should use/define a reference RA/DEC, project the data, and compute sizes on the projected intermediate coordinates.

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




int do_PtNd_Naiv(struct samples samples_struct, double *PNd, std::string outdir, std::vector<std::string> det, long ndet, int orderpoly, int napod, double f_lppix, long ns,
		long long *indpix, long iframe, long *hits)
{


	string field1;
	double *data, *data_lp, *data_out, *bfilter;
	long long *samptopix;
	double aa, bb;
	int *flag;

	data_lp =  new double[ns];
	data_out = new double[ns];
	samptopix = new long long[ns];
	bfilter = new double[ns/2+1];


	for (long idet1=0;idet1<ndet;idet1++){

		field1 = det[idet1];

		//Read pointing data
		if(read_samptopix(samples_struct.dirfile_pointer, ns, samptopix, samples_struct.basevect[iframe], field1))
			return 1;

		if(read_data_from_dirfile(samples_struct.dirfile_pointer, samples_struct.basevect[iframe], field1, data, ns))
			return 1;
		if(read_flag_from_dirfile(samples_struct.dirfile_pointer, samples_struct.basevect[iframe], field1, flag, ns))
			return 1;

		fill(data_out,data_out+ns,0.0);

		//fill gaps with straight line
		fillgaps2(data,ns,data_out,flag,40);
		for (long ii=0;ii<ns;ii++)
			data[ii] = data_out[ii];

		//remove polynomia to correct from time varying calibration
		remove_poly(data,ns,orderpoly,data_out,flag);
		for (long ii=0;ii<ns;ii++)
			data[ii] = data_out[ii];

		//linear prediction
		for (long ii=0;ii<ns;ii++)
			data_lp[ii] = data[ii];


		/// remove a baseline
		aa = (data_lp[ns-1]-data[0])/double(ns);
		bb = data_lp[0];
		for (long ii=0;ii<ns;ii++)
			data_lp[ii] -= aa*(double)ii+bb;


		if (f_lppix > 0.0)
			butterworth(data_lp,ns,f_lppix,8,data_out,bfilter,1,napod,0);


		fill(data,data+ns,0.0);
		//******************* process gaps
		fillgaps2(data_out,ns,data,flag,40);


		for (long ii=0;ii<ns;ii++){
			PNd[indpix[samptopix[ii]]] += data[ii];
		}

		//compute hit counts
		for (long ii=0;ii<ns;ii++){
			hits[indpix[samptopix[ii]]] += 1;
		}

		delete[] data;
		delete[] flag;

	}// end of idet1 loop


	delete[] samptopix;
	delete[] data_out;
	delete[] data_lp;
	delete[] bfilter;


	return 0;
}


