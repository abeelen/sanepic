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
		long iframe_min, long iframe_max, struct wcsprm * & wcs,
		double &lon_min,double &lon_max,double &lat_min,double &lat_max, std::vector<std::vector<std::string> > bolo_vect)
{
	// Compute map extrema by projecting the bolometers offsets back into the sky plane
	// output or update (lon|lat) (min|max)

	string fits_file;
	string field; // test


	// Define default values
	lon_min  =  360;
	lon_max  = -360;
	lat_min =  360;
	lat_max = -360;

	for (long iframe=iframe_min;iframe<iframe_max;iframe++){
		// for each scan
		fits_file=dirfile + samples_struct.fitsvect[iframe];

		std::vector<string> det_vect = bolo_vect[iframe];

		long ndet = (long)det_vect.size();

		double *lon, *lat, *phi, **offsets;

		long ns = samples_struct.nsamples[iframe];

		// read bolo offsets
		// TODO : This function should also return the PRJCODE to be used below...
		if(read_all_bolo_offsets_from_fits(fits_file, det_vect, offsets))
			return 1;

		// read reference position
		long test_ns;
		if(read_ReferencePosition_from_fits(fits_file, lon, lat, phi, test_ns))
			return 1;

		if (test_ns != ns) {
			cerr << "Read position does not correspond to frame position" << endl;
			cerr << "Check !!" << endl;
			return 1;
		}

		// find the pointing solution at each time stamp for each detector
		struct celprm arrayProj;
		celini(&arrayProj);

		// TODO: use the PRJCODE read from the file...
		tanset(&arrayProj.prj);

		for (long ii=0; ii <ns; ii++){

			arrayProj.ref[0] =  lon[ii];
			arrayProj.ref[1] =  lat[ii];

			if(celset(&arrayProj))
				cout << "problem celset\n";

			double * offxx, *offyy, *lon_deg, *lat_deg;
			int * status;

			double *dummy1, *dummy2, *x, *y;

			offxx   = new double[ndet];
			offyy   = new double[ndet];
			lon_deg = new double[ndet];
			lat_deg = new double[ndet];


			dummy1  = new double [ndet];
			dummy2  = new double [ndet];
			x       = new double [ndet];
			y       = new double [ndet];


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

			// Projection away from the detector plane, into the spherical sky...
			if (celx2s(&arrayProj, ndet, 0, 1, 1, offxx, offyy, dummy1, dummy2, lon_deg, lat_deg, status) == 1) {
				printf("   TAN(X2S) ERROR 1: %s\n", prj_errmsg[1]);
				continue;
			}

			// Projection back into the sky plane ...
			if (cels2x(&(wcs->cel), ndet, 0, 1, 1, lon_deg, lat_deg, dummy1, dummy2, x, y, status) == 1) {
				printf("ERROR 1: %s\n", prj_errmsg[1]);
			}

			// find coordinates min and max
			double llon_max  = *max_element(x, x+ndet);
			double llon_min  = *min_element(x, x+ndet);
			double llat_max = *max_element(y, y+ndet);
			double llat_min = *min_element(y, y+ndet);

			if (lon_max < llon_max)    lon_max = llon_max;
			if (lon_min > llon_min)    lon_min = llon_min;
			if (lat_max < llat_max) lat_max = llat_max;
			if (lat_min > llat_min) lat_min = llat_min;



			delete [] offxx;
			delete [] offyy;
			delete [] x;
			delete [] y;
			delete [] lon_deg;
			delete [] lat_deg;
			delete [] dummy1;
			delete [] dummy2;
			delete [] status;

		}



		delete [] lon;
		delete [] lat;
		delete [] phi;

		free_dmatrix(offsets,(long)0,ndet-1,(long)0,2-1);
	}

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
		long iframe_min, long iframe_max, struct wcsprm * & wcs,
		double &lon_min,double &lon_max,double &lat_min,double &lat_max, std::vector<std::vector<std::string> > bolo_vect){

	// Compute map extrema by projecting the bolometers position
	// output (lon|lat)_(min|max)

	string base_file;
	string field;
	int drop_sanepos=0;

	double llon_min, llon_max;
	double llat_min, llat_max;

	// Define default values
	lon_min  =  360;
	lon_max  = -360;
	lat_min =  360;
	lat_max = -360;

	for (long iframe=iframe_min;iframe<iframe_max;iframe++){
		// for each scan
		base_file=samples_struct.basevect[iframe];

		std::vector<string> det_vect = bolo_vect[iframe];

		long ndet = (long)det_vect.size();

		long ns = samples_struct.nsamples[iframe];

		double *phi, *theta, *x, *y;
		int *status;

		phi    = new double [ns];
		theta  = new double [ns];
		x      = new double [ns];
		y      = new double [ns];
		status = new int    [ns];

		for (long idet=0; idet < ndet; idet++){

			field = det_vect[idet];

			double *lon, *lat;
			int *flag=NULL;


			//			if(read_lon_lat_from_fits(fits_file, field, lon, lat, test_ns))
			if(read_LON_from_dirfile(samples_struct.dirfile_pointer, base_file, field, lon, ns))
				return 1;
			if(read_LAT_from_dirfile(samples_struct.dirfile_pointer, base_file, field, lat, ns))
				return 1;

			if(read_flag_from_dirfile(samples_struct.dirfile_pointer, base_file, field, flag, ns))
				return 1;


			if (cels2x(&(wcs->cel), ns, 0, 1, 1, lon, lat, phi, theta, x, y, status) == 1) {
				printf("ERROR 1: %s\n", prj_errmsg[1]);
			}

			if( minmax_flag(x,flag,ns,llon_min,llon_max) ||
					minmax_flag(y,flag,ns,llat_min,llat_max) ){

				cerr << "WARNING - frame : " << base_file << " : " << field << " has no usable data : Check !!" << endl;
				drop_sanepos++;

			} else {

				if (lon_max < llon_max)    lon_max = llon_max;
				if (lon_min > llon_min)    lon_min = llon_min;
				if (lat_max < llat_max) lat_max = llat_max;
				if (lat_min > llat_min) lat_min = llat_min;
			}

			delete [] lon;
			delete [] lat;
			delete [] flag;

		}

		delete [] phi;
		delete [] theta;
		delete [] x;
		delete [] y;

	}

	if(drop_sanepos>0)
		return 1;

	return 0;

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


