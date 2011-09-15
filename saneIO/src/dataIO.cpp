#include <string>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>

#include <cstdlib>
#include <cstring>

#include <math.h>

extern "C" {
#include <nrutil.h>
#include <fitsio.h>
}


#include "dataIO.h"

using namespace std;


int read_all_bolo_offsets_from_fits(string filename, std::vector<string> bolonames, double **& offsets){

	//TODO : Handle angle unit to transform to a common internal known unit

	fitsfile *fptr;
	int status = 0;
	int colnum;
	double * temp_dx,*temp_dy;
	char comment[80];
	int ndet_total;

	long ndet = bolonames.size();
	long indice;

	offsets = dmatrix((long)0,ndet-1,(long)0,2-1);

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// read the Channel List
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "offsets", 0, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	if (fits_read_key(fptr,TLONG, (char *) "NAXIS2", &ndet_total, (char *) &comment, &status)){
		fits_report_error(stderr, status);
		return 1;
	}


	temp_dx = new double[ndet_total];
	temp_dy = new double[ndet_total];

	if(fits_get_colnum(fptr, CASEINSEN, (char *) "dX", &colnum, &status))
		return 1;
	if(fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ndet_total, NULL, temp_dx, 0, &status))
		return 1;

	if(fits_get_colnum(fptr, CASEINSEN, (char *) "dY", &colnum, &status))
		return 1;
	if(fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ndet_total, NULL, temp_dy, 0, &status))
		return 1;

	// match read offsets with requested offsets
	for (long idet=0;idet<ndet;idet++){

		indice = find_channel_index(fptr, (char *) bolonames[idet].c_str());

		// transform arcsec to deg
		offsets[idet][0] = temp_dx[indice-1]/3600;
		offsets[idet][1] = temp_dy[indice-1]/3600;

	}

	delete [] temp_dx;
	delete [] temp_dy;

	if (fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return 1;
	}


	return 0;

}

int read_ReferencePosition_from_fits(string filename, double *&LON, double *&LAT, double *&PHI, long &ns){
	//TODO : Handle angle unit to transform to a common internal known unit

	fitsfile *fptr;
	int status = 0;
	long repeat, width;
	int colnum, typecode;

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// move to the reference position table
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "reference position", 0, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// Retrieve size of the array and ...
	fits_get_num_rows(fptr, &ns, &status);

	// ... allocate corresponding memory
	LON  = new double[ns];
	LAT  = new double[ns];
	PHI  = new double[ns];

	// Read LON
	fits_get_colnum(fptr, CASEINSEN, (char*) "lon", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
	fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ns, NULL, LON, 0, &status);

	// Read LAT
	fits_get_colnum(fptr, CASEINSEN, (char*) "lat", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
	fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ns, NULL, LAT, 0, &status);

	// Read PHI
	fits_get_colnum(fptr, CASEINSEN, (char*) "PHI", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
	fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ns, NULL, PHI, 0, &status);


	if (fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	return 0;

}


int read_flag_from_fits(string filename, string field, int *&mask, long & ns){

	// HIPE like format

	fitsfile *fptr;
	int status = 0;
	int rowIndex = 0, naxis = 0, anynul;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// Retrieve the row of the specified channel
	rowIndex = find_channel_index(fptr, field.c_str());

	// ---------------------------------------------
	// Move ptr to mask hdu
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", 0, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// Retrieve the size of the signal
	if (fits_get_img_dim(fptr, &naxis, &status))
		fits_report_error(stderr, status);
	if(naxis != 2)
		fits_report_error(stderr,BAD_NAXIS);
	if (fits_get_img_size(fptr, 2, naxes, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// Allocate Memory
	ns = naxes[0];
	mask = new int[ns];

	// ---------------------------------------------
	// Retrieve the corresponding row
	fpixel[0] = 1;
	fpixel[1] = rowIndex;
	if (fits_read_pix(fptr, TINT, fpixel, ns, 0, mask, &anynul, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// close file
	if(fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	return 0;
}

int read_signal_from_fits(string filename, string field, double *& signal, long & ns){
	//TODO : Handle unit to transform to a common internal known unit

	// HIPE like format

	fitsfile *fptr;
	int status = 0;
	int rowIndex = 0, naxis = 0, anynul;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// Retrieve the row of the specified channel
	rowIndex = find_channel_index(fptr, field.c_str());

	// ---------------------------------------------
	// Move ptr to signal hdu
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", 0, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// Retrieve the size of the signal
	if (fits_get_img_dim(fptr, &naxis, &status)){
		fits_report_error(stderr, status);
		return 1;
	}
	if(naxis != 2){
		fits_report_error(stderr, status);
		return 1;
	}
	if (fits_get_img_size(fptr, 2, naxes, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// Allocate Memory
	ns = naxes[0];
	signal = new double[ns];

	// ---------------------------------------------
	// Retrieve the corresponding row
	fpixel[0] = 1;
	fpixel[1] = rowIndex;
	if (fits_read_pix(fptr, TDOUBLE, fpixel, ns, 0, signal, &anynul, &status)){
		fits_report_error(stderr, status);
		return 1;
	}


	// ---------------------------------------------
	// close file
	if(fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	return 0;

}

int read_bolo_list(string fname,std::vector<string> &boloname, long &ndet){

	fitsfile *fptr;
	int status = 0;
	char **temp_bolo;

	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	read_channels(fptr,temp_bolo, ndet);

	for (long ii=0; ii < ndet; ii++)
		boloname.push_back(temp_bolo[ii]);

	// close file
	if(fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	for (long ii=0; ii < ndet; ii++)
		delete [] temp_bolo[ii];
	delete [] temp_bolo;

	return 0;
}

int read_channels(fitsfile *fptr, char **& data, long &nBolos){
	//TODO : Handle angle unit to transform to a common internal known unit

	int status = 0;
	long repeat, width;
	int colnum, typecode;

	// ---------------------------------------------
	// read the Channel List
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "channels", 0, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	fits_get_num_rows(fptr, &nBolos, &status);
	fits_get_colnum(fptr, CASEINSEN, (char*) "names", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);

	// Initialize the data container
	data = new char*[nBolos];
	for (int i = 0; i < nBolos; i++) {
		data[i] = new char[repeat];
	}

	fits_read_col(fptr, TSTRING, colnum, 1, 1, nBolos, NULL, data, 0, &status); // TBYTE pour PACS ??

	return 0;
}

long find_channel_index(fitsfile *fptr, const char * field){
	long nBolos;
	char **data;

	// read the channel list
	read_channels(fptr,data,nBolos);

	//find the fits index of the bolo and return it
	//fits index are 1 indexed, so indexes starts at 1
	for (long idet = 1; idet <= nBolos; idet++)
		if (strcmp(field,data[idet-1]) == 0){
			for(long ii=0;ii<nBolos;ii++)
				delete [] data[ii];
			delete [] data;
			return idet;
		}

	cout << "EE - " << field << " not found" << endl;
	for(long ii=0;ii<nBolos;ii++)
		delete [] data[ii];
	delete [] data;
	return -1L;

}


int read_LON_LAT_from_fits(string filename, string field, double *&lon, double *& lat, long & ns){
	//TODO : Handle angle unit to transform to a common internal known unit

	// HIPE like format

	fitsfile *fptr;
	int status = 0;
	int rowIndex = 0, naxis = 0, anynul;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// Retrieve the row of the specified channel
	rowIndex = find_channel_index(fptr, field.c_str());

	// ---------------------------------------------
	// Move ptr to signal hdu
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "lon", 0, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// Retrieve the size of the signal
	if (fits_get_img_dim(fptr, &naxis, &status)){
		fits_report_error(stderr, status);
		return 1;
	}
	if(naxis != 2){
		fits_report_error(stderr, status);
		return 1;
	}
	if (fits_get_img_size(fptr, 2, naxes, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// Allocate Memory
	ns = naxes[0];
	lon = new double[ns];

	// ---------------------------------------------
	// Retrieve the corresponding row
	fpixel[0] = 1;
	fpixel[1] = rowIndex;
	if (fits_read_pix(fptr, TDOUBLE, fpixel, ns, 0, lon, &anynul, &status)){
		fits_report_error(stderr, status);
		return 1;
	}



	// ---------------------------------------------
	// Move ptr to signal hdu
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "lat", 0, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// Retrieve the size of the signal
	if (fits_get_img_dim(fptr, &naxis, &status)){
		fits_report_error(stderr, status);
		return 1;
	}
	if(naxis != 2){
		fits_report_error(stderr, status);
		return 1;
	}
	if (fits_get_img_size(fptr, 2, naxes, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// Allocate Memory
	ns = naxes[0];
	lat = new double[ns];

	// ---------------------------------------------
	// Retrieve the corresponding row
	fpixel[0] = 1;
	fpixel[1] = rowIndex;
	if (fits_read_pix(fptr, TDOUBLE, fpixel, ns, 0, lat, &anynul, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// close file
	if(fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	return 0;
}

int read_time_from_fits(string filename, double *& time, long ns){

	// HIPE like format

	fitsfile *fptr;
	int status = 0;
	int ns_test = 0;
	char comment[80];

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// Move ptr to signal hdu
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "time", 0, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	if (fits_read_key(fptr,TLONG, (char *) "NAXIS1", &ns_test, (char *) &comment, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	if(ns!=ns_test){
		cout << "time image has a wrong size : " << ns_test << " != " << ns << endl;
		return 1;
	}

	// ---------------------------------------------
	// Allocate Memory
	time = new double[ns];

	// ---------------------------------------------
	// Retrieve the corresponding col
	if(fits_read_col(fptr, TDOUBLE, 2, 1, 1, ns_test, NULL, time, 0, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	// ---------------------------------------------
	// close file
	if(fits_close_file(fptr, &status)){
		fits_report_error(stderr, status);
		return 1;
	}

	return 0;

}


int test_format(string fitsname){

	fitsfile *fptr;
	int format=1; // 1 = HIPE, 2 = sanepic
	int status = 0;

	if (fits_open_file(&fptr, fitsname.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	if ((fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "lon", 0, &status)>0) &&
			(fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "lat", 0, &status)>0)) // "lon" and "lat" tables were not found
		format=2;

	status = 0;

	if((fits_movnam_hdu(fptr, BINARY_TBL, (char*) "reference position", 0, &status)>0) &&
			(fits_movnam_hdu(fptr, BINARY_TBL, (char*) "offsets", 0, &status)>0)){ // both tables were not found
		if(format==2)
			format=0;
		else
			format=1;
	}

	status = 0;

	// close file
	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

	return format;
}

void copy_offsets(fitsfile * fptr, fitsfile *outfptr)
/*! copy offsets table from this file to output file */
{

	int status=0; // fits error status

	fits_movnam_hdu(fptr, BINARY_TBL, (char*) "offsets", 0, &status); // move HDU pointer to offsets table
	fits_copy_header(fptr, outfptr, &status); // copy header informations

	for(int col=1;col<4;col++) // copy the 3 columns to output file
		fits_copy_col(fptr, outfptr,  col, col,	0, &status);

}


void copy_channels(fitsfile * fptr, fitsfile *outfptr)
/*! copy channels list from this file to output file */
{

	int status=0; // fits error status


	fits_movnam_hdu(fptr, BINARY_TBL, (char*) "channels", 0, &status); // move HDU pointer to channels table
	fits_copy_header(fptr, outfptr, &status); // copy header information


	fits_copy_col(fptr, outfptr,  1, 1,	0, &status); // copy detectors list

}
