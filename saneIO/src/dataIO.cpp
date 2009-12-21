
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

// #include "getdata.h"

using namespace std;


void read_all_bolo_offsets_from_fits(string filename, std::vector<string> bolonames, double **& offsets){

	//TODO : Handle angle unit to transform to a common internal known unit
	fitsfile *fptr;
	int status = 0;
	int colnum;
	double * temp;
	temp = new double[2];

	long ndet = bolonames.size();

	offsets = dmatrix((long)0,ndet-1,(long)0,2-1);

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// read the Channel List
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "Channels Offsets", NULL, &status))
		fits_report_error(stderr, status);

	// match read offsets with requested offsets
	for (long idet=0;idet<ndet;idet++){


		fits_get_colnum(fptr, CASEINSEN, (char *) bolonames[idet].c_str(), &colnum, &status);
		fits_read_col(fptr, TDOUBLE, colnum, 1, 1, 2, NULL, temp, 0, &status);

		// transform arcsec to deg
		for (int ii=0; ii<2; ii++)	offsets[idet][ii] = temp[ii]/3600;

	}

	delete [] temp;

	if (fits_close_file(fptr, &status))
			fits_report_error(stderr, status);


}

void read_ReferencePosition_from_fits(string filename, double *&RA, double *&DEC, double *&PHI, short *&FLAG, long &ns){
	//TODO : Handle angle unit to transform to a common internal known unit

	fitsfile *fptr;
	int status = 0;
	long repeat, width;
	int colnum, typecode;

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// move to the reference position table
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "Reference Position", NULL, &status))
		fits_report_error(stderr, status);

	// Retrieve size of the array and ...
	fits_get_num_rows(fptr, &ns, &status);

	// ... allocate corresponding memory
	RA   = new double[ns];
	DEC  = new double[ns];
	PHI  = new double[ns];
	FLAG = new short[ns];

	// Read RA
	fits_get_colnum(fptr, CASEINSEN, (char*) "RA", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
	fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ns, NULL, RA, 0, &status);

	//TODO: Remove this step
	// transform RA in hours
	for(long ii = 0; ii<ns; ii++)
		RA[ii]=RA[ii]/15;

	// Read DEC
	fits_get_colnum(fptr, CASEINSEN, (char*) "DEC", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
	fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ns, NULL, DEC, 0, &status);

	// Read PHI
	fits_get_colnum(fptr, CASEINSEN, (char*) "PHI", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
	fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ns, NULL, PHI, 0, &status);

	// TODO: Does not exist anymore...
	// Read FLAG
	fits_get_colnum(fptr, CASEINSEN, (char*) "FLAG", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
	fits_read_col(fptr, TSHORT, colnum, 1, 1, ns, NULL, FLAG, 0, &status);

//	// check for bad values
//	for (long ii = 0; ii<ns; ii++)
//		if (isnan(RA[ii]) || isnan(DEC[ii]))
//			FLAG[ii] = 1;

	if (fits_close_file(fptr, &status))
			fits_report_error(stderr, status);

}


void read_flpoint_from_fits(string filename, short *FLAG){

	fitsfile *fptr;
	int status = 0;
	int colnum;
	long ns;
	//double temp;

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// read the Channel List
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "Reference Position", NULL, &status))
		fits_report_error(stderr, status);

	fits_get_num_rows(fptr, &ns, &status);
	fits_get_colnum(fptr, CASEINSEN, (char*) "FLAG", &colnum, &status);
	fits_read_col(fptr, TSHORT, colnum, 1, 1, ns, NULL, FLAG, 0, &status);

	// close file
	if(fits_close_file(fptr, &status)) // ajout mat 15/09
		fits_report_error(stderr, status);

}


void read_flag_from_fits(string filename, string field, short *& mask, long & ns){

	// HIPE like format

	fitsfile *fptr;
	int status = 0;
	int rowIndex = 0, naxis = 0, anynul;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// Retrieve the row of the specified channel
	rowIndex = find_channel_index(fptr, field.c_str());

	// ---------------------------------------------
	// Move ptr to signal hdu
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", NULL, &status))
		fits_report_error(stderr, status);

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
	mask = new short[ns];

	// ---------------------------------------------
	// Retrieve the corresponding row
	fpixel[0] = 1;
	fpixel[1] = rowIndex;
	if (fits_read_pix(fptr, TSHORT, fpixel, ns, 0, mask, &anynul, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// close file
	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

}



void read_signal_from_fits(string filename, string field, double *& signal, long & ns){
	//TODO : Handle unit to transform to a common internal known unit

	// HIPE like format

	fitsfile *fptr;
	int status = 0;
	int rowIndex = 0, naxis = 0, anynul;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// Retrieve the row of the specified channel
	rowIndex = find_channel_index(fptr, field.c_str());

	// ---------------------------------------------
	// Move ptr to signal hdu
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", NULL, &status))
		fits_report_error(stderr, status);

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
	signal = new double[ns];

	// ---------------------------------------------
	// Retrieve the corresponding row
	fpixel[0] = 1;
	fpixel[1] = rowIndex;
	if (fits_read_pix(fptr, TDOUBLE, fpixel, ns, 0, signal, &anynul, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// close file
	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

}


void read_channels(fitsfile *fptr, char **& data, long &nBolos){
	//TODO : Handle angle unit to transform to a common internal known unit

	int status = 0;
	long repeat, width;
	int colnum, typecode;

	// ---------------------------------------------
	// read the Channel List
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "channels", NULL, &status))
		fits_report_error(stderr, status);

	fits_get_num_rows(fptr, &nBolos, &status);
	fits_get_colnum(fptr, CASEINSEN, (char*) "names", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);

	// Initialize the data container
	data = new char*[nBolos];
	for (int i = 0; i < nBolos; i++) {
		data[i] = new char[repeat];
	}

	fits_read_col(fptr, TSTRING, colnum, 1, 1, nBolos, NULL, data, 0, &status);


}

long find_channel_index(fitsfile *fptr, const char * field){
	long nBolos;
	char **data;

	// read the channel list
	read_channels(fptr,data,nBolos);

	long idet;
	//find the fits index of the bolo and return it
	//fits index are 1 indexed, so indexes starts at 1
	for (idet = 1; idet <= nBolos; idet++)
		if (strcmp(field,data[idet-1]) == 0)
			return idet;

	cout << "EE - " << field << " not found" << endl;
	return -1L;

}

void read_ra_from_fits(string filename, string field, double *& ra, long & ns){
	//TODO : Handle angle unit to transform to a common internal known unit

	// HIPE like format

	fitsfile *fptr;
	int status = 0;
	int rowIndex = 0, naxis = 0, anynul;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// Retrieve the row of the specified channel
	rowIndex = find_channel_index(fptr, field.c_str());

	// ---------------------------------------------
	// Move ptr to signal hdu
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "ra", NULL, &status))
		fits_report_error(stderr, status);

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
	ra = new double[ns];

	// ---------------------------------------------
	// Retrieve the corresponding row
	fpixel[0] = 1;
	fpixel[1] = rowIndex;
	if (fits_read_pix(fptr, TDOUBLE, fpixel, ns, 0, ra, &anynul, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// close file
	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

}

//TODO : read_ra_from_fits and read_dec_from_fits could be merged with an additionnal argument
void read_dec_from_fits(string filename, string field, double *& dec, long & ns){
	//TODO : Handle angle unit to transform to a common internal known unit

	// HIPE like format

	fitsfile *fptr;
	int status = 0;
	int rowIndex = 0, naxis = 0, anynul;
	long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// Retrieve the row of the specified channel
	rowIndex = find_channel_index(fptr, field.c_str());

	// ---------------------------------------------
	// Move ptr to signal hdu
	if (fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "dec", NULL, &status))
		fits_report_error(stderr, status);

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
	dec = new double[ns];

	// ---------------------------------------------
	// Retrieve the corresponding row
	fpixel[0] = 1;
	fpixel[1] = rowIndex;
	if (fits_read_pix(fptr, TDOUBLE, fpixel, ns, 0, dec, &anynul, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// close file
	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

}

