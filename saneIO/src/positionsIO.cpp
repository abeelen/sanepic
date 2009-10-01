/*
 * positionsIO.cpp
 *
 *  Created on: 31 août 2009
 *      Author: abeelen
 */

#include <string>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdlib>

#include <math.h>
#include "positionsIO.h"

extern "C" {
#include <nrutil.h>
#include <fitsio.h>
}

using namespace std;

/*!
 * Reads the detectors offsets in a .txt file
 * Returns an array containing the considered channel offsets + the source offsets
 */
void read_bolo_offsets(string field, string file_BoloOffsets, float *scoffsets, double *offsets){

	double lel, xel;
	long temp1, temp2, temp3;
	int nobolo = 1;

	char boloname[100];
	FILE *fp;


	if ((fp = fopen(file_BoloOffsets.c_str(),"r")) == NULL){
		cerr << "ERROR: Can't find offset file. Exiting. \n";
		exit(1);
	}
	while (fscanf(fp, "%s%ld%ld%ld%lf%lf\n", boloname, &temp1, &temp2, &temp3, &lel, &xel) != EOF) {
		if (field == boloname) {
			nobolo = 0;
			if (temp3 == 250){
				offsets[0] = xel/60.0/60.0 - scoffsets[1];
				offsets[1] = lel/60.0/60.0 + scoffsets[0];
			}
			if (temp3 == 350){
				offsets[0] = xel/60.0/60.0 - scoffsets[3];
				offsets[1] = lel/60.0/60.0 + scoffsets[2];
			}
			if (temp3 == 500){
				offsets[0] = xel/60.0/60.0 - scoffsets[5];
				offsets[1] = lel/60.0/60.0 + scoffsets[4];
			}
		}
	}
	fclose (fp);


	if (nobolo){
		cerr << "Bolometer name not found in offset list" << endl;
		exit(1);
	}


}

void read_bolo_offsets_from_fits(string filename, string field, double * offsets){

	fitsfile *fptr;
	int status = 0;
	//long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
	int colnum;
	double *temp;
	temp = new double[2];

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// read the Channel List
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "Channel Offsets", NULL, &status))
		fits_report_error(stderr, status);

//	fits_get_num_rows(fptr, &nBolos, &status);
//	fits_get_colnum(fptr, CASEINSEN, (char*) "NAME", &colnum, &status);
//	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
//
//	// Initialize the data container
//	char ** data;
//	data = new char* [nBolos+1];
//	for (long i = 0; i < nBolos; i++) {
//		data[i] = new char[repeat];
//	}
//
//	fits_read_col(fptr, TSTRING, colnum, 1, 1, nBolos, NULL, data, 0, &status);
//
//	std::vector<string> bolos;
//	// convert to string vector and free the container
//	bolos.resize(nBolos);
//	for (long i = 0; i < nBolos; i++) {
//		bolos[i] = data[i];
//		free(data[i]);
//	}
//	free(data);
//
//	//	std::vector<string>::const_iterator it = find (bolos.begin(), bolos.end(),field);
//	// +1 because cfitsio is 1 indexed
//	long firstrow = distance(bolos.begin(), find (bolos.begin(), bolos.end(),field))+1;
//
	// Look-up for the Bolometer
	fits_get_colnum(fptr, CASEINSEN, (char *) field.c_str(), &colnum, &status);
	fits_read_col(fptr, TDOUBLE, colnum, 1, 1, 2, NULL, temp, 0, &status);

	// transform arcsec o deg
	for (int ii=0; ii<2; ii++)	offsets[ii] = temp[ii]/3600;

	//	//	std::vector<string>::const_iterator it = find (bolos.begin(), bolos.end(),field);
//	// +1 because cfitsio is 1 indexed
//	long firstrow = distance(bolos.begin(), find (bolos.begin(), bolos.end(),field))+1;
//
//	fits_get_colnum(fptr, CASEINSEN, (char*) "X", &colnum, &status);
//	fits_read_col(fptr, TDOUBLE, colnum, firstrow, 1, 1, NULL, &temp, 0, &status);
//	//offsets[0] = temp/60.0/60.0; // deg
//	offsets[0] = temp/3600;
//
//	fits_get_colnum(fptr, CASEINSEN, (char*) "Y", &colnum, &status);
//	fits_read_col(fptr, TDOUBLE, colnum, firstrow, 1, 1, NULL, &temp, 0, &status);
//	//offsets[1] = temp/60.0/60.0; //deg
//	offsets[1] = temp/3600;

	delete [] temp;
}


void read_all_bolo_offsets_from_fits(string filename, std::vector<string> bolonames, double **& offsets){

	fitsfile *fptr;
	int status = 0;
	int colnum;
	double * temp;
	temp = new double[2];

	unsigned long ndet = bolonames.size();

	offsets = dmatrix((long)0,ndet-1,(long)0,2-1);

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// read the Channel List
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "Channel Offsets", NULL, &status))
		fits_report_error(stderr, status);
	// match read offsets with requested offsets
	for (unsigned long idet=0;idet<ndet;idet++){


		fits_get_colnum(fptr, CASEINSEN, (char *) bolonames[idet].c_str(), &colnum, &status);
		fits_read_col(fptr, TDOUBLE, colnum, 1, 1, 2, NULL, temp, 0, &status);

		// transform arcsec to deg
		for (int ii=0; ii<2; ii++)	offsets[idet][ii] = temp[ii]/3600;

	}

	delete [] temp;
}

// TODO : not used anymore
//void read_data_from_fits(string filename, double *data, double *data2, double *data3, short *data4, bool flag, short *data5, long &ns, string field){
void read_position_from_fits(string filename, double *RA, double *DEC, double *PHI, short *FLAG, bool flag, short *mask, long &ns, string field){

	/*int sizetype;
	char test[2];


	test[0] = type;
	test[1] = '\0';
	string typestr = string(test);
	//  printf("type = %s\n",test);
	free(test);

	if (typestr == "d") sizetype = 8;
	if (typestr == "c") sizetype = 1;
	 */
	////////////////////////////////////

	fitsfile *fptr;
	int status = 0;
	//long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
	long repeat, width;
	int colnum, typecode;
	//double temp;

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// read the Channel List
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "Reference Position", NULL, &status))
		fits_report_error(stderr, status);

	fits_get_num_rows(fptr, &ns, &status);
	fits_get_colnum(fptr, CASEINSEN, (char*) "RA", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
	fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ns, NULL, RA, 0, &status);

	// divide by 15 to get the same value as before
	for(int ii = 0; ii<ns; ii++)
		RA[ii]=RA[ii]/15; // TODO : Pourquoi ??????? 15 ??


	//fits_get_num_rows(fptr, &ns, &status);
	fits_get_colnum(fptr, CASEINSEN, (char*) "DEC", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
	fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ns, NULL, DEC, 0, &status);


	//fits_get_num_rows(fptr, &ns, &status);
	fits_get_colnum(fptr, CASEINSEN, (char*) "PHI", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
	fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ns, NULL, PHI, 0, &status);


	//fits_get_num_rows(fptr, &ns, &status);
	fits_get_colnum(fptr, CASEINSEN, (char*) "FLAG", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
	fits_read_col(fptr, TSHORT, colnum, 1, 1, ns, NULL, FLAG, 0, &status);

	if(flag){
		// read the Channel List
		if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "mask", NULL, &status))
			fits_report_error(stderr, status);

		fits_get_num_rows(fptr, &ns, &status); // optionnel ???
		fits_get_colnum(fptr, CASEINSEN, (char*)field.c_str(), &colnum, &status);
		fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
		fits_read_col(fptr, TSHORT, colnum, 1, 1, ns, NULL, mask, 0, &status);

	}



}

void read_ReferencePosition_from_fits(string filename, double *&RA, double *&DEC, double *&PHI, short *&FLAG, long &ns){

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

	// Read FLAG
	fits_get_colnum(fptr, CASEINSEN, (char*) "FLAG", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
	fits_read_col(fptr, TSHORT, colnum, 1, 1, ns, NULL, FLAG, 0, &status);

	// check for bad values
	for (long ii = 0; ii<ns; ii++)
		if (isnan(RA[ii]) || isnan(DEC[ii]))
			FLAG[ii] = 1;

}


void read_flpoint_from_fits(string filename, short *FLAG){


        fitsfile *fptr;
        int status = 0;
        //long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
        long repeat, width;
        int colnum, typecode;
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
        fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
        fits_read_col(fptr, TSHORT, colnum, 1, 1, ns, NULL, FLAG, 0, &status);

        // close file
        if(fits_close_file(fptr, &status)) // ajout mat 15/09
                fits_report_error(stderr, status);

}

void read_flag_from_fits(string filename, string field, short *&mask, long &ns){


        fitsfile *fptr;
        int status = 0;
        long repeat, width;
        int colnum, typecode;

        if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
                fits_report_error(stderr, status);

        // ---------------------------------------------
        if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "mask", NULL, &status))
                fits_report_error(stderr, status);

        // Get array size ...
        fits_get_num_rows(fptr, &ns, &status);

        // ... and allocate memory
        mask = new short[ns];

        // Read the actual data
        fits_get_colnum(fptr, CASEINSEN, (char*)field.c_str(), &colnum, &status);
        fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
        fits_read_col(fptr, TSHORT, colnum, 1, 1, ns, NULL, mask, 0, &status);

        // close file
        if(fits_close_file(fptr, &status))
                fits_report_error(stderr, status);

}

void read_signal_from_fits(string filename, double *signal, string field){ // tested in sanepos, work fine



        fitsfile *fptr;
        int status = 0;
        //long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
        long repeat, width;
        int colnum, typecode;
        long ns;
        //double temp;

        if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
                fits_report_error(stderr, status);


        // ---------------------------------------------
        // Move ptr to signal hdu
        if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "signal", NULL, &status))
                fits_report_error(stderr, status);

        fits_get_num_rows(fptr, &ns, &status);
        fits_get_colnum(fptr, CASEINSEN, (char*) field.c_str(), &colnum, &status);
        fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
        fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ns, NULL, signal, 0, &status);

        // close file
        if(fits_close_file(fptr, &status)) // ajout mat 15/09
                fits_report_error(stderr, status);
        //print_fits_error(status);



}


