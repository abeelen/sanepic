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


#include "positionsIO.h"


extern "C" {
#include <nrutil.h>
#include <fitsio.h>
}

using namespace std;

/*!
 * Reads a detector list in a .txt file
 * Returns a vector of string containing the name of the considered channels
 */
void read_strings(string fname, std::vector<string> &bolos) {
	string line;

	ifstream inputFile(fname.c_str());
	if (!inputFile.is_open()) {
		cerr << "Error opening bolometer file '" << fname << "'. Exiting.\n";
		exit(1);
	}

	while (!inputFile.eof()) {
		getline(inputFile, line);

		line.erase(0, line.find_first_not_of(" \t")); // remove leading white space
		if (line.empty() || line[0] == '#')
			continue; // skip if empty or commented
		line = line.substr(0, line.find_first_of(" \t")); // pick out first word

		bolos.push_back(line);
	}

	inputFile.close();
}

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

void read_bolo_offsets_from_fits(string filename, string field, double *offsets){

	fitsfile *fptr;
	int status = 0;
	//long naxes[2] = { 1, 1 }, fpixel[2] = { 1, 1 };
	long nBolos, repeat, width;
	int colnum, typecode;
	double temp;

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// read the Channel List
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "Channel Offsets", NULL, &status))
		fits_report_error(stderr, status);

	fits_get_num_rows(fptr, &nBolos, &status);
	fits_get_colnum(fptr, CASEINSEN, (char*) "NAME", &colnum, &status);
	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);

	// Initialize the data container
	char ** data;
	data = new char* [nBolos+1];
	for (long i = 0; i < nBolos; i++) {
		data[i] = new char[repeat];
	}

	fits_read_col(fptr, TSTRING, colnum, 1, 1, nBolos, NULL, data, 0, &status);

	std::vector<string> bolos;
	// convert to string vector and free the container
	bolos.resize(nBolos);
	for (long i = 0; i < nBolos; i++) {
		bolos[i] = data[i];
		free(data[i]);
	}
	free(data);

	//	std::vector<string>::const_iterator it = find (bolos.begin(), bolos.end(),field);
	// +1 because cfitsio is 1 indexed
	long firstrow = distance(bolos.begin(), find (bolos.begin(), bolos.end(),field))+1;

	fits_get_colnum(fptr, CASEINSEN, (char*) "X", &colnum, &status);
	fits_read_col(fptr, TDOUBLE, colnum, firstrow, 1, 1, NULL, &temp, 0, &status);
	offsets[0] = temp/60.0/60.0; // deg

	fits_get_colnum(fptr, CASEINSEN, (char*) "Y", &colnum, &status);
	fits_read_col(fptr, TDOUBLE, colnum, firstrow, 1, 1, NULL, &temp, 0, &status);
	offsets[1] = temp/60.0/60.0; //deg

}


void read_position_from_fits(string filename, double *RA, double *DEC, double *PHI, short *FLAG, long &ns){


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


}
