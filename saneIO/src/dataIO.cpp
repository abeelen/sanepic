
#include <string>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdlib>

#include <math.h>

extern "C" {
#include <nrutil.h>
#include <fitsio.h>
}

#include "dataIO.h"

// #include "getdata.h"

using namespace std;


//int read_data_std(string fname, int frame, int fs, int ns,
//		void* data, string field, char type)
//{
//
//	int sizetype;
//	char test[2];
//	size_t result;
//
//	test[0] = type;
//	test[1] = '\0';
//	string typestr = string(test);
//	//  printf("type = %s\n",test);
//
//
//	FILE *fp;
//
//	if (typestr == "d") sizetype = 8;
//	if (typestr == "c") sizetype = 1;
//
//	string filename = fname + field;
//
//	if((fp = fopen(filename.c_str(),"r"))!=NULL){
//		fseek(fp,(20*frame+fs)*sizetype,SEEK_SET);
//		result = fread(data,sizetype,ns,fp);
//		fclose(fp);
//	}else{
//		printf("Error. Could not open %s. Exiting...\n",filename.c_str());
//		exit(0);
//	}
//
//		return 1;
//}


// int read_data(string fname, int frame, int fs, int ns,
//               void* data, string field, char type)
// {
//   int error_code, nread;

//   char ffname[100];
//   strcpy(ffname,fname.c_str());

//   nread = GetData(ffname, field.c_str(),
//                     frame,fs, /* 1st sframe, 1st samp */
//                     0, ns, /* num sframes, num samps */
//                     type, data,
//                     &error_code);

//   if (error_code != GD_E_OK) {
//     cerr << "    GetData Error while reading "<< field
// 	 << " from " << fname <<":\n";
//     cerr << GD_ERROR_CODES[error_code] << "\n";
//     cerr << " Frame: " << frame << "\n";

//     //exit(0);
//   }

//   if(nread == 0) {
//     cerr << "Warning: nread = 0\n";
//     //exit(0);
//   }

//   return nread;
// }



void read_all_bolo_offsets_from_fits(string filename, std::vector<string> bolonames, double **& offsets){

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
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "Channel Offsets", NULL, &status))
		fits_report_error(stderr, status);
	// match read offsets with requested offsets
	for (long idet=0;idet<ndet;idet++){


		fits_get_colnum(fptr, CASEINSEN, (char *) bolonames[idet].c_str(), &colnum, &status);
		fits_read_col(fptr, TDOUBLE, colnum, 1, 1, 2, NULL, temp, 0, &status);

		// transform arcsec to deg
		for (int ii=0; ii<2; ii++)	offsets[idet][ii] = temp[ii]/3600;

	}

	delete [] temp;

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

//TODO : NOT USED ? => answer : used in saneLib in write_ftrprocess...
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

void read_flag_from_fits(string filename, string field, short *&mask, long &ns){


	fitsfile *fptr;
	int status = 0;
	int colnum;

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "mask", NULL, &status))
		fits_report_error(stderr, status);

	// Get array size ...
	if(fits_get_num_rows(fptr, &ns, &status))
		fits_report_error(stderr, status);

	// ... and allocate memory
	mask = new short[ns];

	// Read the actual data
	fits_get_colnum(fptr, CASEINSEN, (char*)field.c_str(), &colnum, &status);
//	fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
	fits_read_col(fptr, TSHORT, colnum, 1, 1, ns, NULL, mask, 0, &status);

	// close file
	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

}

// TODO: Report size and make allocation here and check in code (see read_ra|dec)
void read_signal_from_fits(string filename, double *signal, string field){

	fitsfile *fptr;
	int status = 0;
	int colnum;
	long ns;

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// Move ptr to signal hdu
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "signal", NULL, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// Retrieve the size of the signal
	if (fits_get_num_rows(fptr, &ns, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// Retrieve the column number corresponding to given channel
	if (fits_get_colnum(fptr, CASEINSEN, (char*) field.c_str(), &colnum, &status))
		fits_report_error(stderr, status);
	// Retrieve all signal
	if (fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ns, NULL, signal, 0, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// close file
	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

}

void read_ra_from_fits(string filename, string field, double *& ra, long & ns){

	// HIPE like format

	fitsfile *fptr;
	int status = 0;
	int colnum;

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// Move ptr to signal hdu
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "ra", NULL, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// Retrieve the size of the signal
	if (fits_get_num_rows(fptr, &ns, &status))
		fits_report_error(stderr, status);

	ra = new double[ns];
	// ---------------------------------------------
	// Retrieve the column number corresponding to given channel
	if (fits_get_colnum(fptr, CASEINSEN, (char*) field.c_str(), &colnum, &status))
		fits_report_error(stderr, status);
	// Retrieve all signal
	if (fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ns, NULL, ra, 0, &status))
		fits_report_error(stderr, status);
	// ---------------------------------------------
	// close file
	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

}

void read_dec_from_fits(string filename, string field, double *& dec, long & ns){

	// HIPE like format

	fitsfile *fptr;
	int status = 0;
	int colnum;

	if (fits_open_file(&fptr, filename.c_str(), READONLY, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// Move ptr to signal hdu
	if (fits_movnam_hdu(fptr, BINARY_TBL, (char*) "dec", NULL, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// Retrieve the size of the signal
	if (fits_get_num_rows(fptr, &ns, &status))
		fits_report_error(stderr, status);

	dec = new double[ns];

	// ---------------------------------------------
	// Retrieve the column number corresponding to given channel
	if (fits_get_colnum(fptr, CASEINSEN, (char*) field.c_str(), &colnum, &status))
		fits_report_error(stderr, status);
	// Retrieve all signal
	if (fits_read_col(fptr, TDOUBLE, colnum, 1, 1, ns, NULL, dec, 0, &status))
		fits_report_error(stderr, status);

	// ---------------------------------------------
	// close file
	if(fits_close_file(fptr, &status))
		fits_report_error(stderr, status);

}
