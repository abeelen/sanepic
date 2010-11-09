#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parser_functions.h"
#include "tools.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()
#include <sysexits.h>


extern "C" {
#include "nrutil.h"
#include <fitsio.h>
}



using namespace std;

void copy_ref_pos(fitsfile * fptr, fitsfile *outfptr, string name, long min_sample, long max_sample)
/*! Copy RA DEC and PHI (reference detector) tables from input fits to output */
{

	long ns_temp; // temporary value of ns, needed only to read input data
	double *RA,*DEC,*PHI;
	double *RA_bis, *DEC_bis, *PHI_bis;
	int status=0; // fits error status
	long ns_final = max_sample - min_sample +1; // total number of samples to copy

	// Read original tables
	read_ReferencePosition_from_fits(name, RA, DEC, PHI, ns_temp);

	RA_bis = new double [ns_final];
	DEC_bis = new double [ns_final];
	PHI_bis = new double [ns_final];

	// copy corresponding data
	for(long ii = 0; ii< ns_final; ii++){
		RA_bis[ii]=RA[min_sample+ii]*15.0;
		DEC_bis[ii]=DEC[min_sample+ii];
		PHI_bis[ii]=PHI[min_sample+ii];
	}

	// copy header to output
	if(fits_movnam_hdu(fptr, BINARY_TBL, (char*) "reference position", NULL, &status))
		fits_report_error(stderr, status);
	if(fits_copy_header(fptr, outfptr, &status))
		fits_report_error(stderr, status);
	if(fits_movnam_hdu(outfptr, BINARY_TBL, (char*) "reference position", NULL, &status))
		fits_report_error(stderr, status);

	// insert columns
	fits_write_col(outfptr, TDOUBLE, 1, 1, 1, ns_final, RA_bis, &status);
	fits_write_col(outfptr, TDOUBLE, 2, 1, 1, ns_final, DEC_bis, &status);
	fits_write_col(outfptr, TDOUBLE, 3, 1, 1, ns_final, PHI_bis, &status);
	if(fits_update_key(outfptr, TLONG, (char*)"NAXIS2", &ns_final, (char*)"Number of rows", &status))
		fits_report_error(stderr, status);

	delete [] RA;
	delete [] DEC;
	delete [] PHI;
	delete [] RA_bis;
	delete [] DEC_bis;
	delete [] PHI_bis;

}

void copy_time(fitsfile * fptr, fitsfile *outfptr, double *time, long min_sample, long max_sample)
/*! Copy resized time table from input fits to output */
{

	int status=0; // fits error status
	double *time_bis;
	long ns_final = max_sample - min_sample +1; // total number of samples to copy

	time_bis = new double [ns_final];

	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "time", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);  // copy header

	// copy corresponding data
	for(long ii = 0; ii< ns_final; ii++){
		time_bis[ii]=time[min_sample+ii];
	}

	// write data in output file and update header
	fits_write_col(outfptr, TDOUBLE, 1, 1, 1, ns_final, time_bis, &status);
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);

	delete [] time_bis;

}

void copy_signal(fitsfile * fptr, fitsfile *outfptr, string name, long min_sample, long max_sample,  std::vector<std::string> det, long ndet)
/*! Copy resized signal table from input fits to output */
{

	int status=0; // fits error status
	long ns_temp; // temporary value of ns, needed only to read input data
	long ns_final = max_sample - min_sample +1; // total number of samples to copy

	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "signal", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);  // copy header
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);

	double *signal, *signal_bis;
	signal_bis = new double [ns_final];
	for(long jj=0;jj<ndet;jj++){

		read_signal_from_fits(name, det[jj], signal, ns_temp); // read input data

		// copy corresponding data
		for(long ii = 0; ii< ns_final; ii++)
			signal_bis[ii]=signal[min_sample+ii];

		string field= det[jj];
		long rowIndex = find_channel_index(fptr, field.c_str()); // find the row index for the channel named boloname[jj]
		long fpixel[2]={1,rowIndex}; // write a row which row number is roxIndex
		fits_write_pix(outfptr, TDOUBLE, fpixel, ns_final, signal_bis, &status);
	}

	delete [] signal_bis;
	delete [] signal;

}

void copy_mask(fitsfile * fptr, fitsfile *outfptr,  string name, long min_sample, long max_sample,  std::vector<std::string> det, long ndet)
/*! Copy resized mask table from input fits to output */
{

	int status=0; // fits error status
	long ns_temp; // temporary value of ns, needed only to read input data
	long ns_final = max_sample - min_sample +1; // total number of samples to copy

	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "mask", NULL, &status);
	fits_copy_header(fptr, outfptr, &status); // copy header
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status); // update output header


	int *mask, *mask_bis;
	mask_bis = new int [ns_final];
	for(long jj=0;jj<ndet;jj++){

		read_flag_from_fits(name, det[jj], mask, ns_temp); // read detector flag row
		for(long ii = 0; ii< ns_final; ii++)
			mask_bis[ii]=mask[min_sample+ii];

		string field= det[jj];
		long rowIndex = find_channel_index(fptr, field.c_str()); // find the row index for the channel named boloname[jj]
		long fpixel[2]={1,rowIndex}; // write the mask which row number is rowIndex
		fits_write_pix(outfptr, TINT, fpixel, ns_final, mask_bis, &status);
	}

	delete [] mask_bis;
	delete [] mask;
}


void copy_RA_DEC(fitsfile * fptr, fitsfile *outfptr, string name, long min_sample, long max_sample,  std::vector<std::string> det, long ndet)
/*! Copy resized RA and DEC tables from input fits to output */
// Only for HIPE format
{

	int status=0; // fits error status
	long ns_temp; // temporary value of ns, needed only to read input data
	long ns_final = max_sample - min_sample +1; // total number of samples to copy

	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "ra", NULL, &status);
	fits_copy_header(fptr, outfptr, &status);  // copy RA header
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);

	fits_movnam_hdu(fptr, IMAGE_HDU, (char*) "dec", NULL, &status);
	fits_copy_header(fptr, outfptr, &status); // copy DEC header
	fits_update_key(outfptr, TLONG, (char*)"NAXIS1", &ns_final, (char*)"Number of rows", &status);

	double *RA, *RA_bis;
	double *DEC, *DEC_bis;
	RA_bis = new double [ns_final];
	DEC_bis = new double [ns_final];
	for(long jj=0;jj<ndet;jj++){

		read_ra_dec_from_fits(name, det[jj], RA, DEC, ns_temp);  // read data from input
		for(long ii = 0; ii< ns_final; ii++){
			RA_bis[ii]=RA[min_sample+ii]; // copy corresponding data
			DEC_bis[ii]=DEC[min_sample+ii];
		}

		string field= det[jj];
		long rowIndex = find_channel_index(fptr, field.c_str()); // find the correct row number in output table
		long fpixel[2]={1,rowIndex};
		fits_movnam_hdu(outfptr, IMAGE_HDU, (char*) "ra", NULL, &status); // move to RA table
		fits_write_pix(outfptr, TDOUBLE, fpixel, ns_final, RA_bis, &status); // Write RA row
		fits_movnam_hdu(outfptr, IMAGE_HDU, (char*) "dec", NULL, &status); // move to DEC table
		fits_write_pix(outfptr, TDOUBLE, fpixel, ns_final, DEC_bis, &status); // Write DEC row
	}

	delete [] RA_bis;
	delete [] RA;
	delete [] DEC_bis;
	delete [] DEC;
}
