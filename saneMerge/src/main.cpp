
#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parse_saneMerge.h"
#include "tools.h"
#include "struct_definition.h"

#include <iostream>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()
#include <sysexits.h>


extern "C" {
#include "nrutil.h"
}

using namespace std;


//! Usage function, print to stdout a standard command line that explains how to run correctly saneMerge
/*!
 * \param name A char* containing program name
 */
void usage(char *name)
{
	cerr << "USAGE: " << name << " ini_file fitsfile_number1.fits fitsfile_number2.fits [any other fitsfile that you want to merge with file1 and 2...]" << endl;

}

/*!
 *  This is organized as :
 *
 *  - parse the command line
 *  - check for existence of directory/files given in the command line
 *  - Print parser output to screen
 *
 *  - for each file, get:
 *      - fits file format (HIPE or SANEPIC)
 *      - Check whole file have the same format
 *
 *  - Generate output fits filename
 *  - Check files compatibility (same channel list) and that time table are crescent (user has to give inputs in a crescent time order)
 *  - Create a new fits file to merge the fits file in a single new one
 *  - Merge files HDU by HDU, creating a new table in the output fits for each
 *
 */

int main(int argc, char *argv[]) {

	int parsed=1;  /* parser error status */

	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */
//	struct detectors det;  /* A structure that contains everything about the detectors names and number */
	std::vector<string> det;
	long ndet;
	string outname =""; /* output filename */
	string output = ""; /* parser output string */
	struct param_common dir;

	printf("\nBeginning of saneMerge:\n\n");

	if (argc<4)
		parsed=1; /* wrong number of arguments */
	else{
		parsed=parse_saneMerge_ini_file(argv, output, argc, dir, samples_struct); /* parse ini file */

		// print parser warning and/or errors
		cout << endl << output << endl;
	}

	if (parsed>0){ /* check for error during parsing phase */
		switch (parsed){

		case 1: printf("Incorrect number of arguments\n");
		usage(argv[0]);
		break;

		case 2 : printf("Please Correct your command line. Exiting !\n");
		break;

		default :
			break;
		}

		return EX_CONFIG;
	}

	// Read file size once for all
	readFramesFromFits(samples_struct, 0);

	string fname = samples_struct.fitsvect[0]; /* first input fits file name */


	int status=0; /* fits error status */
	fitsfile *fptr; /* input fits file pointer */
	fitsfile *outfptr; /* output fits file pointer */

	// generate output filename using input file names
	int format_fits=0;
	// 0 = Unknown format,
	// 1 = RefPos & offsets format (sanepic),
	// 2 = lon/lat format (HIPE),
	// 3 = both sanepic & HIPE

	for(long ii=0; ii<samples_struct.ntotscan-1;ii++){ // for each scan that has to be merged
		format_fits+=test_format(samples_struct.fitsvect[ii]); // Check each scan fits format
		outname += FitsBasename(samples_struct.fitsvect[ii]) + "_merged_with_"; // generate output filename !
	}
	format_fits+=test_format(samples_struct.fitsvect[samples_struct.ntotscan-1]); // add last file informations
	outname += FitsBasename(samples_struct.fitsvect[samples_struct.ntotscan-1]) + ".fits";
	outname = "!" + dir.output_dir + outname; // complete output path + file name

#ifdef DEBUG
	cout << "Creating : \n" << outname << endl << endl;
#endif

	switch(format_fits/samples_struct.ntotscan){ // check compatibility between each input files (same format)

	case 1: cout << "RefPos/offsets format found" << endl;;
	format_fits=1;
	break;

	case 2: cout << "lon/lat format found" << endl;;
	format_fits=2;
	break;

	case 3: cout << "hybrid format found" << endl;;
	format_fits=3;
	break;

	default : cout << "The files you are trying to merge have not the same format or format is unknown. Exiting...\n";
	return EX_CONFIG;
	}

	long ns_total=0; // compute total number of samples in the ouput file
	for(long ii=0; ii<samples_struct.ntotscan;ii++){
		ns_total+=samples_struct.nsamples[ii];
	}

	/* this function tests that the detectors lists are the same in whole fits files,
	 *  also tests that the files have a crescent time reference */
	file_compatibility_verification(dir.data_dir, samples_struct);

	// read the first file detector list (which is the same in every scans)
	read_bolo_list(samples_struct.fitsvect[0], det, ndet);

	cout << "Merging files..." << endl;

	// open final fits file
	if (fits_create_file(&outfptr, outname.c_str(), &status))
		fits_report_error(stderr, status);


	// 1 signal
	copy_signal(outfptr, dir.data_dir, samples_struct, det, ndet, ns_total); // copy signal tables from each file to output file

	// 4 mask
	copy_mask(outfptr, dir.data_dir, samples_struct, det,  ndet, ns_total); // copy flag tables from each file to output file

	// 4 time
	copy_time(outfptr, dir.data_dir, samples_struct, ns_total); // copy time tables from each file to output file


	if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)) // open first input file
		fits_report_error(stderr, status);

	// 6 channels
	copy_channels(fptr, outfptr); // copy channels list from this file to output file

	switch(format_fits) {
	// 0 = Unknown format,
	// 1 = RefPos & offsets format (sanepic),
	// 2 = lon/lat format (HIPE),
	// 3 = both sanepic & HIPE

	case 1:
		copy_ref_pos(outfptr, dir.data_dir, samples_struct, ns_total); // copy reference position tables from each file to output file
		copy_offsets(fptr, outfptr); // copy offsets table from this file to output file
		break;
	case 3:
		copy_ref_pos(outfptr, dir.data_dir, samples_struct, ns_total); // copy reference position tables from each file to output file
		copy_offsets(fptr, outfptr); // copy offsets table from this file to output file
		// and...
	case 2:
		copy_LON_LAT(outfptr, dir.data_dir, samples_struct, det, ndet, ns_total);
		break;
	}

	if (fits_close_file(fptr, &status)) // close first input file
		fits_report_error(stderr, status);


	if (fits_close_file(outfptr, &status)) // close output file
		fits_report_error(stderr, status);

	cout << "done.\n";

	return EXIT_SUCCESS;

}
