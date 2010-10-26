
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


extern "C" {
#include "nrutil.h"
}

using namespace std;



void usage(char *name)
{
	cerr << "USAGE: " << name << " Output_path fitsfile_number1.fits fitsfile_number2.fits [any other fitsfile that you want to merge with file1 and 2...]" << endl;

}



int main(int argc, char *argv[]) {

	int parsed=1;  /* parser error status */

	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */
	struct detectors det;  /*! A structure that contains everything about the detectors names and number */
	string outdir; /*! output directory */
	string outname =""; /*! output filename */
	string output = "";

	printf("\nBeginning of saneMerge:\n\n");

	if (argc<4)
		parsed=1; /* wrong number of arguments */
	else{
		parsed=parse_saneMerge_ini_file(argv, output, argc, outdir, samples_struct); /* parse ini file */

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

		default :;
		}

		exit(EXIT_FAILURE);
	}

	string fname = samples_struct.fitsvect[0]; /* input fits file name */


	int status=0; /* fits error status */
	fitsfile *fptr; /*! input fits file pointer */
	fitsfile *outfptr; /*! output fits file pointer */


	int format_fits=0; // 1= HIPE, 2 = SANEPIC
	for(long ii=0; ii<samples_struct.ntotscan-1;ii++){ // for each scan that has to be merged
		format_fits+=test_format(samples_struct.fitsvect[ii]); // Check each scan fits format
		outname += FitsBasename(samples_struct.fitsvect[ii]) + "_merged_with_"; // generate output filename !
	}

	format_fits+=test_format(samples_struct.fitsvect[samples_struct.ntotscan-1]); // add last file informations
	outname += FitsBasename(samples_struct.fitsvect[samples_struct.ntotscan-1]) + ".fits";
	outname = "!" + outdir + outname; // complete output path + file name
	cout << outname << endl;

	switch(format_fits/samples_struct.ntotscan){ // check compatibility between each input files (same format)

	case 1: cout << "HIPE format found\n";
	format_fits=1;
	break;

	case 2: cout << "SANEPIC format found\n";
	format_fits=2;
	break;

	default : cout << "The files you are trying to merge have not the same format or format is unknown. Exiting...\n";
	exit(EXIT_FAILURE);
	}

	long ns_total=0; // compute total number of samples in the ouput file
	for(long ii=0; ii<samples_struct.ntotscan;ii++){
		ns_total+=samples_struct.nsamples[ii];
	}

	/* this function tests that the detectors lists are the same in whole fits files,
	 *  also tests that the files have a crescent time reference */
	file_compatibility_verification(samples_struct);

	// read the first file detector list (which is the same in every scans)
	read_bolo_list(samples_struct.fitsvect[0], det);

	// open final fits file
	if (fits_create_file(&outfptr, outname.c_str(), &status))
		fits_report_error(stderr, status);

	if(format_fits==1){ // HIPE format

		// 1 signal
		copy_signal(outfptr, samples_struct, det, ns_total); // copy signal tables from each file to output file

		// 2 RA 3 DEC
		copy_RA_DEC(outfptr, samples_struct, det, ns_total); // copy RA and DEC tables from each file to output file

		// 4 mask
		copy_mask(outfptr, samples_struct, det, ns_total); // copy flag tables from each file to output file

		// 5 time
		copy_time(outfptr, samples_struct, ns_total); // copy time tables from each file to output file


		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)) // open first input file
			fits_report_error(stderr, status);

		// 6 channels
		copy_channels(fptr, outfptr); // copy channels list from this file to output file

		if (fits_close_file(fptr, &status)) // close first input file
			fits_report_error(stderr, status);

		// 7 ref pos
		copy_ref_pos(outfptr, samples_struct, ns_total); // copy reference position tables from each file to output file

		// 8 offsets
		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)) // open first input file
			fits_report_error(stderr, status);

		copy_offsets(fptr, outfptr); // copy offsets table from this file to output file

		if (fits_close_file(fptr, &status)) // close first input file
			fits_report_error(stderr, status);

	}else{ // sanepic format

		// 1 ref pos
		copy_ref_pos(outfptr, samples_struct, ns_total); // copy reference position tables from each file to output file


		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status)) // open first input file
			fits_report_error(stderr, status);

		// 2 offsets
		copy_offsets(fptr, outfptr); // copy offsets table from this file to output file

		// 3 channels
		copy_channels(fptr, outfptr); // copy channels list from this file to output file

		if (fits_close_file(fptr, &status)) // close first input file
			fits_report_error(stderr, status);

		// 4 time
		copy_time(outfptr, samples_struct, ns_total); // copy time tables from each file to output file

		// 5 signal
		copy_signal(outfptr, samples_struct, det, ns_total); // copy signal tables from each file to output file

		// 6 mask
		copy_mask(outfptr, samples_struct, det, ns_total); // copy flag tables from each file to output file
	}


	if (fits_close_file(outfptr, &status)) // close output file
		fits_report_error(stderr, status);

	//clean up
	delete [] samples_struct.nsamples;

	cout << "End of saneMerge\n";

}
