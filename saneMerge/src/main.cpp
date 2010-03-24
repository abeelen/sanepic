
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
#include "nrcode.h"
}

using namespace std;



void usage(char *name)
{
	cerr << "USAGE: " << name << " Output_path fitsfile_number1.fits fitsfile_number2.fits [any other fitsfile that you want to merge with file1 and 2...]" << endl;

}



int main(int argc, char *argv[]) {

	int parsed=1;

	struct samples samples_struct;
	struct detectors det;
	string outdir;
	string outname ="";

	printf("\nBeginning of saneFix:\n\n");

	if (argc<4)
		parsed=1;
	else
		parsed=parse_saneMerge_ini_file(argv, argc, outdir, samples_struct);

	if (parsed>0){
		switch (parsed){

		case 1: printf("Incorrect number of arguments\n");
		usage(argv[0]);
		break;

		case 2 : printf("Please give an existing output directory. Exiting !\n");
		break;

		default :;
		}

		exit(EXIT_FAILURE);
	}


	//	cout << samples_struct.ntotscan << endl;
	//	cout << samples_struct.fitsvect[0] << endl << samples_struct.fitsvect[1] << endl;
	string fname = samples_struct.fitsvect[0];


	int status=0;
	fitsfile *fptr;
	fitsfile *outfptr;
	std::ostringstream oss;


	int format_fits=0;
	for(long ii=0; ii<samples_struct.ntotscan-1;ii++){
		format_fits+=test_format(samples_struct.fitsvect[ii]);
		oss << samples_struct.fitsvect[ii];
		string filename = oss.str();
		outname += Basename(filename) + "_merged_with_";
		oss.str("");
	}

	format_fits+=test_format(samples_struct.fitsvect[samples_struct.ntotscan-1]);
	oss << samples_struct.fitsvect[samples_struct.ntotscan-1];
	string filename = oss.str();
	outname += Basename(filename) + ".fits";
	oss.str("");
	outname = "!" + outdir + outname;
	cout << outname << endl;

	switch(format_fits/samples_struct.ntotscan){

	case 1: cout << "HIPE format found\n";
	format_fits=1;
	break;

	case 2: cout << "SANEPIC format found\n";
	format_fits=2;
	break;

	default : cout << "The files you are trying to merge have not the same format. Exiting...\n";
	exit(EXIT_FAILURE);
	}

	long ns_total=0;
	for(long ii=0; ii<samples_struct.ntotscan;ii++){
		ns_total+=samples_struct.nsamples[ii];
	}

	file_compatibility_verification(samples_struct);

	read_bolo_list(samples_struct.fitsvect[0], det);

	// open final fits file
	if (fits_create_file(&outfptr, outname.c_str(), &status))
		fits_report_error(stderr, status);

	if(format_fits==1){ // HIPE format
		// 1 signal
		copy_signal(outfptr, samples_struct, det, ns_total);

		// 2 RA 3 DEC
		copy_RA_DEC(outfptr, samples_struct, det, ns_total);

		// 4 mask
		copy_mask(outfptr, samples_struct, det, ns_total);

		// 5 time
		copy_time(outfptr, samples_struct, ns_total);


		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
			fits_report_error(stderr, status);
		// 6 channels
		copy_channels(fptr, outfptr);
		if (fits_close_file(fptr, &status))
			fits_report_error(stderr, status);

		// 7 ref pos
		copy_ref_pos(outfptr, samples_struct, ns_total);

		// 8 offsets
		if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
			fits_report_error(stderr, status);
		copy_offsets(fptr, outfptr);
		if (fits_close_file(fptr, &status))
			fits_report_error(stderr, status);

	}else{ // sanepic format

		// 1 ref pos
		copy_ref_pos(outfptr,samples_struct, ns_total);

		// 2 offsets
		copy_offsets(fptr, outfptr);

		//		if(ii==0)
		// 3 channels
		copy_channels(fptr, outfptr);

		// 4 time
		copy_time(outfptr, samples_struct, ns_total);

		// 5 signal
		copy_signal(outfptr, samples_struct, det, ns_total);

		// 6 mask
		copy_mask(outfptr, samples_struct, det, ns_total);
	}


	if (fits_close_file(outfptr, &status))
		fits_report_error(stderr, status);

	cout << "End of saneMerge\n";

}
