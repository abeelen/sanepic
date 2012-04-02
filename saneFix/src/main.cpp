#include <iostream>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()
#include <sysexits.h>

#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parser_functions.h"
#include "tools.h"
#include "struct_definition.h"
#include "error_code.h"

extern "C" {
#include "nrutil.h"
}

using namespace std;

#ifdef PARA_FRAME
#include "mpi.h"
#endif

/*!
 *  This is organized as :
 *
 *  - parse the input ini file
 *  - check for existence of directory/files pointed from the ini file
 *  - Print parser output to screen
 *
 *  - for each file :
 *  	- Assign a rank to treat the fits file
 *      - Read saneCheck information file
 *      - Get fits file format (HIPE or SANEPIC)
 *      - Suppress time gaps indices that points to beginning or ending flagged data (that will not be copied in the ouput fixed file)
 *      - Create a new fits file
 *      - Fix the input fits file HDU by HDU, by adding flagged data to fill time gaps
 *
 */

int main(int argc, char *argv[]) {

	int rank, size; /* MPI processor rank and MPI total number of used processors */

#ifdef PARA_FRAME

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#else
	size = 1;
	rank = 0;
#endif

	if(rank==0)
		cout << "\nBeginning of saneFix:\n\n";

	struct samples samples_struct; /* fits file list + number of scans */
	struct param_common dir; /* directories : temporary input and output */
	struct param_sanePos pos_param;
	struct param_saneProc proc_param;
	struct param_sanePic sanePic_struct;
	struct param_saneInv saneInv_struct;
	struct param_sanePS structPS;
	struct param_saneCheck check_struct;

	long iframe_min=0, iframe_max=0; /* frame number min and max each processor has to deal with */

	std::vector<long> indice; /* gap indexes (sample index) */
	std::vector <long> add_sample; /* number of samples to add per gap */
	double fsamp; /* sampling frequency */
	string output ="";

	uint16_t mask_sanefix = INI_NOT_FOUND | DATA_INPUT_PATHS_PROBLEM | OUPUT_PATH_PROBLEM | TMP_PATH_PROBLEM |
			FSAMP_PROBLEM | FITS_FILELIST_NOT_FOUND; // 0x410f

	//	if (rank==0){ // root parse ini file and fill the structures. Also print warnings or errors

	uint16_t parsed=0x0000; // parser error status
	uint16_t compare_to_mask; // parser error status


	if (argc<2)/* not enough argument */
		compare_to_mask=0x001;
	else {
		/* parse ini file and fill structures */
		parsed=parser_function(argv[1], output, dir, samples_struct, pos_param, proc_param,
				structPS, saneInv_struct, sanePic_struct, size, rank);

		compare_to_mask = parsed & mask_sanefix;

		// print parser warning and/or errors
		if (rank == 0)
			cout << endl << output << endl;
	}


	// in case there is a parsing error or the dirfile format file was not created correctly
	if(compare_to_mask>0x0000){

		switch (compare_to_mask){/* error during parsing phase */

		case 0x0001: cout << "Please run " << argv[0] << " using a correct *.ini file\n";
		break;

		default : cout << "Wrong program options or argument. Exiting !\n";
		break;


		}

#ifdef PARA_FRAME
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return EX_CONFIG;
	}
	//	}


#ifdef PARA_FRAME

	MPI_Barrier(MPI_COMM_WORLD);

	if(configure_PARA_FRAME_samples_struct(dir.output_dir, samples_struct, rank, size, iframe_min, iframe_max)){
		MPI_Abort(MPI_COMM_WORLD, 1);
		return EX_IOERR;
	}

	MPI_Barrier(MPI_COMM_WORLD);


#else
	iframe_min = 0;
	iframe_max = samples_struct.ntotscan;
#endif


	if(rank==0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, pos_param,  proc_param,
				structPS, sanePic_struct, saneInv_struct);

		cout << "\nFixing Files..." << endl << endl;
	}

	for (long iframe=iframe_min;iframe<iframe_max;iframe++){

			int format_fits; // 1 = HIPE, 2 = Sanepic
			std::vector <long> suppress_time_sample;
			long init_num_delete =0;
			long end_num_delete =0;

#ifdef DEBUG
			cout  << "[ " << rank << " ] " << "Fixing file : " << samples_struct.fitsvect[iframe] << endl;
#endif
			// fits files pointer
			fitsfile * fptr;
			fitsfile *outfptr;

			// lookfor saneCheck log files
			if((read_indices_file(samples_struct.fitsvect[iframe],dir,  indice, fsamp, init_num_delete, end_num_delete))){
				cout << "[ " << rank << " ] " << "Skipping file : " << samples_struct.fitsvect[iframe] << ". Please run saneCheck on this file before\n";
				continue; // log file were not found for the 'iframe'th fits file
			}


			format_fits=test_format(samples_struct.fitsvect[iframe]); // get fits file format
			if(format_fits==0){
				cerr << "[ " << rank << " ] " << "input fits file format is undefined : " << samples_struct.fitsvect[iframe] << " . Skipping file...\n";
				continue;
			}

			refresh_indice(fsamp, init_num_delete, end_num_delete, indice, samples_struct.nsamples[iframe]);

			// compute the number of sample that must be added to fill the gaps and have a continous timeline
			long samples_to_add=how_many(samples_struct.fitsvect[iframe], samples_struct.nsamples[iframe] ,indice, fsamp, add_sample,suppress_time_sample);
#ifdef DEBUG
			cout << "[ " << rank << " ] " << "samptoadd : " << samples_to_add << endl;
			cout << "[ " << rank << " ] " << "samples to delete (begin/end) : (" << init_num_delete << "/" <<  end_num_delete << ")" << endl;
			cout << "[ " << rank << " ] " << "total : " << samples_struct.nsamples[iframe] + samples_to_add - init_num_delete - end_num_delete << endl;
#endif
			// total number of samples in the fixed fits file
			long ns_total = samples_struct.nsamples[iframe] + samples_to_add - init_num_delete - end_num_delete;

			if(samples_to_add+init_num_delete+end_num_delete==0){
				cout << "[ " << rank << " ]" << samples_to_add << " : " << init_num_delete << " : " << end_num_delete << endl;
				cout << "[ " << rank << " ] " << "Nothing to do for : " << samples_struct.fitsvect[iframe] << " . Skipping file...\n";
				continue;
			}

			string fname2 = "!" + dir.output_dir + FitsBasename(samples_struct.fitsvect[iframe]) + ".fits"; // output fits filename
			string fname=samples_struct.fitsvect[iframe]; // input fits filename

			int status=0; // fits error status
			std::vector<string> det;
			long ndet;

			// open original fits file
			if (fits_open_file(&fptr, fname.c_str(), READONLY, &status))
				fits_report_error(stderr, status);

			// create ouput fixed fits file
			if (fits_create_file(&outfptr, fname2.c_str(), &status))
				fits_report_error(stderr, status);

			// read channels list
			read_bolo_list(samples_struct.fitsvect[iframe], det, ndet);

			// Copy primary Header
			fits_copy_header(fptr, outfptr, &status);

			if(format_fits==1){ // HIPE format
#ifdef DEBUG
				cout << "[ " << rank << " ] " << "HIPE format found\n";
#endif

				// 1 signal
				fix_signal(fptr, outfptr, samples_struct.fitsvect[iframe], ns_total, det, ndet, indice, add_sample, suppress_time_sample, init_num_delete);

				// 2 LON 3 LAT
				fix_LON_LAT(fptr, outfptr, samples_struct.fitsvect[iframe], ns_total, det, ndet, indice, add_sample, suppress_time_sample, init_num_delete);

				// 4 mask
				fix_mask(fptr, outfptr, samples_struct.fitsvect[iframe], ns_total, det, ndet, indice, add_sample, suppress_time_sample, init_num_delete);

				// 5 time
				fix_time_table(fptr, outfptr, samples_struct.fitsvect[iframe], ns_total, indice, add_sample, samples_struct.nsamples[iframe], fsamp,suppress_time_sample, init_num_delete);

				// 6 channels
				copy_channels(fptr, outfptr);

				// 7 reference positions
				fix_ref_pos(fptr, outfptr, samples_struct.fitsvect[iframe], ns_total, indice, add_sample, suppress_time_sample, init_num_delete);

				// 8 offsets
				copy_offsets(fptr, outfptr);

			}else{ // format =2
#ifdef DEBUG
				cout << "[ " << rank << " ] " << "SANEPIC format found\n";
#endif

				// 1 ref pos
				fix_ref_pos(fptr, outfptr, samples_struct.fitsvect[iframe], ns_total, indice, add_sample, suppress_time_sample, init_num_delete);

				// 2 offsets
				copy_offsets(fptr, outfptr);

				// 3 channels
				copy_channels(fptr, outfptr);

				// 4 time
				fix_time_table(fptr, outfptr, samples_struct.fitsvect[iframe], ns_total, indice, add_sample, samples_struct.nsamples[iframe], fsamp, suppress_time_sample, init_num_delete);

				// 5 signal
				fix_signal(fptr, outfptr, samples_struct.fitsvect[iframe], ns_total, det, ndet, indice, add_sample, suppress_time_sample, init_num_delete);

				// 6 mask
				fix_mask(fptr, outfptr, samples_struct.fitsvect[iframe], ns_total, det, ndet, indice, add_sample, suppress_time_sample, init_num_delete);

			}

			// close both fits files
			if (fits_close_file(fptr, &status))
				fits_report_error(stderr, status);

			if (fits_close_file(outfptr, &status))
				fits_report_error(stderr, status);

			indice.clear();
			add_sample.clear();

	}

	if(rank==0)
		cout << "done." << endl << endl;

#ifdef PARA_FRAME
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif

	if(rank==0)
		cout << "END OF SANEFIX\n";

	return EXIT_SUCCESS;


}
