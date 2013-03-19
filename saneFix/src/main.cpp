#include <iostream>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib> // for exit()
#include <cstdio>  // for printf()
#include <sysexits.h>

#include "InputFileIO.h"
#include "MPIConfiguration.h"
#include "DataIO.h"
#include "SaneFixParse.h"
#include "ParserFunctions.h"
#include "SaneFixTools.h"
#include "StructDefinition.h"
#include "ErrorCode.h"

extern "C" {
#include "nrutil.h"
}

using namespace std;

#ifdef USE_MPI
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


	int      rank,      size; /* MPI processor rank and MPI total number of used processors */
	int  bolo_rank,  bolo_size; /* As for parallel scheme */
	int node_rank, node_size; /* On a node basis, same as *sub* but for frame scheme */

#ifdef USE_MPI

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	MPI_Comm MPI_COMM_NODE, MPI_COMM_MASTER_NODE;
#else
	size = 1;
	rank = 0;
	bolo_size  = 1;
	bolo_rank  = 0;
	node_size = 1;
	node_rank = 0;
#endif


	if(rank==0)
		cout << endl << "Beginning of saneFix: fixing fits files" << endl << endl;

	struct samples samples_struct; /* fits file list + number of scans */
	struct param_common dir; /* directories : temporary input and output */

	struct param_sanePos pos_param;
	struct param_saneProc proc_param;
	struct param_sanePic Pic_param;
	struct param_saneInv Inv_param;
	struct param_sanePS PS_param;
	struct param_saneFix Fix_param;

	string parser_output = "";

	std::vector<long> indice; /* gap indexes (sample index) */
	std::vector <long> add_sample; /* number of samples to add per gap */
	double fsamp; /* sampling frequency */

	uint32_t mask_saneFix = INI_NOT_FOUND | DATA_INPUT_PATHS_PROBLEM | OUPUT_PATH_PROBLEM | TMP_PATH_PROBLEM |
			FSAMP_PROBLEM | FITS_FILELIST_NOT_FOUND; // 0x410f

	uint32_t parsed=0x0000; // parser error status
	uint32_t compare_to_mask; // parser error status

	if (argc<2) {/* not enough argument */
		if (rank == 0)
			cerr << "EE - Please run  " << StringOf(argv[0]) << " with a .ini file" << endl;
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD );
		MPI_Finalize();
#endif
		exit(EXIT_FAILURE);
	} else {
		/* parse ini file and fill structures */
		parsed=parse_saneFix_ini_file(argv[1], parser_output, dir, samples_struct, pos_param, proc_param,
				PS_param, Inv_param, Pic_param, Fix_param, size, rank);

		compare_to_mask = parsed & mask_saneFix;

		// print parser warning and/or errors
		if (rank == 0)
			cout << endl << parser_output << endl;

		// in case there is a parsing error or the dirfile format file was not created correctly
		if(compare_to_mask>0x0000){


			switch (compare_to_mask){/* error during parsing phase */

			case 0x0001:
				if (rank==0)
					cerr << " EE - Please run " << StringOf(argv[0]) << " using a correct *.ini file" << endl;
				break;

			default :
				if (rank==0)
					cerr << "EE - Wrong program options or argument. Exiting ! " <<  "("<< hex << compare_to_mask << ")" << endl;
				break;

			}

#ifdef USE_MPI
			MPI_Finalize();
#endif
			return EX_CONFIG;
		}
	}


	/* ------------------------------------------------------------------------------------*/
	// Start...

	if(rank==0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, pos_param,  proc_param,
				PS_param, Pic_param, Inv_param);
		print_param_saneFix(Fix_param);

	}


#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);

	if(configureMPI(dir.output_dir, samples_struct, rank, size,
			bolo_rank,  bolo_size, node_rank, node_size,
			MPI_COMM_NODE, MPI_COMM_MASTER_NODE)){
		if (rank==0)
			cerr << endl << endl << "Exiting..." << endl;

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return EX_CONFIG;
	}
	MPI_Barrier(MPI_COMM_WORLD);

#endif

	/* ------------------------------------------------------------------------------------*/

	// Read file size once for all
	if (rank == 0)
		readFramesFromFits(samples_struct);
#ifdef USE_MPI
	MPI_Bcast_vector_long(samples_struct.nsamples, 0, MPI_COMM_WORLD);
#endif

	long nFrames = samples_struct.iframe_max-samples_struct.iframe_min;

	for (long iframe = samples_struct.iframe_min +  floor(bolo_rank*nFrames*1.0/bolo_size); iframe < samples_struct.iframe_min + floor((bolo_rank+1)*nFrames*1.0/bolo_size); iframe++){

		int format_fits;
		std::vector <long> suppress_time_sample;
		long init_num_delete = 0;
		long end_num_delete  = 0;

#ifdef DEBUG
		cout  << "[ " << rank << " ] " << "Fixing file : " << samples_struct.fitsvect[iframe] << endl;
#endif
		// fits files pointer
		fitsfile * fptr;
		fitsfile *outfptr;

		// lookfor saneCheck log files
		long ns_dummy;
		int * speedFlag;

		if ( importCheckFits(dir.tmp_dir, samples_struct.fitsvect[iframe], init_num_delete, end_num_delete, fsamp, ns_dummy, speedFlag, indice) ){
			cout <<  "WW - Skipping file : " << samples_struct.fitsvect[iframe] << ". Please run saneCheck on this file before\n";
			continue;
		}


		format_fits=testExtensions(samples_struct.fitsvect[iframe]); // get fits file format
		if( ((format_fits & HIPE_FORMAT) != HIPE_FORMAT) && ((format_fits & SANEPIC_FORMAT) != SANEPIC_FORMAT) ){
			cerr << "[ " << rank << " ] " << "input fits file format is undefined : " << samples_struct.fitsvect[iframe] << " . Skipping file...\n";
			continue;
		}

		refresh_indice(fsamp, init_num_delete, end_num_delete, indice, samples_struct.nsamples[iframe]);

		// compute the number of sample that must be added to fill the gaps and have a continuous timeline
		long samples_to_add=how_many(samples_struct.fitsvect[iframe], samples_struct.nsamples[iframe] ,indice, fsamp, add_sample,suppress_time_sample);
#ifdef DEBUG
		cout << "[ " << rank << " ] " << "samptoadd : " << samples_to_add << endl;
		cout << "[ " << rank << " ] " << "samples to delete (begin/end) : (" << init_num_delete << "/" <<  end_num_delete << ")" << endl;
		cout << "[ " << rank << " ] " << "total : " << samples_struct.nsamples[iframe] + samples_to_add - init_num_delete - end_num_delete << endl;
#endif
		// total number of samples in the fixed fits file
		long ns_total = samples_struct.nsamples[iframe] + samples_to_add - init_num_delete - end_num_delete;

		if(samples_to_add== 0 && init_num_delete == 0 && end_num_delete==0){
			//				cout << "WW - [ " << rank << " ] " << samples_to_add << " : " << init_num_delete << " : " << end_num_delete << endl;
			cout << "WW - [ " << rank << " ] " << "Nothing to do for : " << samples_struct.fitsvect[iframe] << " . Skipping file...\n";
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


		fix_signal(fptr, outfptr, samples_struct.fitsvect[iframe], ns_total, det, ndet, indice, add_sample, suppress_time_sample, init_num_delete);
		fix_mask(  fptr, outfptr, samples_struct.fitsvect[iframe], ns_total, det, ndet, indice, add_sample, suppress_time_sample, init_num_delete, speedFlag);
		fix_time_table(fptr, outfptr, samples_struct.fitsvect[iframe], ns_total, indice, add_sample, samples_struct.nsamples[iframe], fsamp,suppress_time_sample, init_num_delete);
		copy_channels(fptr, outfptr);

		switch(format_fits & EXT_POS) {
		case SANEPIC_FORMAT:
			fix_ref_pos(fptr, outfptr, samples_struct.fitsvect[iframe], ns_total, indice, add_sample, suppress_time_sample, init_num_delete);
			copy_offsets(fptr, outfptr);
			break;
		case BOTH_FORMAT:
			fix_ref_pos(fptr, outfptr, samples_struct.fitsvect[iframe], ns_total, indice, add_sample, suppress_time_sample, init_num_delete);
			copy_offsets(fptr, outfptr);
			// and ...
			/* no break */
		case HIPE_FORMAT:
			fix_LON_LAT(fptr, outfptr, samples_struct.fitsvect[iframe], ns_total, det, ndet, indice, add_sample, suppress_time_sample, init_num_delete);
			break;
		}

		// close both fits files
		if (fits_close_file(fptr, &status))
			fits_report_error(stderr, status);

		if (fits_close_file(outfptr, &status))
			fits_report_error(stderr, status);

		indice.clear();
		add_sample.clear();

	}



#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	// Close previously openened dirfile
	for (long iframe = samples_struct.iframe_min; iframe < samples_struct.iframe_max; iframe++){
		if (samples_struct.dirfile_pointers[iframe]) {
			if (gd_close(samples_struct.dirfile_pointers[iframe])){
				cerr << "EE - error closing dirfile...";
			} else {
				samples_struct.dirfile_pointers[iframe] = NULL;
			}
		}
	}

	if(rank==0)
		cout << "done." << endl << endl;

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Comm_free(&MPI_COMM_NODE);
	MPI_Comm_free(&MPI_COMM_MASTER_NODE);

	MPI_Finalize();
#endif


	if(rank==0)
		cout << endl << "End of "<< StringOf(argv[0]) << endl;

	return EXIT_SUCCESS;


}
