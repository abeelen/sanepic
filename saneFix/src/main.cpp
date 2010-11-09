#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parse_saneFix.h"
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

#if defined(PARA_FRAME)
#define USE_MPI
#endif

#ifdef USE_MPI
#include "mpi.h"
#endif



int main(int argc, char *argv[]) {


	int rank, size; /* MPI processor rank and MPI total number of used processors */

#ifdef USE_MPI

	// setup MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

#else
	size = 1;
	rank = 0;
	cout << "Mpi is not used for this step" << endl;
#endif

	int parsed=1;

	struct samples samples_struct; /* fits file list + number of scans */
	struct param_common dir; /* directories : temporary input and output */
	std::vector<long> indice; /*! gap indexes (sample index) */
	std::vector <long> add_sample; /*! number of samples to add per gap */
	double fsamp; /*! sampling frequency */
	string output ="";

	if(rank==0)
		printf("\nBeginning of saneFix:\n\n");

	// get directories and fits file list
	parsed=parse_saneFix_ini_file(argv[1], output, dir,
			samples_struct, rank);

	// print parser warning and/or errors
	cout << endl << output << endl;

	if(parsed==-1){ /* error during parsing phase */
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return EX_CONFIG;
	}

	for(long ii=0; ii<samples_struct.ntotscan;ii++){ // for each input scan

		int do_it=who_do_it(size, rank, ii); /* which rank do the job ? */

		if(rank==do_it){ /* if this rank has to do the job ... */

			int format_fits; // 1 = HIPE, 2 = Sanepic
			std::vector <long> suppress_time_sample;

			cout  << "[ " << rank << " ] " << "Fixing file : " << samples_struct.fitsvect[ii] << endl;

			// fits files pointer
			fitsfile * fptr;
			fitsfile *outfptr;

			// lookfor saneCheck log files
			if((read_indices_file(samples_struct.fitsvect[ii],dir,  indice, fsamp))){
				cout << "[ " << rank << " ] " << "Skipping file : " << samples_struct.fitsvect[ii] << ". Please run saneCheck on this file before\n";
				continue; // log file were not found for the 'ii'th fits file
			}

			format_fits=test_format(samples_struct.fitsvect[ii]); // get fits file format
			if(format_fits==0){
				cerr << "[ " << rank << " ] " << "input fits file format is undefined : " << samples_struct.fitsvect[ii] << " . Skipping file...\n";
				continue;
			}


			// compute the number of sample that must be added to fill the gaps and have a continous timeline
			long samples_to_add=how_many(samples_struct.fitsvect[ii], samples_struct.nsamples[ii] ,indice, fsamp, add_sample,suppress_time_sample);
			cout << "[ " << rank << " ] " << "samptoadd : " << samples_to_add << endl;
			cout << "[ " << rank << " ] " << "total : " << samples_struct.nsamples[ii] + samples_to_add << endl;

			// total number of samples in the fixed fits file
			long ns_total = samples_struct.nsamples[ii] + samples_to_add;

			string fname2 = dir.output_dir + FitsBasename(samples_struct.fitsvect[ii]) + "_fixed.fits"; // output fits filename
			string fname=samples_struct.fitsvect[ii]; // input fits filename

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
			read_bolo_list(samples_struct.fitsvect[ii], det, ndet);

			// Copy primary Header
			fits_copy_header(fptr, outfptr, &status);

			if(format_fits==1){ // HIPE format
				cout << "[ " << rank << " ] " << "HIPE format found\n";

				// 1 signal
				fix_signal(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, ndet, indice, add_sample, suppress_time_sample);

				// 2 RA 3 DEC
				fix_RA_DEC(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, ndet, indice, add_sample, suppress_time_sample);

				// 4 mask
				fix_mask(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, ndet, indice, add_sample, suppress_time_sample);

				// 5 time
				fix_time_table(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, indice, add_sample, samples_struct.nsamples[ii], fsamp,suppress_time_sample);

				// 6 channels
				copy_channels(fptr, outfptr);

				// 7 reference positions
				fix_ref_pos(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, indice, add_sample, suppress_time_sample);

				// 8 offsets
				copy_offsets(fptr, outfptr);

			}else{ // format =2
				cout << "[ " << rank << " ] " << "SANEPIC format found\n";

				// 1 ref pos
				fix_ref_pos(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, indice, add_sample, suppress_time_sample);

				// 2 offsets
				copy_offsets(fptr, outfptr);

				// 3 channels
				copy_channels(fptr, outfptr);

				// 4 time
				fix_time_table(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, indice, add_sample, samples_struct.nsamples[ii], fsamp,suppress_time_sample);

				// 5 signal
				fix_signal(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, ndet, indice, add_sample, suppress_time_sample);

				// 6 mask
				fix_mask(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, ndet, indice, add_sample, suppress_time_sample);

			}

			// close both fits files
			if (fits_close_file(fptr, &status))
				fits_report_error(stderr, status);

			if (fits_close_file(outfptr, &status))
				fits_report_error(stderr, status);

			indice.clear();
			add_sample.clear();

		}

	}

	//clean up
	delete [] samples_struct.nsamples;

	if(rank==0)
		cout << "\nEnd of saneFix\n";

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif


	return EXIT_SUCCESS;


}
