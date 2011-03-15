#include "inputFileIO.h"
#include "mpi_architecture_builder.h"
#include "dataIO.h"
#include "parser_functions.h"
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

#ifdef PARA_FRAME
#include "mpi.h"
#endif


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
	cout << "Mpi is not used for this step" << endl;
#endif

	struct samples samples_struct; /* fits file list + number of scans */
	struct param_common dir; /* directories : temporary input and output */
	struct param_sanePos pos_param;
	struct param_sanePre proc_param;
	struct param_sanePic sanePic_struct;
	struct param_saneInv saneInv_struct;
	struct param_sanePS structPS;
	struct saneCheck check_struct;
	std::vector<long> indice; /*! gap indexes (sample index) */
	std::vector <long> add_sample; /*! number of samples to add per gap */
	double fsamp; /*! sampling frequency */
	string output ="";

	uint16_t mask_sanefix = INI_NOT_FOUND | DATA_INPUT_PATHS_PROBLEM | OUPUT_PATH_PROBLEM | TMP_PATH_PROBLEM |
				FSAMP_WRONG_VALUE | FITS_FILELIST_NOT_FOUND; // 0x410f


	if(rank==0)
		printf("\nBeginning of saneFix:\n\n");



	if (rank==0){ // root parse ini file and fill the structures. Also print warnings or errors

		uint16_t parsed=0x0000; // parser error status
		uint16_t compare_to_mask; // parser error status


		if (argc<2)/* not enough argument */
			compare_to_mask=0x001;
		else {
			/* parse ini file and fill structures */
			parsed=parser_function(argv[1], output, dir, samples_struct, pos_param, proc_param,
					structPS, saneInv_struct, sanePic_struct, size, rank);

			compare_to_mask = parsed & mask_sanefix;
			//			cout << "compare to mask : " << (int)compare_to_mask << endl;

			//	// get directories and fits file list
			//	parsed=parse_saneFix_ini_file(argv[1], output, dir,
			//			samples_struct, rank);

		}

		// print parser warning and/or errors
		cout << endl << output << endl;

		// in case there is a parsing error or the dirfile format file was not created correctly
		if(compare_to_mask>0x0000){

			switch (compare_to_mask){/* error during parsing phase */

			case 0x0001: printf("Please run %s using a correct *.ini file\n",argv[0]);
			break;

			default : printf("Wrong program options or argument. Exiting !\n");
			break;


			}

#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return EX_CONFIG;
		}
	}


#ifdef PARA_FRAME

	MPI_Datatype message_type;
	struct ini_var_strings ini_v;
	int ntotscan;

	if(rank==0){
		fill_var_sizes_struct(dir, pos_param, proc_param,
				saneInv_struct, structPS, samples_struct, ini_v);

		ntotscan = ini_v.ntotscan;
	}

	MPI_Bcast(&ntotscan, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(rank!=0){
		ini_v.fitsvect=new int[ntotscan];
		ini_v.noisevect=new int[ntotscan];
		ini_v.bolovect=new int[ntotscan];
	}

	ini_v.ntotscan=ntotscan;

	Build_derived_type_ini_var (&ini_v,	&message_type);

	MPI_Bcast(&ini_v, 1, message_type, 0, MPI_COMM_WORLD);

	commit_struct_from_root(dir, pos_param, proc_param, saneInv_struct, sanePic_struct, structPS, samples_struct, ini_v, rank);

	MPI_Barrier(MPI_COMM_WORLD);

#endif

	if(rank==0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, pos_param,  proc_param,
				structPS, sanePic_struct, saneInv_struct);
	}

	for(long ii=0; ii<samples_struct.ntotscan;ii++){ // for each input scan

		int do_it=who_do_it(size, rank, ii); /* which rank do the job ? */

		if(rank==do_it){ /* if this rank has to do the job ... */

			int format_fits; // 1 = HIPE, 2 = Sanepic
			std::vector <long> suppress_time_sample;
			long init_num_delete =0;
			long end_num_delete =0;

			cout  << "[ " << rank << " ] " << "Fixing file : " << samples_struct.fitsvect[ii] << endl;

			// fits files pointer
			fitsfile * fptr;
			fitsfile *outfptr;

			// lookfor saneCheck log files
			if((read_indices_file(samples_struct.fitsvect[ii],dir,  indice, fsamp, init_num_delete, end_num_delete))){
				cout << "[ " << rank << " ] " << "Skipping file : " << samples_struct.fitsvect[ii] << ". Please run saneCheck on this file before\n";
				continue; // log file were not found for the 'ii'th fits file
			}

			format_fits=test_format(samples_struct.fitsvect[ii]); // get fits file format
			if(format_fits==0){
				cerr << "[ " << rank << " ] " << "input fits file format is undefined : " << samples_struct.fitsvect[ii] << " . Skipping file...\n";
				continue;
			}


			//			cout << "readed : " << fsamp << " " << init_num_delete  << " " << end_num_delete << " " << indice[0] << " " << indice[1] << endl;

			refresh_indice(fsamp, init_num_delete, end_num_delete, indice, samples_struct.nsamples[ii]);

			//			cout << "refresh : " << fsamp << " " << init_num_delete  << " " << end_num_delete << " " << indice.size() << endl;

			// compute the number of sample that must be added to fill the gaps and have a continous timeline
			long samples_to_add=how_many(samples_struct.fitsvect[ii], samples_struct.nsamples[ii] ,indice, fsamp, add_sample,suppress_time_sample);
			cout << "[ " << rank << " ] " << "samptoadd : " << samples_to_add << endl;
			cout << "[ " << rank << " ] " << "samples to delete (begin/end) : (" << init_num_delete << "/" <<  end_num_delete << ")" << endl;
			cout << "[ " << rank << " ] " << "total : " << samples_struct.nsamples[ii] + samples_to_add - init_num_delete - end_num_delete << endl;

			// total number of samples in the fixed fits file
			long ns_total = samples_struct.nsamples[ii] + samples_to_add - init_num_delete - end_num_delete;

			if(samples_to_add+init_num_delete+end_num_delete==0){
				cout << "[ " << rank << " ] " << "Nothing to do for : " << samples_struct.fitsvect[ii] << " . Skipping file...\n";
				continue;
			}

			string fname2 = "!" + dir.output_dir + FitsBasename(samples_struct.fitsvect[ii]) + "_fixed.fits"; // output fits filename
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
				fix_signal(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, ndet, indice, add_sample, suppress_time_sample, init_num_delete);

				// 2 RA 3 DEC
				fix_RA_DEC(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, ndet, indice, add_sample, suppress_time_sample, init_num_delete);

				// 4 mask
				fix_mask(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, ndet, indice, add_sample, suppress_time_sample, init_num_delete);

				// 5 time
				fix_time_table(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, indice, add_sample, samples_struct.nsamples[ii], fsamp,suppress_time_sample, init_num_delete);

				// 6 channels
				copy_channels(fptr, outfptr);

				// 7 reference positions
				fix_ref_pos(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, indice, add_sample, suppress_time_sample, init_num_delete);

				// 8 offsets
				copy_offsets(fptr, outfptr);

			}else{ // format =2
				cout << "[ " << rank << " ] " << "SANEPIC format found\n";

				// 1 ref pos
				fix_ref_pos(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, indice, add_sample, suppress_time_sample, init_num_delete);

				// 2 offsets
				copy_offsets(fptr, outfptr);

				// 3 channels
				copy_channels(fptr, outfptr);

				// 4 time
				fix_time_table(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, indice, add_sample, samples_struct.nsamples[ii], fsamp, suppress_time_sample, init_num_delete);

				// 5 signal
				fix_signal(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, ndet, indice, add_sample, suppress_time_sample, init_num_delete);

				// 6 mask
				fix_mask(fptr, outfptr, samples_struct.fitsvect[ii], ns_total, det, ndet, indice, add_sample, suppress_time_sample, init_num_delete);

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
	//	delete [] samples_struct.nsamples;

	if(rank==0)
		cout << "\nEnd of saneFix\n";

#ifdef PARA_FRAME
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif


	return EXIT_SUCCESS;


}
