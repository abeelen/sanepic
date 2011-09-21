#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <sysexits.h>
#include <cmath>

#include "imageIO.h"
#include "temporary_IO.h"
#include "mpi_architecture_builder.h"
#include "parser_functions.h"
#include "estimPS.h"
#include "struct_definition.h"
#include "inputFileIO.h"
#include "crc.h"

extern "C" {
#include "wcslib/wcs.h"
#include "wcslib/wcshdr.h"
}


#ifdef PARA_FRAME
#include "mpi.h"
#include <algorithm>
#include <fstream>
#endif

using namespace std;

/*!
 *  This is organized as :
 *
 *  - parse the input ini file and verify his validity
 *  - check for existence of directory/files pointed from the ini file
 *  - Print parser output to screen
 *
 *  - for each file :
 *      - Generate or clear the dirfile parts that will be filled : fData
 *
 *  - Read all channel files, store it into a vector<vector> (and commit to other ranks if needed)
 *
 *	For each scan :
 *
 *
 */
int main(int argc, char *argv[])
{

	int size;
	int rank;

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
		cout << endl << "Beginning of sanePS : " << endl;

	struct param_saneProc proc_param; /*! A structure that contains user options about preprocessing properties */
	struct samples samples_struct; /* A structure that contains everything about frames, noise files and frame processing order */
	struct param_sanePos pos_param; /*! A structure that contains user options about map projection and properties */
	struct param_common dir; /*! structure that contains output input temp directories */

	// map making parameters
	int flagon; /*!  if one sample is rejected, flagon=1 */
	long long indpix_size; /*! indpix read size */
	//	long long indpsrc_size; /*! indpsrc read size */
	long long *indpix; // map index
	long long npix; // npix = number of filled pixels

	int nwcs=1;             // We will only deal with one wcs....
	struct wcsprm * wcs;    // wcs structure of the image
	long NAXIS1, NAXIS2;  // size of the image

	char * subheader;       // Additionnal header keywords
	int nsubkeys;           //


	string field; // actual boloname in the bolo loop
	string prefixe; // prefix used for temporary name file creation

	long iframe_min = 0, iframe_max = 0;

	// main loop variables
	double *S = NULL; // signal
	struct param_sanePS structPS;
	struct param_saneInv saneInv_struct;
	// those variables will not be used by saneProc but they are read in ini file (to check his conformity)
	//	std::vector<double> fcut_vector;
	struct param_sanePic struct_sanePic;
	string output = "";

	std::vector<std::vector<std::string> > bolo_list; // this vector contains all bolonames for all the scans

	uint16_t mask_sanePS = INI_NOT_FOUND | DATA_INPUT_PATHS_PROBLEM | OUPUT_PATH_PROBLEM | TMP_PATH_PROBLEM |
			BOLOFILE_NOT_FOUND | FSAMP_WRONG_VALUE | NCOMP_WRONG_VALUE | ELL_FILE_NOT_FOUND | MIX_FILE_NOT_FOUND |
			FITS_FILELIST_NOT_FOUND | FCUT_FILE_PROBLEM; // 0xdd1f

	//	ncomp = number of noise component to estimate
	//	fcut = cut-off freq : dont focus on freq larger than fcut for the estimation !

	// parser variables
	int indice_argv = 1;
	uint16_t parsed=0x0000; // parser error status
	uint16_t compare_to_mask; // parser error status

	//	if(rank==0){

	if ((argc < 2) || (argc > 3)) // no enough or too many arguments
		compare_to_mask = 0x0001;
	else {
		structPS.restore = 0; //default
		if (argc == 3) {

			structPS.restore = 1;
			if (strcmp(argv[1], (char*) "--restore") != 0) {
				if (strcmp(argv[2], (char*) "--restore") != 0)
					indice_argv = -1;
				else
					indice_argv = 1;
			} else {
				indice_argv = 2;
			}
		}

		if (indice_argv > 0){
			/* parse ini file and fill structures */
			parsed=parser_function(argv[indice_argv], output, dir, samples_struct, pos_param, proc_param,
					structPS, saneInv_struct, struct_sanePic, size, rank);


			compare_to_mask = parsed & mask_sanePS;

			// print parser warning and/or errors
			if (rank == 0)
				cout << endl << output << endl;
		}else
			compare_to_mask = 0x0001;

	}

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
	//	}

#ifdef PARA_FRAME

	MPI_Barrier(MPI_COMM_WORLD);

	if(configure_PARA_FRAME_samples_struct(dir.output_dir, samples_struct, rank, size, iframe_min, iframe_max)){
		MPI_Abort(MPI_COMM_WORLD, 1);
		MPI_Finalize();
		return EX_IOERR;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (iframe_max==iframe_min){ // ifram_min=iframe_max => This processor will not do anything
		cout << "Warning. Rank " << rank << " will not do anything ! please run saneFrameorder\n";
	}
#else
	iframe_min = 0;
	iframe_max = samples_struct.ntotscan;
#endif

	if(channel_list_to_vect_list(samples_struct, bolo_list, rank)){
		cout << "error in channel_list_to_vect_list" << endl;
		return EX_CONFIG;
	}

	if(rank==0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, pos_param,  proc_param,
				structPS, struct_sanePic, saneInv_struct);

		cleanup_dirfile_fdata(dir.tmp_dir, samples_struct, bolo_list);
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD); // other procs wait untill rank 0 has created dirfile architecture.
#endif


	// Open the dirfile to read temporary files
	string filedir = dir.tmp_dir + "dirfile";
	samples_struct.dirfile_pointer = gd_open((char *) filedir.c_str(), GD_RDWR | GD_VERBOSE
			| GD_UNENCODED);

	if (gd_error(samples_struct.dirfile_pointer) != 0) {
		cout << "error opening dirfile : " << filedir << endl;
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return 1;
	}

	fill_sanePS_struct(structPS, samples_struct, dir); // all ranks do it !

	long long addnpix = 0;
	int factdupl = 1;
	long long npixsrc = 0;
	long long *indpsrc;


	if(read_keyrec(dir.tmp_dir, wcs, &NAXIS1, &NAXIS2, &subheader, &nsubkeys, rank)){ // read keyrec file
#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#endif
		return (EX_IOERR);
	}


	if (pos_param.flgdupl)
		factdupl = 2; // default 0 : if flagged data are put in a duplicated map

	if (rank == 0){
		cout << "Map Size         : " << NAXIS1 << " x " << NAXIS2 << " pixels\n" << endl; // print map size

		//		if(pos_param.maskfile!=""){ // in case a mask have been used
		long long test_size;

		if (read_indpsrc(test_size, npixsrc, indpsrc, dir.tmp_dir)) { // read mask index
#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return (EX_IOERR);
		}

		if (test_size != NAXIS1 * NAXIS2) { // check size compatibility
			if (rank == 0)
				cout
				<< "indpsrc size is not the right size : Check indpsrc.bin file or run sanePos"
				<< endl;
#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return (EX_IOERR);
		}

		addnpix = samples_struct.ntotscan * npixsrc;

		if (read_indpix(indpix_size, npix, indpix, dir.tmp_dir, flagon)) { // read map indexes
#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return (EX_IOERR);
		}

		if (indpix_size != (factdupl * NAXIS1 * NAXIS2 + 2 + addnpix + 1)) { // check size compatibility
			if (rank == 0)
				cout
				<< "indpix size is not the right size : Check Indpix_*.bi file or run sanePos"
				<< endl;
#ifdef PARA_FRAME
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return (EX_IOERR);
		}

	}
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD); // other procs wait untill rank 0 has created dirfile architecture.
	MPI_Bcast(&npix,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(&npixsrc,1,MPI_LONG_LONG,0,MPI_COMM_WORLD);
	MPI_Bcast(&addnpix,1,MPI_LONG_LONG, 0, MPI_COMM_WORLD);
	MPI_Bcast(&indpix_size, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
#endif

	//First time run S=0, after sanepic, S = Pure signal
	if (structPS.signame != "") {

		S = new double[npix]; // pure signal
		fill(S,S+npix,0.0);

		// if map argument build S from map
		if (rank == 0){
#ifdef DEBUG
			cout << "Reading map.     : " << structPS.signame << endl;
#endif

			// read pure signal
			if (read_fits_signal(structPS.signame, S, indpix, NAXIS1, NAXIS2, wcs)) {
#ifdef PARA_FRAME
				MPI_Abort(MPI_COMM_WORLD, 1);
#else
				return (EX_IOERR);
#endif
			}else{

				// 	Check for unobserved Pixels
				int badPix = 0;
				for (long ii=0; ii<npix; ii++) {
					if ( isnan(S[ii]) || isinf(S[ii]) ){
						badPix++;
						S[ii] = 0.0;
					}
				}
				if (badPix > 0)
					cout << "WW - Some observed pixels fall outside the given map... set to 0" << endl;
			}

		} //rank ==0

#ifdef PARA_FRAME

		MPI_Barrier(MPI_COMM_WORLD);

		if(rank!=0){
			indpix = new long long[indpix_size];
			S = new double[npix];
			//			fill(S,S+npix,0.0);
		}

		MPI_Bcast(indpix,indpix_size,MPI_LONG_LONG,0,MPI_COMM_WORLD);
		MPI_Bcast(S,npix,MPI_DOUBLE,0,MPI_COMM_WORLD); // broadcast it to the other procs

#endif
		wcsvfree(&nwcs, &wcs);

	} // structPS.signame != ""

	if (structPS.restore) { // restore incomplete work with previous saved data
		if (rank == 0){
			cout << "Checking previous session\n";
			struct checksum chk_t, chk_t2;
			compute_checksum(dir, pos_param, proc_param, saneInv_struct, structPS, struct_sanePic, samples_struct, npix,
					indpix, indpsrc, NAXIS1 * NAXIS2, chk_t);
			read_checksum(dir.tmp_dir, chk_t2, "sanePS"); // read previous checksum
			if (compare_checksum(chk_t, chk_t2)) { // compare them
				cout << "Checksums are different !!! Exiting..." << endl;
#ifdef USE_MPI
				MPI_Abort(MPI_COMM_WORLD, 1);
#endif
				return EX_CONFIG;
			}
		}
	}

	if (structPS.save_data) {
		if (rank == 0) {
			struct checksum chk_t;
			/* Compute Checsum for crash recovery ! */
			compute_checksum(dir, pos_param, proc_param, saneInv_struct, structPS, struct_sanePic, samples_struct, npix,
					indpix, indpsrc, NAXIS1 * NAXIS2, chk_t);
			if(write_checksum(dir.tmp_dir, chk_t, "sanePS")){ // write down on disk the checksum values
#ifdef USE_MPI
				MPI_Abort(MPI_COMM_WORLD, 1);
#endif
				return EX_CANTCREAT;
			}
		}
	}

#ifdef PARA_FRAME
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(rank==0)
		cout << "Noise Power Spectra Estimation started : " << endl;

	for (long iframe = iframe_min; iframe < iframe_max; iframe++) { // proceed scan by scan

		std::vector<string> det=bolo_list[iframe];

		//		string output_read = "";
		//		if(read_channel_list(output_read, samples_struct.bolovect[iframe], det)){
		//			cout << output_read << endl;
		//			cerr << "input channel list could not be read for : " << samples_struct.fitsvect[iframe] << " . Skipping...\n";
		//			continue; // skip the file if channel list (was not found) / (is incorrect) !
		//		}

		if(EstimPowerSpectra(det, proc_param, dir, pos_param, structPS,
				samples_struct, NAXIS1, NAXIS2, npix, iframe, indpix, S, rank)){
			cout << "Error in EstimPowerSpectra procedure. Exiting ...\n";
			break;
		}
		// ns = number of samples in the "iframe" scan
		// npix = total number of filled pixels
		// iframe = scan number
		// indpix = pixels index
		// MixMatfile = this file contains the number of components in the common-mode component of the noise
		// and the value of alpha, the amplitude factor which depends on detectors but not on time (see formulae (3) in "Sanepic:[...], Patanchon et al.")
	}

	if (rank == 0)
		cout << endl << "done." << endl;

	fftw_cleanup();

	if (structPS.signame != "") {
		delete[] S;
		delete[] indpix;
	}else
		if(rank==0)
			delete[] indpix;


	if(rank==0)
		delete[] indpsrc;

	if (gd_close(samples_struct.dirfile_pointer))
		cout << "error closing dirfile : " << filedir << endl;

#ifdef PARA_FRAME

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif

	if(rank==0)
		cout << endl << "END OF SANEPS" << endl;

	return EXIT_SUCCESS;
}
