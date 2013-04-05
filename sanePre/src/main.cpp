#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <time.h>
#include <cmath>
#include <vector>
#include <string>
#include <fftw3.h>
#include <sysexits.h>


#include "ImageIO.h"
#include "TemporaryIO.h"
#include "MPIConfiguration.h"
#include "ParserFunctions.h"
#include "StructDefinition.h"
#include "InputFileIO.h"
#include "ErrorCode.h"


extern "C" {
#include "nrutil.h"
#include "wcslib/wcs.h"     // to export wcs header keys
#include "wcslib/wcshdr.h"
#include <fftw3.h>          // for wisdom
}


#ifdef USE_MPI
#include "mpi.h"
#include <algorithm>
#include <fstream>
#endif

using namespace std;


//**********************************************************************************//
//**********************************************************************************//
//*************************** Beginning of main program ****************************//
//**********************************************************************************//
//**********************************************************************************//


/*! \mainpage SanePre for SPIRE
 *
 * \section intro_sec Matthieu HUSSON & Alexandre Beelen
 *
 * Sanepic for SPIRE Manual
 */

/*!
 *  This is organized as :
 *
 *  - parse the input ini file and verify his validity
 *  - check for existence of directory/files pointed from the ini file
 *  - Print parser output to screen
 *
 *  - Generate, or clear dirfile tree (folders and format files)
 *  - for each file :
 *  	- Generate or clear the dirfile parts that will be filled : data, flag, LON, LAT
 *
 *  - Read all channel files, store it into a vector<vector> (and commit to other ranks if needed)
 *
 *  - Copy from scans (fits) to disk (getdata binaries in a dirfile tree) : data, flag, LON, LAT
 *
 */

int main(int argc, char *argv[]){


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
		cout << endl << "sanePre :  data distribution" << endl;

	struct param_saneProc Proc_param; /* contains user options about preprocessing properties */
	struct samples samples_struct;  /*  everything about frames, noise files and frame processing order */
	struct param_sanePos Pos_param; /* contains user options about map projection and properties */
	struct param_common dir; /* contains output input temp directories */

	// those variables will not be used by sanePre but they are read in ini file (to check his conformity)
	struct param_sanePS   PS_param;
	struct param_saneInv Inv_param;
	struct param_sanePic Pic_param;
	string parser_output = "";

	//	uint32_t mask_sanePre = 0x405f;
	uint32_t mask_sanePre = INI_NOT_FOUND | DATA_INPUT_PATHS_PROBLEM | OUTPUT_PATH_PROBLEM | TMP_PATH_PROBLEM |
			BOLOFILE_NOT_FOUND | FILEFORMAT_NOT_FOUND | FITS_FILELIST_NOT_FOUND; // 0x405f

#ifdef DEBUG
	// Processing time estimation
	time_t t2, t3;
#endif

	uint32_t parsed=0x0000; // parser error status
	uint32_t compare_to_mask; // parser error status

	if (argc<2) {/* not enough argument */
		if (rank == 0)
			cerr << "EE - Please run  " << argv[0] << " with a .ini file" << endl;
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD );
		MPI_Finalize();
#endif
		exit(EXIT_FAILURE);
	}  else {

		/* parse ini file and fill structures */
		parsed=parser_function(argv[1], parser_output, dir, samples_struct, Pos_param, Proc_param, PS_param, Inv_param, Pic_param, size, rank);

		compare_to_mask = parsed & mask_sanePre;

		// print parser warning and/or errors
		if (rank==0)
			cout << endl << parser_output << endl;

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

#ifdef DEBUG
	// processing begins here
	t2=time(NULL);
#endif

	// Start...

	if(rank==0){
		// parser print screen function
		parser_printOut(argv[0], dir, samples_struct, Pos_param,  Proc_param,
				PS_param, Pic_param, Inv_param);
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
	parser_output.clear();

	if (rank==0){
		cout << endl << "Initialization... " << endl;
	}
	// Create tmpDir if needed (bolo_rank_dummy here to prevent concurrence in the frame case on same node
	if (init_tmpdir(parser_output, samples_struct, dir.tmp_dir, node_rank) ){
		cerr << parser_output;
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return EX_CONFIG;
	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_NODE);
#endif

	// Create empty dirfile structure
	if ( init_dirfile(dir.tmp_dir, samples_struct, bolo_rank)) {
		cerr << "EE - Error in initializing dirfile" << endl;
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
#endif
		return EX_CONFIG;
	}


	// Read file size once for all
	if (rank == 0)
		readFramesFromFits(samples_struct);
#ifdef USE_MPI
	MPI_Bcast_vector_long(samples_struct.nsamples, 0, MPI_COMM_WORLD);
#endif

	if(rank==0)
		cout << "Exporting signal and flags... " << endl;

	if (bolo_rank == 0) {
		if(writeDataFlagToDirfile(dir, samples_struct)){
			cerr << "EE - write_data_flag_to_dirfile !! Exiting ..." << endl;
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return EXIT_FAILURE;
		}
	}


#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_NODE);
#endif


	if(rank==0)
		cout << "Exporting positions... " << endl;

	if (bolo_rank == 0 ) {
		switch (Pos_param.fileFormat) {
		case 0:
			if(exportLonLatToDirfile(dir, samples_struct)){
				cerr << "EE - exportLON_LAT_to_dirfile !! Exiting ..." << endl;
#ifdef USE_MPI
				MPI_Abort(MPI_COMM_WORLD, 1);
#endif
				return EXIT_FAILURE;
			}
			break;
		case 1:
			if(writeLonLatToDirfile(dir, samples_struct)){
				cerr << "EE - write_LON_LAT_to_dirfile !! Exiting ..." << endl;
#ifdef USE_MPI
				MPI_Abort(MPI_COMM_WORLD, 1);
#endif
				return EXIT_FAILURE;
			}
			break;
		}
	}


	// Create a fake WCS image header and populate it with info from the first fits file
	if (rank == 0){

		cout << "Exporting header keys... " << endl;

		int nwcs=1;                      // We will only deal with one wcs....
		struct wcsprm * wcs;             // wcs structure of the image

		char * subheader;               // Additionnal header keywords
		int nsubkeys;                   //


		if (get_fits_META(samples_struct.fitsvect[0], wcs, &subheader, &nsubkeys))
			cerr << "WW - Problem getting fits META" << endl;

		// Updated wcs with value from the ini file if any.
		if (Pos_param.equinox != 0.0)
			wcs->equinox = Pos_param.equinox;
		if (Pos_param.restwav != 0.0)
			wcs->restwav = Pos_param.restwav;
		if (Pos_param.radesys != "")
			strcpy(wcs->radesys, Pos_param.radesys.c_str());


		if(save_keyrec(dir.tmp_dir,wcs, 0, 0, subheader, nsubkeys)){
#ifdef USE_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
#endif
			return(EX_CANTCREAT);
		}

		wcsvfree(&nwcs, &wcs);

	}

	if (rank == 0 && Proc_param.wisdom)
		cout << "Building some wisdom..." << endl;

	if (bolo_rank == 0 && Proc_param.wisdom ){

		fftw_complex * fdata;
		double *        data;
		long ns;
		fftw_plan plan;

		for (long iframe=samples_struct.iframe_min ;iframe<samples_struct.iframe_max;iframe++){

			ns = samples_struct.nsamples[iframe];

			data  = (double *) fftw_malloc(sizeof(double)*ns);
			fdata = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (ns/2+1));

			plan = fftw_plan_dft_r2c_1d(ns, data, fdata, FFTW_PATIENT);
			//			fftw_print_plan(plan);
			fftw_destroy_plan(plan);

			plan = fftw_plan_dft_c2r_1d(ns, fdata, data, FFTW_PATIENT);
			//			fftw_print_plan(plan);
			fftw_destroy_plan(plan);

			fftw_free(data);
			fftw_free(fdata);
		}

		string filename = dir.tmp_dir + "fftw.wisdom";
		FILE * pFilename;

		pFilename = fopen((const char*) filename.c_str(), "w+");

		if ( pFilename != NULL ){
			fftw_export_wisdom_to_file( pFilename );
		} else {
			cerr << "EE - Problem while writing FFTW wisdom in " << filename << endl;
		}

		fclose(pFilename);

	}

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_NODE);
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

#ifdef DEBUG
	//Get processing time
	t3=time(NULL);

	if(rank==0)
		cout << "Total Time : " << t3-t2 << " sec\n\n";

#endif

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
